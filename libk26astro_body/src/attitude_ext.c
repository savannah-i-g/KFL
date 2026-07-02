/* attitude_ext.c — full 3×3-inertia attitude state + integration.
 *
 * For spacecraft whose mass distribution carries off-diagonal
 * inertia coupling and whose inertia changes mid-simulation (stage
 * events, propellant depletion, deployable extension). Sibling of
 * the diagonal-inertia path in attitude.c.
 *
 * Equations of motion (body frame):
 *
 *   ω̇ = I⁻¹ (τ − ω × (I ω))
 *   q̇ = (1/2) q ⊗ ω_quat
 *
 * Implementation: I⁻¹ is precomputed at init and recomputed by
 * k26astro_attitude_update_inertia. The 3×3 inverse is inlined here
 * (cofactor / adjugate form); K26M3 algebra has no consumers in
 * libk26m3d yet so the matrix utility stays local until enough
 * call sites accrue to motivate promotion.
 *
 * References:
 *   - Markley & Crassidis, Fundamentals of Spacecraft Attitude
 *     Determination and Control (2014), §3.7 — quaternion
 *     exponential map; §3.6.2 — Euler's rotational equation in
 *     non-principal axes.
 *   - Hughes, Spacecraft Attitude Dynamics (1986), §4.5 — full
 *     inertia tensor form of Euler's equation.
 *   - Wertz (ed.), Spacecraft Attitude Determination and Control
 *     (1978), §16.3 — inertia tensor properties (symmetric
 *     positive-definite). */
#include "k26astro_body/attitude.h"

#include <math.h>
#include <stdlib.h>

/* ---- Local 3×3 algebra helpers ------------------------------- */

static K26V3 m3_mul_v3_(const K26M3 *m, K26V3 v)
{
    K26V3 out;
    out.x = m->m[0][0] * v.x + m->m[0][1] * v.y + m->m[0][2] * v.z;
    out.y = m->m[1][0] * v.x + m->m[1][1] * v.y + m->m[1][2] * v.z;
    out.z = m->m[2][0] * v.x + m->m[2][1] * v.y + m->m[2][2] * v.z;
    return out;
}

/* Cofactor inverse for a 3×3 row-major matrix. Returns the inverse
 * via `out` and 1 on success; on singular input (|det| < 1e-30)
 * writes the zero matrix and returns 0. */
static int m3_inverse_(const K26M3 *in, K26M3 *out)
{
    double a = in->m[0][0], b = in->m[0][1], c = in->m[0][2];
    double d = in->m[1][0], e = in->m[1][1], f = in->m[1][2];
    double g = in->m[2][0], h = in->m[2][1], i = in->m[2][2];

    double cof00 = e * i - f * h;
    double cof01 = c * h - b * i;
    double cof02 = b * f - c * e;
    double cof10 = f * g - d * i;
    double cof11 = a * i - c * g;
    double cof12 = c * d - a * f;
    double cof20 = d * h - e * g;
    double cof21 = b * g - a * h;
    double cof22 = a * e - b * d;

    double det = a * cof00 + b * cof10 + c * cof20;
    if (fabs(det) < 1.0e-30) {
        for (int r = 0; r < 3; r++)
            for (int cc = 0; cc < 3; cc++)
                out->m[r][cc] = 0.0;
        return 0;
    }
    double inv = 1.0 / det;
    out->m[0][0] = cof00 * inv;
    out->m[0][1] = cof01 * inv;
    out->m[0][2] = cof02 * inv;
    out->m[1][0] = cof10 * inv;
    out->m[1][1] = cof11 * inv;
    out->m[1][2] = cof12 * inv;
    out->m[2][0] = cof20 * inv;
    out->m[2][1] = cof21 * inv;
    out->m[2][2] = cof22 * inv;
    return 1;
}

/* ---- Ext torque-source registry ------------------------------ */

struct K26AstroTorqueListExt {
    K26AstroTorqueFnExt *fns;
    void               **ctxs;
    int                  count;
    int                  capacity;
};

/* The K26AstroTorqueList forward-declaration in attitude.h is
 * opaque, but both the diagonal-path and Ext-path registries use
 * pointer-sized slots; the Ext path stores its own list type via
 * the same `torques` pointer slot. A cast at clear/dispatch sites
 * disambiguates. The two types never alias because each state struct
 * is created and torn down by its own init/destroy pair. */

int k26astro_attitude_register_torque_ext(K26AstroAttitudeStateExt *a,
                                          K26AstroTorqueFnExt fn,
                                          void *ctx)
{
    if (!a || !fn) return 1;
    if (!a->torques) {
        a->torques = (K26AstroTorqueList *)
            calloc(1, sizeof(struct K26AstroTorqueListExt));
        if (!a->torques) return 2;
    }
    struct K26AstroTorqueListExt *p = (struct K26AstroTorqueListExt *)a->torques;
    if (p->count >= p->capacity) {
        int new_cap = p->capacity ? p->capacity * 2 : 4;
        K26AstroTorqueFnExt *new_fns = realloc(p->fns,
            (size_t)new_cap * sizeof(K26AstroTorqueFnExt));
        void **new_ctxs = realloc(p->ctxs,
            (size_t)new_cap * sizeof(void *));
        if (!new_fns || !new_ctxs) {
            free(new_fns);
            free(new_ctxs);
            return 2;
        }
        p->fns  = new_fns;
        p->ctxs = new_ctxs;
        p->capacity = new_cap;
    }
    p->fns [p->count] = fn;
    p->ctxs[p->count] = ctx;
    p->count++;
    return 0;
}

void k26astro_attitude_clear_torques_ext(K26AstroAttitudeStateExt *a)
{
    if (!a || !a->torques) return;
    struct K26AstroTorqueListExt *p = (struct K26AstroTorqueListExt *)a->torques;
    free(p->fns);
    free(p->ctxs);
    free(p);
    a->torques = NULL;
}

/* ---- Lifecycle ------------------------------------------------ */

void k26astro_attitude_init_ext(K26AstroAttitudeStateExt *a,
                                K26M3 inertia)
{
    if (!a) return;
    a->q = k26m3d_quat_identity();
    a->omega_body.x = a->omega_body.y = a->omega_body.z = 0.0;
    a->inertia = inertia;
    (void)m3_inverse_(&a->inertia, &a->inertia_inverse);
    a->torques = NULL;
}

void k26astro_attitude_destroy_ext(K26AstroAttitudeStateExt *a)
{
    k26astro_attitude_clear_torques_ext(a);
}

void k26astro_attitude_update_inertia(K26AstroAttitudeStateExt *a,
                                      K26M3 new_inertia)
{
    if (!a) return;
    a->inertia = new_inertia;
    (void)m3_inverse_(&a->inertia, &a->inertia_inverse);
}

/* ---- Integration --------------------------------------------- */

void k26astro_attitude_step_free_ext(K26AstroAttitudeStateExt *a, double dt)
{
    if (!a) return;
    /* Free-body Euler: ω̇ = I⁻¹ (−ω × Iω). */
    K26V3 Iw   = m3_mul_v3_(&a->inertia, a->omega_body);
    K26V3 cross = k26m3d_v3_cross(a->omega_body, Iw);
    K26V3 wdot = m3_mul_v3_(&a->inertia_inverse,
                            (K26V3){ -cross.x, -cross.y, -cross.z });
    a->omega_body.x += wdot.x * dt;
    a->omega_body.y += wdot.y * dt;
    a->omega_body.z += wdot.z * dt;
    K26V3 theta = {
        a->omega_body.x * dt,
        a->omega_body.y * dt,
        a->omega_body.z * dt
    };
    K26Quat dq = k26astro_quat_exp_half(theta);
    a->q = k26m3d_quat_norm(k26m3d_quat_mul(a->q, dq));
}

void k26astro_attitude_step_torque_ext(K26AstroAttitudeStateExt *a,
                                       K26V3 torque_body, double dt)
{
    if (!a) return;
    /* ω̇ = I⁻¹ (τ − ω × (I ω)). */
    K26V3 Iw    = m3_mul_v3_(&a->inertia, a->omega_body);
    K26V3 cross = k26m3d_v3_cross(a->omega_body, Iw);
    K26V3 rhs = {
        torque_body.x - cross.x,
        torque_body.y - cross.y,
        torque_body.z - cross.z
    };
    K26V3 wdot = m3_mul_v3_(&a->inertia_inverse, rhs);
    a->omega_body.x += wdot.x * dt;
    a->omega_body.y += wdot.y * dt;
    a->omega_body.z += wdot.z * dt;
    K26V3 theta = {
        a->omega_body.x * dt,
        a->omega_body.y * dt,
        a->omega_body.z * dt
    };
    K26Quat dq = k26astro_quat_exp_half(theta);
    a->q = k26m3d_quat_norm(k26m3d_quat_mul(a->q, dq));
}

void k26astro_attitude_step_torque_registry_ext(K26AstroAttitudeStateExt *a,
                                                double t, double dt)
{
    if (!a) return;
    K26V3 sum = { 0.0, 0.0, 0.0 };
    if (a->torques) {
        struct K26AstroTorqueListExt *p =
            (struct K26AstroTorqueListExt *)a->torques;
        for (int i = 0; i < p->count; i++) {
            p->fns[i](a, t, &sum, p->ctxs[i]);
        }
    }
    k26astro_attitude_step_torque_ext(a, sum, dt);
}
