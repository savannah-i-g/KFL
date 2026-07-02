/* test_attitude_step_ext.c — full-3×3-inertia attitude integration
 * gate.
 *
 * Validates:
 *
 *   (1) k26astro_attitude_init_ext computes the inertia inverse;
 *       I · I⁻¹ ≈ I_3 within FP tolerance.
 *   (2) Free-body spin around a principal axis is stationary —
 *       ω̇ = I⁻¹(−ω × Iω) vanishes when ω is collinear with an
 *       eigenvector of I.
 *   (3) Free-body spin off a principal axis exhibits Eulerian
 *       nutation; for a symmetric top the polhode period
 *       T = 2π · I_t / (ωz · |I_s − I_t|) is exact and ω returns
 *       to its starting value within accumulated integration error.
 *   (4) A constant torque about a principal axis with ω₀ = 0 gives
 *       the closed-form Δω = τ / I · dt (cross-coupling term zero).
 *   (5) Inertia update via k26astro_attitude_update_inertia
 *       recomputes the inverse correctly.
 *   (6) Torque registry on the Ext path dispatches independently of
 *       the diagonal path's registry. */
#include "k26astro_body/attitude.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

static K26M3 diag_(double xx, double yy, double zz)
{
    K26M3 m;
    m.m[0][0] = xx; m.m[0][1] = 0;  m.m[0][2] = 0;
    m.m[1][0] = 0;  m.m[1][1] = yy; m.m[1][2] = 0;
    m.m[2][0] = 0;  m.m[2][1] = 0;  m.m[2][2] = zz;
    return m;
}

/* I · I⁻¹ entry (r,c). */
static double mat_mul_at_(const K26M3 *a, const K26M3 *b, int r, int c)
{
    return a->m[r][0] * b->m[0][c]
         + a->m[r][1] * b->m[1][c]
         + a->m[r][2] * b->m[2][c];
}

/* Constant body-frame torque source for the registry test. */
static void torque_const_x_(const K26AstroAttitudeStateExt *state, double t,
                            K26V3 *out, void *ctx)
{
    (void)state; (void)t; (void)ctx;
    out->x += 1.5;
}

int main(void)
{
    /* ---- (1) Init computes inverse ----------------------- */
    K26AstroAttitudeStateExt a;
    K26M3 I = diag_(2.0, 3.0, 5.0);
    k26astro_attitude_init_ext(&a, I);
    /* Diagonal inverse: diag(1/2, 1/3, 1/5). */
    assert(near_(a.inertia_inverse.m[0][0], 0.5,       1e-12));
    assert(near_(a.inertia_inverse.m[1][1], 1.0 / 3.0, 1e-12));
    assert(near_(a.inertia_inverse.m[2][2], 0.2,       1e-12));
    /* I · I⁻¹ = I_3. */
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            double e = (r == c) ? 1.0 : 0.0;
            assert(near_(mat_mul_at_(&a.inertia, &a.inertia_inverse, r, c),
                         e, 1e-12));
        }
    }

    /* ---- Non-diagonal inverse: I · I⁻¹ = I_3 ------------- */
    K26M3 J;
    J.m[0][0] = 4.0; J.m[0][1] = 1.0; J.m[0][2] = 0.5;
    J.m[1][0] = 1.0; J.m[1][1] = 5.0; J.m[1][2] = 0.2;
    J.m[2][0] = 0.5; J.m[2][1] = 0.2; J.m[2][2] = 3.0;
    K26AstroAttitudeStateExt nd;
    k26astro_attitude_init_ext(&nd, J);
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            double e = (r == c) ? 1.0 : 0.0;
            assert(fabs(mat_mul_at_(&nd.inertia, &nd.inertia_inverse, r, c) - e)
                   < 1e-12);
        }
    }
    k26astro_attitude_destroy_ext(&nd);

    /* ---- (2) Principal-axis spin is stationary ----------- */
    K26AstroAttitudeStateExt p;
    k26astro_attitude_init_ext(&p, diag_(1.0, 2.0, 3.0));
    p.omega_body.x = 1.0;
    p.omega_body.y = 0.0;
    p.omega_body.z = 0.0;
    for (int n = 0; n < 1000; n++) {
        k26astro_attitude_step_free_ext(&p, 0.001);
    }
    /* ω should still be (1, 0, 0) ± numerical noise. */
    assert(near_(p.omega_body.x, 1.0, 1e-10));
    assert(fabs(p.omega_body.y) < 1e-10);
    assert(fabs(p.omega_body.z) < 1e-10);
    k26astro_attitude_destroy_ext(&p);

    /* ---- (3) Off-principal-axis Eulerian nutation -------- *
     * Symmetric top: Ixx = Iyy = 1, Izz = 2.
     * ω₀ = (0.1, 0, 1.0). Polhode period
     *   T = 2π · I_t / (ωz · |I_s − I_t|) = 2π · 1 / (1 · 1) = 2π s.
     * After one full period the body-frame ω should return to its
     * starting value within the first-order integrator's drift. */
    K26AstroAttitudeStateExt s;
    k26astro_attitude_init_ext(&s, diag_(1.0, 1.0, 2.0));
    s.omega_body.x = 0.1;
    s.omega_body.y = 0.0;
    s.omega_body.z = 1.0;
    double T = 2.0 * 3.14159265358979323846;
    int    N = 100000;
    double dt = T / (double)N;
    for (int n = 0; n < N; n++) {
        k26astro_attitude_step_free_ext(&s, dt);
    }
    /* ωz is exactly conserved analytically (decoupled axis). */
    assert(fabs(s.omega_body.z - 1.0) < 1e-9);
    /* ωx returns to 0.1 within first-order accumulated error
     * (∝ dt · N · ω_perp² ≈ 6.3e-6 with these settings). */
    assert(fabs(s.omega_body.x - 0.1) < 5e-4);
    assert(fabs(s.omega_body.y - 0.0) < 5e-4);
    /* Quaternion stays on the unit sphere. */
    double qn = sqrt(s.q.x * s.q.x + s.q.y * s.q.y + s.q.z * s.q.z + s.q.w * s.q.w);
    assert(near_(qn, 1.0, 1e-10));
    k26astro_attitude_destroy_ext(&s);

    /* ---- (4) Constant torque, principal axis ------------- */
    K26AstroAttitudeStateExt t1;
    k26astro_attitude_init_ext(&t1, diag_(2.0, 3.0, 5.0));
    K26V3 tau = { 4.0, 0.0, 0.0 };
    /* ω₀ = 0; ω̇x = τx / Ixx = 2. After dt = 0.01: Δωx = 0.02. */
    k26astro_attitude_step_torque_ext(&t1, tau, 0.01);
    assert(near_(t1.omega_body.x, 0.02, 1e-12));
    assert(fabs(t1.omega_body.y) < 1e-15);
    assert(fabs(t1.omega_body.z) < 1e-15);
    k26astro_attitude_destroy_ext(&t1);

    /* ---- (5) Update inertia recomputes inverse ----------- */
    K26AstroAttitudeStateExt u;
    k26astro_attitude_init_ext(&u, diag_(1.0, 1.0, 1.0));
    assert(near_(u.inertia_inverse.m[0][0], 1.0, 1e-12));
    k26astro_attitude_update_inertia(&u, diag_(10.0, 20.0, 40.0));
    assert(near_(u.inertia.m[0][0],         10.0,    1e-12));
    assert(near_(u.inertia_inverse.m[0][0],  0.1,    1e-12));
    assert(near_(u.inertia_inverse.m[1][1],  0.05,   1e-12));
    assert(near_(u.inertia_inverse.m[2][2],  0.025,  1e-12));
    k26astro_attitude_destroy_ext(&u);

    /* ---- (6) Ext registry dispatch ----------------------- */
    K26AstroAttitudeStateExt r;
    k26astro_attitude_init_ext(&r, diag_(1.0, 1.0, 1.0));
    assert(r.torques == NULL);
    assert(k26astro_attitude_register_torque_ext(&r, torque_const_x_, NULL) == 0);
    /* Step from ω₀ = 0; expect Δωx = 1.5 · 0.01 / 1.0 = 0.015. */
    k26astro_attitude_step_torque_registry_ext(&r, 0.0, 0.01);
    assert(near_(r.omega_body.x, 0.015, 1e-12));
    k26astro_attitude_clear_torques_ext(&r);
    assert(r.torques == NULL);
    k26astro_attitude_destroy_ext(&r);

    /* ---- Null guards ------------------------------------- */
    k26astro_attitude_init_ext(NULL, diag_(1, 1, 1));
    k26astro_attitude_destroy_ext(NULL);
    k26astro_attitude_step_free_ext(NULL, 0.0);
    k26astro_attitude_step_torque_ext(NULL, (K26V3){0, 0, 0}, 0.0);
    k26astro_attitude_update_inertia(NULL, diag_(1, 1, 1));
    assert(k26astro_attitude_register_torque_ext(NULL, torque_const_x_, NULL) == 1);
    k26astro_attitude_clear_torques_ext(NULL);
    k26astro_attitude_step_torque_registry_ext(NULL, 0.0, 0.0);

    printf("test_attitude_step_ext: OK\n");
    return 0;
}
