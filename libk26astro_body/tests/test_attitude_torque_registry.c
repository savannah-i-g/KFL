/* test_attitude_torque_registry.c — multi-source torque-registry
 * gate.
 *
 * Registers three torque sources whose closed-form sum is known by
 * construction, verifies:
 *
 *   (1) k26astro_attitude_step_torque_registry sums the contributions
 *       and advances ω by the analytic Δω = (Σ τ / I) · dt;
 *   (2) the registry survives clear + re-register cycles;
 *   (3) k26astro_attitude_destroy is a no-op on a state that never
 *       registered, and frees cleanly on a populated one.
 *
 * The three synthetic sources stand in for gravity-gradient + SRP +
 * thrust-vector-control contributions. The numerical content is the
 * registry, not the physics — the per-source torques are constants
 * chosen for an unambiguous sum. */
#include "k26astro_body/attitude.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

/* Source A — fixed body-frame torque along +x. */
static void torque_a(const K26AstroAttitudeState *state, double t,
                     K26V3 *out, void *ctx)
{
    (void)state; (void)t; (void)ctx;
    out->x += 1.0;
}

/* Source B — fixed body-frame torque along +y. */
static void torque_b(const K26AstroAttitudeState *state, double t,
                     K26V3 *out, void *ctx)
{
    (void)state; (void)t; (void)ctx;
    out->y += 2.0;
}

/* Source C — fixed body-frame torque along -z, ctx is the magnitude.
 * Exercises the ctx pointer end-to-end. */
static void torque_c(const K26AstroAttitudeState *state, double t,
                     K26V3 *out, void *ctx)
{
    (void)state; (void)t;
    double *mag = (double *)ctx;
    out->z -= *mag;
}

int main(void)
{
    /* ---- Lifecycle: destroy on never-registered is a no-op ---- */
    K26AstroAttitudeState empty;
    k26astro_attitude_init(&empty, 1.0, 1.0);
    assert(empty.torques == NULL);
    k26astro_attitude_destroy(&empty);
    assert(empty.torques == NULL);

    /* ---- Three-source registration + dispatch ---------------- */
    K26AstroAttitudeState a;
    k26astro_attitude_init(&a, 1.0, 1.0);
    /* Override solid-sphere inertia for a clean Δω calc:
     *   I = diag(1, 2, 4) kg·m². */
    a.inertia_diag.x = 1.0;
    a.inertia_diag.y = 2.0;
    a.inertia_diag.z = 4.0;
    a.omega_body.x = 0.0;
    a.omega_body.y = 0.0;
    a.omega_body.z = 0.0;

    double c_mag = 3.0;
    assert(k26astro_attitude_register_torque(&a, torque_a, NULL) == 0);
    assert(k26astro_attitude_register_torque(&a, torque_b, NULL) == 0);
    assert(k26astro_attitude_register_torque(&a, torque_c, &c_mag) == 0);
    assert(a.torques != NULL);

    /* Expected summed torque: (1, 2, -3) N·m.
     * Closed-form Δω over dt with ω₀ = 0 and diagonal I:
     *   Δωx = τx / Ixx · dt = 1 / 1 · dt
     *   Δωy = τy / Iyy · dt = 2 / 2 · dt
     *   Δωz = τz / Izz · dt = -3 / 4 · dt
     * The Euler cross-coupling term (Iyy - Izz) ωy ωz vanishes at
     * ω₀ = 0 so this matches the registry sum exactly for one step. */
    const double dt = 0.01;
    k26astro_attitude_step_torque_registry(&a, 0.0, dt);

    assert(near_(a.omega_body.x,  1.0  / 1.0 * dt, 1e-12));
    assert(near_(a.omega_body.y,  2.0  / 2.0 * dt, 1e-12));
    assert(near_(a.omega_body.z, -3.0  / 4.0 * dt, 1e-12));

    /* ---- Clear + re-register: registry is reusable ----------- */
    k26astro_attitude_clear_torques(&a);
    assert(a.torques == NULL);

    /* Re-register only sources A and B; advance from the post-first-
     * step ω. The third axis must NOT change because torque_c is no
     * longer in the registry. */
    K26V3 omega_before_second = a.omega_body;
    assert(k26astro_attitude_register_torque(&a, torque_a, NULL) == 0);
    assert(k26astro_attitude_register_torque(&a, torque_b, NULL) == 0);

    k26astro_attitude_step_torque_registry(&a, dt, dt);

    /* x and y advance by their per-source torques plus the Euler
     * cross-coupling (Iyy - Izz)·ωy·ωz / Ixx and (Izz - Ixx)·ωz·ωx /
     * Iyy at the post-first-step ω ≈ (1e-2, 1e-2, -7.5e-3). The
     * cross-coupling bound is |ΔI| · |ω · ω| · dt / I ≈ 1.5e-6.
     * z is unchanged by the registry but receives the
     * (Ixx - Iyy)·ωx·ωy / Izz cross-coupling. */
    double dx = a.omega_body.x - omega_before_second.x;
    double dy = a.omega_body.y - omega_before_second.y;
    double dz = a.omega_body.z - omega_before_second.z;
    assert(fabs(dx - 1.0 / 1.0 * dt) < 1.0e-5);
    assert(fabs(dy - 2.0 / 2.0 * dt) < 1.0e-5);
    /* |Δωz| ≤ |Ixx - Iyy| · |ωx · ωy| · dt / Izz
     *      ≈ |1 - 2| · 1e-2 · 1e-2 · 1e-2 / 4 = 2.5e-7. */
    assert(fabs(dz) < 1.0e-6);

    /* ---- Capacity growth: register beyond the initial cap (4) - */
    K26AstroAttitudeState g;
    k26astro_attitude_init(&g, 1.0, 1.0);
    for (int i = 0; i < 10; i++) {
        assert(k26astro_attitude_register_torque(&g, torque_a, NULL) == 0);
    }
    /* Step with 10 copies of source A summing to (10, 0, 0). */
    g.inertia_diag.x = g.inertia_diag.y = g.inertia_diag.z = 1.0;
    g.omega_body.x = g.omega_body.y = g.omega_body.z = 0.0;
    k26astro_attitude_step_torque_registry(&g, 0.0, 1.0);
    assert(near_(g.omega_body.x, 10.0, 1e-12));
    k26astro_attitude_destroy(&g);

    /* ---- Null inputs -------------------------------------- */
    assert(k26astro_attitude_register_torque(NULL, torque_a, NULL) == 1);
    assert(k26astro_attitude_register_torque(&a, NULL, NULL) == 1);
    k26astro_attitude_clear_torques(NULL); /* no crash */
    k26astro_attitude_destroy(NULL);       /* no crash */
    k26astro_attitude_step_torque_registry(NULL, 0.0, 0.0); /* no crash */

    /* ---- Final teardown ----------------------------------- */
    k26astro_attitude_destroy(&a);
    assert(a.torques == NULL);

    printf("test_attitude_torque_registry: OK\n");
    return 0;
}
