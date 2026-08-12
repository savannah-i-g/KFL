/* test_infostate_light_time.c — retarded-time fixed-point gate.
 *
 * Covers the canonical light-time scenarios:
 *   (1) Stationary target at 1 ls (~3e8 m): τ = 1 s exactly.
 *   (2) Stationary target at 1000 km: τ ≈ 3.34e-3 s exactly.
 *   (3) Inbound target moving at 0.01 c: τ_self-consistent satisfies
 *       |x_t(t - τ) - x_o(t)| = c · τ, within fixed-point tolerance.
 *   (4) History-too-short returns valid=0.
 *
 * The convergence-iter count is recorded; closed-form geometry
 * should converge in 2-3 iterations regardless of regime. */

#include "k26astro_infostate/infostate.h"
#include "k26astro_infostate/infostate_consts.h"

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static K26AstroEpoch epoch_j2000_plus_s_(double dt_s)
{
    K26AstroEpoch e = k26astro_epoch_j2000_tt();
    k26astro_epoch_add_seconds(&e, dt_s);
    return e;
}

static void set_body_xyz_(K26AstroBody *b, double x, double y, double z)
{
    b->pos = k26astro_pos_from_m(x, y, z);
}

/* Push N+1 samples evenly spanning [t_min, t_max] with the final
 * sample forced to exactly t_max — avoids FP drift in dt += step
 * loops, which would otherwise skip the t_max boundary and leave
 * the retarded-time lookup outside the history window. */
static void push_uniform_(K26AstroInfostate *s,
                          const K26AstroVehicle *target,
                          int n_steps,
                          double t_min_s,
                          double t_max_s,
                          K26V3 r0,
                          K26V3 v_inbound,
                          int moving)
{
    for (int i = 0; i <= n_steps; i++) {
        double dt = (i == n_steps)
            ? t_max_s
            : (t_min_s + (t_max_s - t_min_s) *
                         (double)i / (double)n_steps);
        K26V3 r = r0;
        K26V3 v = { 0.0, 0.0, 0.0 };
        if (moving) {
            r.x = r0.x + v_inbound.x * dt;
            r.y = r0.y + v_inbound.y * dt;
            r.z = r0.z + v_inbound.z * dt;
            v = v_inbound;
        }
        k26astro_infostate_target_push(s, target,
            epoch_j2000_plus_s_(dt), r, v);
    }
}

int main(void)
{
    K26AstroBody obs_body, tgt_body;
    memset(&obs_body, 0, sizeof(obs_body));
    memset(&tgt_body, 0, sizeof(tgt_body));

    K26AstroVehicle *observer = k26astro_vehicle_new();
    K26AstroVehicle *target   = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(observer, &obs_body);
    k26astro_vehicle_bind_body(target,   &tgt_body);

    K26AstroInfostate *s = k26astro_infostate_new(observer, 256);
    assert(s != NULL);

    const double c = K26A_C;
    const double tol_range = 1.0;  /* metres — well below ulp-floor at AU scale */

    /* ---- Scenario 1: stationary target at 1 ls. ----------------- */
    /* Light-second distance: τ should be exactly 1 s. */
    set_body_xyz_(&obs_body, 0.0, 0.0, 0.0);
    K26V3 r_target_1ls = { c, 0.0, 0.0 };
    K26V3 v_zero       = { 0.0, 0.0, 0.0 };

    push_uniform_(s, target, 20, -2.0, 0.0, r_target_1ls, v_zero, 0);

    K26AstroInfostateObservation o1 =
        k26astro_infostate_observe(s, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_IR);
    assert(o1.valid == 1);
    assert(fabs(o1.age_s - 1.0) < 1.0e-9);
    assert(fabs(o1.range_m - c) < tol_range);

    /* ---- Scenario 2: stationary target at 1000 km. ------------- */
    K26V3 r_target_1000km = { 1.0e6, 0.0, 0.0 };
    /* Fresh ring per scenario; tracks live inside the infostate. */
    K26AstroInfostate *s2 = k26astro_infostate_new(observer, 256);
    push_uniform_(s2, target, 100, -1.0, 0.0, r_target_1000km, v_zero, 0);
    K26AstroInfostateObservation o2 =
        k26astro_infostate_observe(s2, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_RADAR);
    assert(o2.valid == 1);
    const double tau_expected = 1.0e6 / c;
    assert(fabs(o2.age_s - tau_expected) < 1.0e-12);
    assert(fabs(o2.range_m - 1.0e6) < tol_range);

    /* ---- Scenario 3: inbound target at 0.01 c. ----------------- */
    /* Target at r0 = 1e8 m at t=0, moving inbound at 0.01c along x:
     *   r(t) = (1e8 + v_x * t, 0, 0), with v_x = -0.01 c.
     * Light-time τ satisfies (1e8 + v_x * (-τ)) = c τ at t=0:
     *   1e8 + 0.01c · τ = c τ
     *   τ = 1e8 / (0.99 c). */
    const double v_inbound_scalar = -0.01 * c;
    K26V3 v_inbound = { v_inbound_scalar, 0.0, 0.0 };
    K26V3 r0_inbound = { 1.0e8, 0.0, 0.0 };
    K26AstroInfostate *s3 = k26astro_infostate_new(observer, 5001);
    push_uniform_(s3, target, 5000, -5.0, 0.0, r0_inbound, v_inbound, 1);
    K26AstroInfostateObservation o3 =
        k26astro_infostate_observe(s3, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_IR);
    assert(o3.valid == 1);
    const double tau_inbound = 1.0e8 / (0.99 * c);
    /* Push samples are at 1 ms resolution; linear interp introduces
     * O(v_inbound * 0.5 ms) ≈ 0.5 * 0.01c * 1e-3 ≈ 1.5 km positional
     * uncertainty per sample interval. Tolerance reflects that
     * sampling floor, not the solver's own convergence. */
    assert(fabs(o3.age_s - tau_inbound) < 1.0e-5);
    /* Range at retarded time matches c · τ to ulp. */
    assert(fabs(o3.range_m - c * o3.age_s) <
           K26ASTRO_INFOSTATE_FIXED_POINT_TOL_S * c + 1.0);

    /* ---- Scenario 4: history too short, returns invalid. ------ */
    K26AstroInfostate *s4 = k26astro_infostate_new(observer, 16);
    /* Push only samples from t = -0.01 .. 0; target is 1 ls away so
     * τ = 1 s but oldest sample is -0.01 s. Should return invalid. */
    push_uniform_(s4, target, 10, -0.01, 0.0, r_target_1ls, v_zero, 0);
    K26AstroInfostateObservation o4 =
        k26astro_infostate_observe(s4, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_IR);
    assert(o4.valid == 0);

    /* ---- Scenario 5: bit-identical replay across runs. -------- *
     * Two independent infostates fed the same push sequence must
     * return bit-identical observations. Tests the determinism
     * gate at the unit level. */
    K26AstroInfostate *sA = k26astro_infostate_new(observer, 256);
    K26AstroInfostate *sB = k26astro_infostate_new(observer, 256);
    /* Lateral motion via a non-axis-aligned vector. */
    K26V3 v_inbound_lat = { v_inbound_scalar,
                            v_inbound_scalar * 0.1, 0.0 };
    K26V3 r0_lat        = { 1.0e7, 0.0, 0.0 };
    push_uniform_(sA, target, 40, -2.0, 0.0, r0_lat, v_inbound_lat, 1);
    push_uniform_(sB, target, 40, -2.0, 0.0, r0_lat, v_inbound_lat, 1);
    K26AstroInfostateObservation oA =
        k26astro_infostate_observe(sA, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_IR);
    K26AstroInfostateObservation oB =
        k26astro_infostate_observe(sB, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_IR);
    assert(oA.valid == 1 && oB.valid == 1);
    /* Bit-identical: both observations must memcmp-equal in their
     * numerical fields. */
    assert(oA.age_s    == oB.age_s);
    assert(oA.range_m  == oB.range_m);
    assert(oA.position.x == oB.position.x);
    assert(oA.position.y == oB.position.y);
    assert(oA.position.z == oB.position.z);
    assert(oA.velocity.x == oB.velocity.x);
    assert(oA.velocity.y == oB.velocity.y);
    assert(oA.velocity.z == oB.velocity.z);

    k26astro_infostate_destroy(s);
    k26astro_infostate_destroy(s2);
    k26astro_infostate_destroy(s3);
    k26astro_infostate_destroy(s4);
    k26astro_infostate_destroy(sA);
    k26astro_infostate_destroy(sB);
    k26astro_vehicle_destroy(observer);
    k26astro_vehicle_destroy(target);

    printf("test_infostate_light_time: OK\n");
    return 0;
}
