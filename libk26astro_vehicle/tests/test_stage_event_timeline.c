/* test_stage_event_timeline.c — cross-lib integration gate for the
 * stage-event timeline + libk26astro_grav event-time root-finder.
 *
 * Builds a minimal 2-body gravity state (massive central + small
 * spacecraft body bound to the vehicle), schedules three stage
 * events out of epoch order, flushes them to the gravity event
 * registry, integrates across the firing window, and asserts:
 *
 *   - events sorted ascending by epoch on insertion;
 *   - mass_at(t) walks the timeline correctly;
 *   - all three events fired during integration;
 *   - vehicle basic mass equals the cumulative drops;
 *   - bound body gm and mass updated via k26astro_body_set_mass;
 *   - latest inertia tensor diagonal matches the latest event's
 *     post-event values, off-diagonals preserved from the initial
 *     set. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/stage_event.h"
#include "k26astro_grav/grav.h"
#include "k26astro_body/body.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static K26AstroEpoch ep_plus_(K26AstroEpoch base, double delta_s)
{
    k26astro_epoch_add_seconds(&base, delta_s);
    return base;
}

int main(void)
{
    /* ---- Vehicle setup ---------------------------------------- */
    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);
    k26astro_vehicle_set_dry_mass(v, 1000.0);
    /* Full inertia with off-diagonals so we can verify they're
     * preserved across the stage events' diagonal-only updates. */
    K26M3 J;
    J.m[0][0] = 100.0; J.m[0][1] =  10.0; J.m[0][2] =   5.0;
    J.m[1][0] =  10.0; J.m[1][1] = 200.0; J.m[1][2] =  -3.0;
    J.m[2][0] =   5.0; J.m[2][1] =  -3.0; J.m[2][2] = 300.0;
    k26astro_vehicle_set_inertia_full(v, J);

    /* ---- Two bodies ------------------------------------------- *
     * bodies[0] : central — small mass, fixed at origin.
     * bodies[1] : spacecraft — bound to the vehicle. Tiny but
     *             non-zero orbital velocity so the integrator has
     *             something to do; the exact orbit shape doesn't
     *             matter for the gate. */
    K26AstroBody bodies[2];
    k26astro_body_init(&bodies[0]);
    k26astro_body_init(&bodies[1]);

    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].mass = 1.0;
    bodies[0].gm   = K26A_G;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_SPACECRAFT;
    bodies[1].mass = 1000.0;
    bodies[1].gm   = K26A_G * 1000.0;
    bodies[1].pos  = k26astro_pos_from_m(1.0e6, 0.0, 0.0);
    bodies[1].vel.x = 0.0;
    bodies[1].vel.y = 1.0e-3;
    bodies[1].vel.z = 0.0;
    bodies[1].parent_body_idx = 0;

    k26astro_vehicle_bind_body(v, &bodies[1]);

    K26AstroGravState st;
    int rc = k26astro_grav_state_init(&st, bodies, 2);
    assert(rc == K26ASTRO_E_OK);
    /* RK4 has no special integrator carry so it's a clean baseline
     * for the event-bisection wrapper. */
    rc = k26astro_grav_set_integrator(&st, K26ASTRO_INTEGRATOR_RK4);
    assert(rc == K26ASTRO_E_OK);

    /* ---- Schedule three events out of epoch order ------------- */
    K26AstroEpoch base = k26astro_epoch_j2000_tt();
    K26AstroEpoch t_a = ep_plus_(base, 1.0);
    K26AstroEpoch t_b = ep_plus_(base, 2.0);
    K26AstroEpoch t_c = ep_plus_(base, 3.0);

    int ib = k26astro_vehicle_schedule_stage_event(v, t_b, 100.0,
                                                    80.0, 180.0, 280.0);
    int ia = k26astro_vehicle_schedule_stage_event(v, t_a,  50.0,
                                                    90.0, 190.0, 290.0);
    int ic = k26astro_vehicle_schedule_stage_event(v, t_c, 200.0,
                                                    50.0, 150.0, 250.0);

    /* The schedule call returns the post-insertion sorted index. */
    assert(ib == 0);   /* first scheduled — goes to index 0 */
    assert(ia == 0);   /* earlier epoch — inserted at front */
    assert(ic == 2);   /* latest epoch — appended */

    assert(k26astro_vehicle_event_count(v) == 3);
    K26AstroEpoch e0 = k26astro_vehicle_event_epoch_at(v, 0);
    K26AstroEpoch e1 = k26astro_vehicle_event_epoch_at(v, 1);
    K26AstroEpoch e2 = k26astro_vehicle_event_epoch_at(v, 2);
    assert(k26astro_epoch_diff_seconds(&e0, &t_a) == 0.0);
    assert(k26astro_epoch_diff_seconds(&e1, &t_b) == 0.0);
    assert(k26astro_epoch_diff_seconds(&e2, &t_c) == 0.0);

    /* ---- mass_at query walks the sorted timeline ------------- */
    K26AstroEpoch t_after_a = ep_plus_(base, 1.5);
    assert(k26astro_vehicle_mass_at(v, t_after_a) == 950.0);
    K26AstroEpoch t_after_b = ep_plus_(base, 2.5);
    assert(k26astro_vehicle_mass_at(v, t_after_b) == 850.0);
    K26AstroEpoch t_after_c = ep_plus_(base, 3.5);
    assert(k26astro_vehicle_mass_at(v, t_after_c) == 650.0);
    K26AstroEpoch t_before = ep_plus_(base, 0.5);
    assert(k26astro_vehicle_mass_at(v, t_before) == 1000.0);

    /* ---- Flush to the gravity event registry ----------------- */
    rc = k26astro_vehicle_register_stage_events_with(v, &st);
    assert(rc == K26ASTRO_E_OK);

    /* ---- Integrate across [base, base + 4 s] in 0.5 s steps --- */
    for (int n = 0; n < 8; n++) {
        int sc = k26astro_grav_step(&st, 0.5);
        assert(sc == K26ASTRO_E_OK);
    }

    /* ---- Assert all three events fired ----------------------- *
     *
     * Final basic mass = 1000 - 50 - 100 - 200 = 650. */
    assert(k26astro_vehicle_mass_now(v) == 650.0);

    /* Body GM tracks via k26astro_body_set_mass. */
    assert(bodies[1].mass == 650.0);
    assert(fabs(bodies[1].gm - K26A_G * 650.0) < 1e-25);

    /* Latest inertia tensor: diagonal from event C; off-diagonals
     * preserved from the initial set_inertia_full. */
    K26AstroEpoch t_end = ep_plus_(base, 4.0);
    K26M3 I_final = k26astro_vehicle_inertia_at(v, t_end);
    assert(I_final.m[0][0] == 50.0);
    assert(I_final.m[1][1] == 150.0);
    assert(I_final.m[2][2] == 250.0);
    assert(I_final.m[0][1] == 10.0);
    assert(I_final.m[1][0] == 10.0);
    assert(I_final.m[1][2] == -3.0);

    /* Live Ext attitude state's inertia matches the post-event state. */
    K26AstroAttitudeStateExt *a = k26astro_vehicle_attitude_ext(v);
    assert(a->inertia.m[0][0] == 50.0);
    assert(a->inertia.m[1][1] == 150.0);
    assert(a->inertia.m[2][2] == 250.0);
    /* Off-diagonals preserved through every stage event. */
    assert(a->inertia.m[0][1] == 10.0);
    assert(a->inertia.m[1][2] == -3.0);
    /* Inverse recomputed via k26astro_attitude_update_inertia —
     * verify I . I_inv = identity. */
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) {
                s += a->inertia.m[r][k] * a->inertia_inverse.m[k][c];
            }
            double expected = (r == c) ? 1.0 : 0.0;
            assert(fabs(s - expected) < 1e-12);
        }
    }

    k26astro_grav_state_destroy(&st);
    k26astro_vehicle_destroy(v);

    printf("test_stage_event_timeline: OK\n");
    return 0;
}
