/* test_ias15_wall_budget.c — IAS15 wall-time budget safety net.
 *
 * Verifies that k26astro_grav_ias15_set_wall_budget bounds an
 * otherwise-stalling integration call. Without the budget,
 * pathological chaotic regions (e.g. Burrau-style triple close
 * encounters at impossibly tight tol) can spend minutes accepting
 * picosecond substeps in the IEEE-754 truncation-noise regime;
 * the budget exists as a developer-facing safety net.
 *
 * Strategy: construct a near-degenerate 3-body close encounter,
 * set ias15_tol absurdly tight (1e-16) so the controller is forced
 * to shrink dt aggressively, set a 100 ms wall budget, and assert
 * the step returns K26ASTRO_E_TIME_BUDGET (not OK, not
 * NO_CONVERGE).
 *
 * The budget is NOT a determinism gate; wall-clock measurement is
 * platform-dependent. Deterministic runs must keep budget = 0.0. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroBody bodies[3];
    memset(bodies, 0, sizeof(bodies));

    /* A central mass with two test particles in close encounter. The
     * specific configuration mimics Burrau-style chaos: nearly
     * radial trajectories with a near-tangent intercept. */
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = (K26V3){ 0.0, 0.0, 0.0 };
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH * 100.0;   /* Jupiter-ish mass */
    bodies[1].pos  = k26astro_pos_from_m(1.0 * K26A_AU_M, 0.0, 0.0);
    bodies[1].vel  = (K26V3){ 0.0, 4.0e4, 0.0 };
    bodies[1].parent_body_idx = 0;

    bodies[2].kind = K26ASTRO_BODY_PLANET;
    bodies[2].gm   = K26A_GM_EARTH * 100.0;
    /* Co-located within 1000 km — forces close-encounter dynamics. */
    bodies[2].pos  = k26astro_pos_from_m(1.0 * K26A_AU_M + 1.0e6, 0.0, 0.0);
    bodies[2].vel  = (K26V3){ 0.0, 4.0e4, 0.0 };
    bodies[2].parent_body_idx = 0;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 3) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15) == 0);
    /* Absurdly tight tolerance forces dt to collapse. */
    k26astro_grav_ias15_set_tol(&state, 1.0e-16);
    /* 100 ms wall budget — well below the chaos-stall timescale. */
    k26astro_grav_ias15_set_wall_budget(&state, 0.1);

    /* Request 1 year. Without the budget, this would stall for many
     * seconds in chaotic-region picosecond stepping. With the
     * budget, returns TIME_BUDGET cleanly within ~100 ms. */
    int rc = k26astro_grav_step(&state, 365.25 * 86400.0);
    fprintf(stderr,
            "test_ias15_wall_budget: rc=%d (TIME_BUDGET=%d) rejects=%u\n",
            rc, K26ASTRO_E_TIME_BUDGET,
            k26astro_grav_ias15_rejected_steps(&state));
    assert(rc == K26ASTRO_E_TIME_BUDGET);

    /* Sanity: with budget = 0 the call falls back to normal
     * behaviour (NO_CONVERGE or completion). Re-init state so the
     * carry-over from the budget-exit doesn't bias this leg. */
    K26AstroBody bodies2[3];
    memcpy(bodies2, bodies, sizeof(bodies));
    K26AstroGravState state2;
    assert(k26astro_grav_state_init(&state2, bodies2, 3) == 0);
    assert(k26astro_grav_set_integrator(&state2, K26ASTRO_INTEGRATOR_IAS15) == 0);
    k26astro_grav_ias15_set_tol(&state2, 1.0e-9);
    k26astro_grav_ias15_set_wall_budget(&state2, 0.0);
    /* Short call so we don't exercise the chaos stall — just verify
     * budget = 0 doesn't short-circuit normal integration. */
    rc = k26astro_grav_step(&state2, 86400.0);
    fprintf(stderr,
            "test_ias15_wall_budget: budget=0 1-day step rc=%d\n", rc);
    assert(rc == K26ASTRO_E_OK || rc == K26ASTRO_E_NO_CONVERGE);
    assert(rc != K26ASTRO_E_TIME_BUDGET);

    k26astro_grav_state_destroy(&state);
    k26astro_grav_state_destroy(&state2);

    printf("test_ias15_wall_budget: OK\n");
    return 0;
}
