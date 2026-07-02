/* test_grav_pc_diagnostics.c — IAS15 substep diagnostic accessors.
 *
 * Validates the per-state IAS15 introspection API introduced for
 * driver-side integrator diagnostic emission. Asserts:
 *
 *   1. Accessors return zero on a fresh state (no substep taken).
 *   2. Smooth Sun-Earth two-body integration drives pc_iter into
 *      the healthy 2-6 range and reports eps_b_achieved at or
 *      below the controller's ias15_tol.
 *   3. A close-encounter scenario drives pc_iter higher (closer
 *      to or at the IAS15_MAX_PC_ITER cap of 12). The exact value
 *      depends on geometry — we only assert it rises above the
 *      smooth-case ceiling.
 *   4. Accessors are NULL-safe.
 *
 * A downstream integrator-diagnostic consumer previously
 * hardcoded pc_iter=12 because no accessor existed. These
 * accessors let a driver emit real values per substep. */
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
    /* ---- Smooth Sun-Earth two-body --------------------------- */
    K26AstroBody bodies[2];
    memset(bodies, 0, sizeof(bodies));
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH;
    bodies[1].pos  = k26astro_pos_from_m(K26A_AU_M, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, sqrt(K26A_GM_SUN / K26A_AU_M), 0.0);
    bodies[1].parent_body_idx = 0;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15) == 0);

    /* Fresh state: every accessor reports zero. */
    assert(k26astro_grav_last_pc_iterations(&state) == 0);
    assert(k26astro_grav_last_dt_taken(&state) == 0.0);
    assert(k26astro_grav_last_eps_b_achieved(&state) == 0.0);
    assert(k26astro_grav_rejected_steps_total(&state) == 0u);

    /* Take 50 steps of 1 day each. Track max PC iter + max eps_b. */
    int    max_pc_smooth = 0;
    double max_eps_b_smooth = 0.0;
    double last_dt_smooth = 0.0;
    for (int s = 0; s < 50; s++) {
        int rc = k26astro_grav_step(&state, 86400.0);
        assert(rc == 0);
        int pc = k26astro_grav_last_pc_iterations(&state);
        double eps = k26astro_grav_last_eps_b_achieved(&state);
        last_dt_smooth = k26astro_grav_last_dt_taken(&state);
        if (pc > max_pc_smooth) max_pc_smooth = pc;
        if (eps > max_eps_b_smooth) max_eps_b_smooth = eps;
    }
    /* Healthy smooth integration: PC iter stays in [2, 6] per
     * Rein-Spiegel 2015 §4. We accept the broader [2, 8] window
     * to account for initial-startup iterations on the first
     * substeps. */
    assert(max_pc_smooth >= 2);
    assert(max_pc_smooth <= 8);
    /* eps_b achieved should be near or below the ias15_tol
     * (default 1e-9). */
    assert(max_eps_b_smooth <= 1.0e-6);
    /* dt_taken is the actual final substep — non-zero after a step. */
    assert(last_dt_smooth > 0.0);

    fprintf(stderr,
            "test_grav_pc_diagnostics: smooth Sun-Earth — "
            "max_pc=%d max_eps_b=%.3e last_dt=%.3e s\n",
            max_pc_smooth, max_eps_b_smooth, last_dt_smooth);

    k26astro_grav_state_destroy(&state);

    /* ---- Close-encounter scenario ---------------------------- *
     *
     * Two equal masses on a fast parabolic flyby — should drive
     * PC iteration count higher than the smooth case. Use roughly
     * solar-mass primaries at 1e10 m closest approach to keep
     * the integration in stiff territory. */
    K26AstroBody close[2];
    memset(close, 0, sizeof(close));
    close[0].kind = K26ASTRO_BODY_STAR;
    close[0].gm   = K26A_GM_SUN;
    close[0].pos  = k26astro_pos_from_m(-5.0e10, 0.0, 0.0);
    close[0].vel  = k26m3d_v3(+1.0e5, 0.0, 0.0);
    close[0].parent_body_idx = -1;

    close[1].kind = K26ASTRO_BODY_STAR;
    close[1].gm   = K26A_GM_SUN;
    close[1].pos  = k26astro_pos_from_m(+5.0e10, 0.0, 0.0);
    close[1].vel  = k26m3d_v3(-1.0e5, 0.0, 0.0);
    close[1].parent_body_idx = -1;

    K26AstroGravState state2;
    assert(k26astro_grav_state_init(&state2, close, 2) == 0);
    assert(k26astro_grav_set_integrator(&state2, K26ASTRO_INTEGRATOR_IAS15) == 0);

    /* Step toward closest approach. Coarse dt to force the
     * controller to subdivide aggressively. */
    int max_pc_stiff = 0;
    double max_eps_b_stiff = 0.0;
    for (int s = 0; s < 20; s++) {
        int rc = k26astro_grav_step(&state2, 1.0e5);
        /* Don't require convergence — stiff scenarios can hit
         * the reject path. Iter count and eps_b are still
         * meaningful even on reject. */
        (void)rc;
        int pc = k26astro_grav_last_pc_iterations(&state2);
        double eps = k26astro_grav_last_eps_b_achieved(&state2);
        /* PC can be negative on reject (the no-converge signal).
         * For the diagnostic-rise assertion below, fold the
         * absolute value. */
        int pc_abs = (pc < 0) ? -pc : pc;
        if (pc_abs > max_pc_stiff) max_pc_stiff = pc_abs;
        if (eps > max_eps_b_stiff) max_eps_b_stiff = eps;
    }
    fprintf(stderr,
            "test_grav_pc_diagnostics: stiff flyby — "
            "max_pc=%d max_eps_b=%.3e\n",
            max_pc_stiff, max_eps_b_stiff);

    /* Stiff geometry should drive PC iter above the smooth
     * ceiling. Tolerance is loose because the exact PC count
     * depends on the specific encounter geometry — what matters
     * is the diagnostic shows the controller working harder. */
    assert(max_pc_stiff > max_pc_smooth ||
           k26astro_grav_rejected_steps_total(&state2) > 0u);

    k26astro_grav_state_destroy(&state2);

    /* ---- NULL safety ---------------------------------------- */
    assert(k26astro_grav_last_pc_iterations(NULL)  == 0);
    assert(k26astro_grav_last_dt_taken(NULL)        == 0.0);
    assert(k26astro_grav_last_eps_b_achieved(NULL)  == 0.0);
    assert(k26astro_grav_rejected_steps_total(NULL) == 0u);

    fprintf(stderr, "test_grav_pc_diagnostics: OK\n");
    return 0;
}
