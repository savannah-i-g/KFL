/* burrau_gen.c — generate the Aarseth 1973 Burrau 3-4-5 Pythagorean
 * 3-body reference trajectory via REBOUND-IAS15.
 *
 * Output: CSV to stdout, in the exact format `tests/astro/burrau_ias15`
 * expects:
 *
 *     # t, x0, y0, x1, y1, x2, y2
 *     0.000000, 1.0, 3.0, -2.0, -1.0, 1.0, -1.0
 *     0.100000, 0.9999..., 2.9999..., ...
 *     ...
 *     70.000000, ...
 *
 * Sampling: every 0.1 normalised time units, T=0..70 inclusive
 * (701 rows + 1 header). G=1, all masses dimensionless, all
 * positions dimensionless.
 *
 * Build (REBOUND >= 4.0):
 *     gcc -O2 -o burrau_gen burrau_gen.c -lrebound -lm
 *
 * For pre-v4 REBOUND, the API names drop the "_simulation" infix:
 * change `reb_simulation_create` → `reb_create_simulation`,
 *        `reb_simulation_add`    → `reb_add`,
 *        `reb_simulation_integrate` → `reb_integrate`,
 *        `reb_simulation_free`   → `reb_free_simulation`.
 *
 * Run:
 *     ./burrau_gen > burrau_rebound_T70.csv
 *
 * Then the astro test suite flips burrau_ias15 from
 * SKIP to PASS (or FAIL if the K26 + REBOUND IAS15 trajectories
 * diverge beyond the gate's tolerance — itself a useful signal). */
#include "rebound.h"
#include <stdio.h>

int main(void)
{
    struct reb_simulation *sim = reb_simulation_create();
    if (!sim) {
        fprintf(stderr, "burrau_gen: reb_simulation_create failed\n");
        return 1;
    }

    sim->G = 1.0;
    /* REBOUND ≥ 4.x switched the integrator selector from the
     * REB_INTEGRATOR_* enum to a string-based setter:
     *     reb_simulation_set_integrator(sim, "ias15");
     * Older REBOUND keeps the enum form:
     *     sim->integrator = REB_INTEGRATOR_IAS15;
     * Pick the right one for your installed version. */
    reb_simulation_set_integrator(sim, "ias15");

    /* Aarseth 1973 §I Pythagorean 3-4-5 initial conditions:
     *   Body 0: m = 3, r = ( 1,  3)
     *   Body 1: m = 4, r = (-2, -1)
     *   Body 2: m = 5, r = ( 1, -1)
     * Velocities all zero. Classical chaotic 3-body with multiple
     * close approaches in T = [0, 70]. */
    struct reb_particle p0 = {0}; p0.m = 3.0; p0.x =  1.0; p0.y =  3.0;
    struct reb_particle p1 = {0}; p1.m = 4.0; p1.x = -2.0; p1.y = -1.0;
    struct reb_particle p2 = {0}; p2.m = 5.0; p2.x =  1.0; p2.y = -1.0;
    reb_simulation_add(sim, p0);
    reb_simulation_add(sim, p1);
    reb_simulation_add(sim, p2);

    printf("# t, x0, y0, x1, y1, x2, y2\n");
    for (int k = 0; k <= 700; k++) {
        double t = (double)k * 0.1;
        reb_simulation_integrate(sim, t);
        /* %.17g is round-trip-safe: a CSV row reparsed by strtod
         * yields the exact same double. Necessary for cross-libc
         * gate determinism. */
        printf("%.6f, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g\n",
               t,
               sim->particles[0].x, sim->particles[0].y,
               sim->particles[1].x, sim->particles[1].y,
               sim->particles[2].x, sim->particles[2].y);
    }

    reb_simulation_free(sim);
    return 0;
}
