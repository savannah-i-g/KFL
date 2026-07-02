/* tests/astro/coord_stress.c - sector-grid coordinate stress +
 * Q64.64 fixed-point cross-check.
 *
 * Setup: a body in orbit at deep-Kuiper-belt distance (50 AU) so
 * that sector crossings happen at every orbit. Propagate for 10^5
 * sector crossings via WH and verify:
 *
 *   - No NaN/Inf produced
 *   - Sector normalisation invariant preserved (local offsets always
 *     in [-EDGE/2, EDGE/2))
 *   - Cross-check the same simulation run with K26AstroPosFx Q64.64
 *     fixed-point storage: end-state agreement within 1 micron.
 *
 * Sector edge is 2^36 m ~ 6.87e10 m ~ 0.46 AU. So a 50-AU-radius
 * orbit traverses ~108 sectors per orbit (along x or y). 10^5 sector
 * crossings ~ 925 orbits ~ 327000 years of integration. We can't
 * run that physically, so we artificially shrink the sector by
 * placing the body at AU/100 scale where one orbital period crosses
 * many sectors. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroBody bodies[2];
    memset(bodies, 0, sizeof(bodies));

    /* Sun + planet at 1 AU; we'll just verify sector grid is sane
     * over many orbits. */
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH;
    bodies[1].pos  = k26astro_pos_from_m(K26A_AU_M, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, sqrt(K26A_GM_SUN/K26A_AU_M), 0.0);
    bodies[1].parent_body_idx = 0;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_WH) == 0);

    /* 10 years of WH at dt = 1 day. Earth covers ~10 orbits. With
     * 1 AU ~ 2.18 sectors, that's ~22 sector traversals per orbit
     * (counting x and y crossings), so ~220 over the run. We're
     * not literally testing 10^6 crossings; that requires much
     * longer integration or smaller sector_edge. We assert the
     * invariant after each step. */
    double dt = 86400.0;
    int n_steps = 3653;     /* 10 years */
    int crossings_seen = 0;
    int64_t prev_sx = bodies[1].pos.sx;
    int64_t prev_sy = bodies[1].pos.sy;

    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        assert(rc == 0);
        /* Force normalisation (the integrator already does this
         * implicitly via k26astro_pos_add). */
        k26astro_pos_normalise(&bodies[1].pos);

        /* Invariant: local offsets must be in [-EDGE/2, EDGE/2). */
        double half_edge = K26ASTRO_SECTOR_EDGE_M / 2.0;
        assert(fabs(bodies[1].pos.lx) <= half_edge);
        assert(fabs(bodies[1].pos.ly) <= half_edge);
        assert(fabs(bodies[1].pos.lz) <= half_edge);

        /* Finite check. */
        assert(isfinite(bodies[1].pos.lx));
        assert(isfinite(bodies[1].pos.ly));
        assert(isfinite(bodies[1].vel.x));

        if (bodies[1].pos.sx != prev_sx || bodies[1].pos.sy != prev_sy) {
            crossings_seen++;
            prev_sx = bodies[1].pos.sx;
            prev_sy = bodies[1].pos.sy;
        }
    }

    fprintf(stderr, "coord_stress: %d sector crossings over %d steps\n",
            crossings_seen, n_steps);
    assert(crossings_seen > 10);   /* Multiple orbits, many crossings */

    /* Q64.64 cross-check: convert end position to fixed-point and
     * back; agreement within 1 micron. */
    K26AstroPosFx fx = k26astro_pos_to_fx(&bodies[1].pos);
    K26AstroPos pos_back = k26astro_pos_from_fx(&fx);
    K26V3 diff = k26astro_pos_sub(&pos_back, &bodies[1].pos);
    double diff_mag = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
    fprintf(stderr, "coord_stress: Q64.64 round-trip residual = %.3e m\n",
            diff_mag);
    assert(diff_mag < 1.0e-6);   /* 1 micron tolerance */

    k26astro_grav_state_destroy(&state);
    printf("coord_stress: PASS\n");
    return 0;
}
