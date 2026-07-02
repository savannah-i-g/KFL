/* test_wh_two_body.c — WH on Sun-Earth two-body.
 *
 * For exactly two bodies (central + planet), WH's interaction kick is
 * zero (no other planets to perturb). The KDK step reduces to pure
 * Kepler drift, which is exact for the universal-variable propagator.
 * So WH(2-body) at any dt = exact-Kepler at that dt.
 *
 * Acceptance:
 *   - After 100 steps of dt=1 day, Earth's position matches what
 *     k26astro_kepler_propagate produces in one big step of 100 days.
 *   - Semi-major axis preserved to <1e-9 fractional drift.
 *   - Orbit stays in the xy-plane (|z| < 1 m).
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/wisdom_holman.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* Earth-like circular orbit at 1 AU. */
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
    double v_circ  = sqrt(K26A_GM_SUN / K26A_AU_M);
    bodies[1].vel  = k26m3d_v3(0.0, v_circ, 0.0);
    bodies[1].parent_body_idx = 0;

    /* Stash initial state for the Kepler-comparison run. */
    K26V3 r0 = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
    K26V3 v0 = bodies[1].vel;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_WH) == 0);

    /* Run WH 100 steps at 1 day. */
    double dt = 86400.0;
    int n_steps = 100;
    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        assert(rc == 0);
    }

    /* Compare with one big Kepler step over 100 days. */
    K26V3 r_kepler, v_kepler;
    int rc = k26astro_kepler_propagate(&r_kepler, &v_kepler,
                                        r0, v0, K26A_GM_SUN,
                                        n_steps * dt, 128);
    assert(rc == 0);

    K26V3 r_wh = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
    double dx = r_wh.x - r_kepler.x;
    double dy = r_wh.y - r_kepler.y;
    double dz = r_wh.z - r_kepler.z;
    double miss = sqrt(dx*dx + dy*dy + dz*dz);
    printf("test_wh_two_body: WH vs Kepler miss=%.3e m (%d days)\n",
           miss, n_steps);
    /* For 2-body WH = exact Kepler, miss should be ~ULP-noise-level,
     * well under 1 mm at 1 AU. */
    assert(miss < 1.0e-2);

    /* Semi-major axis preservation. */
    double r_mag = sqrt(r_wh.x*r_wh.x + r_wh.y*r_wh.y + r_wh.z*r_wh.z);
    double v_mag = sqrt(bodies[1].vel.x*bodies[1].vel.x
                       + bodies[1].vel.y*bodies[1].vel.y
                       + bodies[1].vel.z*bodies[1].vel.z);
    double a_final = 1.0 / (2.0 / r_mag - v_mag * v_mag / K26A_GM_SUN);
    double a_drift = fabs(a_final - K26A_AU_M) / K26A_AU_M;
    printf("test_wh_two_body: a drift  = %.3e (relative)\n", a_drift);
    assert(a_drift < 1e-9);

    /* Orbital plane preserved. */
    assert(fabs(r_wh.z) < 1.0);
    assert(fabs(bodies[1].vel.z) < 1.0e-6);

    k26astro_grav_state_destroy(&state);

    printf("test_wh_two_body: OK\n");
    return 0;
}
