/* test_ias15_two_body.c — IAS15 on Sun-Earth two-body.
 *
 * Strategy: IAS15 is a full N-body integrator — the Sun moves
 * under Earth's pull. The total system energy in the inertial
 * frame should be conserved to machine precision per step (IAS15
 * isn't symplectic but achieves 15th-order accuracy at each step).
 *
 * We DON'T compare against k26astro_kepler_propagate's fixed-Sun
 * trajectory because that's a different problem — heliocentric
 * Kepler vs N-body IAS15 differ by the Sun's barycentric motion
 * (~6 km over 100 days at Earth-Sun mass ratio).
 *
 * Acceptance:
 *   - Predictor-corrector iteration converges (rc=0 every step)
 *   - Total system (Sun+Earth) energy in inertial frame conserved
 *     to |dE/E| < 1e-10
 *   - Total angular momentum conserved to |dL/L| < 1e-10
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
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

    double m_S = K26A_GM_SUN  / K26A_G;
    double m_E = K26A_GM_EARTH / K26A_G;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15) == 0);

    /* Initial energy (full barycentric: both bodies move). */
    K26V3 r0 = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
    double r0_mag = sqrt(r0.x*r0.x + r0.y*r0.y + r0.z*r0.z);
    double E0 = 0.5 * m_S * (bodies[0].vel.x*bodies[0].vel.x
                            + bodies[0].vel.y*bodies[0].vel.y
                            + bodies[0].vel.z*bodies[0].vel.z)
              + 0.5 * m_E * (bodies[1].vel.x*bodies[1].vel.x
                            + bodies[1].vel.y*bodies[1].vel.y
                            + bodies[1].vel.z*bodies[1].vel.z)
              - K26A_G * m_S * m_E / r0_mag;
    /* Initial angular momentum (z-component since orbit is in xy-plane). */
    K26V3 p_sun_0   = k26astro_pos_to_m_approx(&bodies[0].pos);
    K26V3 p_earth_0 = k26astro_pos_to_m_approx(&bodies[1].pos);
    double L0 = m_S * (p_sun_0.x   * bodies[0].vel.y - p_sun_0.y   * bodies[0].vel.x)
              + m_E * (p_earth_0.x * bodies[1].vel.y - p_earth_0.y * bodies[1].vel.x);

    /* 100 steps at 1 day. */
    extern int    k26_ias15_last_iters;
    extern double k26_ias15_last_resid;
    double dt = 86400.0;
    int n_steps = 100;
    double dE_max = 0.0;
    int convergence_failures = 0;
    int max_iters_seen = 0;
    double max_resid_seen = 0.0;
    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        if (k26_ias15_last_iters > max_iters_seen)
            max_iters_seen = k26_ias15_last_iters;
        if (k26_ias15_last_resid > max_resid_seen)
            max_resid_seen = k26_ias15_last_resid;
        if (rc != 0) {
            convergence_failures++;
            if (convergence_failures > 0) {
                fprintf(stderr, "step %d rc=%d iters=%d resid=%.3e\n",
                        s, rc, k26_ias15_last_iters, k26_ias15_last_resid);
                break;
            }
        }
        K26V3 r = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
        double rm = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        double E = 0.5 * m_S * (bodies[0].vel.x*bodies[0].vel.x
                              + bodies[0].vel.y*bodies[0].vel.y
                              + bodies[0].vel.z*bodies[0].vel.z)
                 + 0.5 * m_E * (bodies[1].vel.x*bodies[1].vel.x
                              + bodies[1].vel.y*bodies[1].vel.y
                              + bodies[1].vel.z*bodies[1].vel.z)
                 - K26A_G * m_S * m_E / rm;
        double dE = fabs((E - E0) / E0);
        if (dE > dE_max) dE_max = dE;
    }

    assert(convergence_failures == 0);

    /* Final angular momentum drift. */
    K26V3 p_sun_f   = k26astro_pos_to_m_approx(&bodies[0].pos);
    K26V3 p_earth_f = k26astro_pos_to_m_approx(&bodies[1].pos);
    double L_final = m_S * (p_sun_f.x   * bodies[0].vel.y - p_sun_f.y   * bodies[0].vel.x)
                   + m_E * (p_earth_f.x * bodies[1].vel.y - p_earth_f.y * bodies[1].vel.x);
    double dL = fabs((L_final - L0) / L0);

    fprintf(stderr, "test_ias15_two_body: dE_max=%.3e  dL=%.3e (over %d steps; PC max iters=%d max resid=%.3e)\n",
            dE_max, dL, n_steps, max_iters_seen, max_resid_seen);
    /* Conservation tests: IAS15 should hold both energy and
     * angular momentum to high precision over 100 days. */
    assert(dE_max < 1.0e-10);
    assert(dL < 1.0e-10);

    k26astro_grav_state_destroy(&state);

    printf("test_ias15_two_body: OK\n");
    return 0;
}
