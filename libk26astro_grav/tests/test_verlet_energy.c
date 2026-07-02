/* test_verlet_energy.c — DKD Verlet symplectic integrator.
 *
 * Acceptance:
 *   - Sun + Earth circular orbit, 100 days at dt = 1 hour.
 *     Energy oscillation bounded: |ΔE/E_0| peak < 1e-6.
 *     Symplectic = no secular drift, just oscillation.
 *   - Verlet preserves the orbital plane (z-momentum stays at 0).
 *   - Integrator dispatch via k26astro_grav_step routes correctly.
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/verlet.h"
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
    /* Two-body system: Sun at origin, Earth circular at 1 AU. */
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

    K26AstroGravState state;
    int rc = k26astro_grav_state_init(&state, bodies, 2);
    assert(rc == 0);
    rc = k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_VERLET);
    assert(rc == 0);

    /* Barycentric total energy:
     *   E_tot = 0.5·(m_S·|v_S|² + m_E·|v_E|²) - G·m_S·m_E / r_SE
     * With G·m_S = μ_S (so m_S = μ_S/G) etc. */
    double m_S = K26A_GM_SUN  / K26A_G;
    double m_E = K26A_GM_EARTH / K26A_G;
    double G   = K26A_G;
    K26V3 r_rel = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
    double r_mag = sqrt(r_rel.x*r_rel.x + r_rel.y*r_rel.y + r_rel.z*r_rel.z);
    double E0 = 0.5 * m_S * (bodies[0].vel.x*bodies[0].vel.x
                            + bodies[0].vel.y*bodies[0].vel.y
                            + bodies[0].vel.z*bodies[0].vel.z)
              + 0.5 * m_E * (bodies[1].vel.x*bodies[1].vel.x
                            + bodies[1].vel.y*bodies[1].vel.y
                            + bodies[1].vel.z*bodies[1].vel.z)
              - G * m_S * m_E / r_mag;

    /* dt = 1 hour, 100 days. */
    double dt = 3600.0;
    int n_steps = 100 * 24;
    double dE_max = 0.0;

    for (int step = 0; step < n_steps; step++) {
        rc = k26astro_grav_step(&state, dt);
        assert(rc == 0);

        K26V3 r = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
        double rm = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        double E = 0.5 * m_S * (bodies[0].vel.x*bodies[0].vel.x
                              + bodies[0].vel.y*bodies[0].vel.y
                              + bodies[0].vel.z*bodies[0].vel.z)
                 + 0.5 * m_E * (bodies[1].vel.x*bodies[1].vel.x
                              + bodies[1].vel.y*bodies[1].vel.y
                              + bodies[1].vel.z*bodies[1].vel.z)
                 - G * m_S * m_E / rm;
        double dE = fabs((E - E0) / E0);
        if (dE > dE_max) dE_max = dE;
    }

    printf("test_verlet_energy: dE_max=%.3e (over %d steps)\n",
           dE_max, n_steps);
    /* DKD Verlet for Sun + Earth at dt=1hr: ε around 1e-9 to 1e-7. */
    assert(dE_max < 1.0e-6);

    /* Orbital plane preserved: z-position of Earth should stay at 0. */
    K26V3 r_final = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
    assert(fabs(r_final.z) < 1.0);   /* < 1 m out of plane */

    /* Body 1's velocity z-component stays 0 too. */
    assert(fabs(bodies[1].vel.z) < 1.0e-6);

    /* Force balance: the Sun's COM should not be moving (it has no
     * external forces; only Earth pulls it). Sun's velocity should
     * accumulate by a tiny amount over the run. */
    assert(fabs(bodies[0].vel.x) < 1.0e3);
    assert(fabs(bodies[0].vel.y) < 1.0e3);

    k26astro_grav_state_destroy(&state);

    printf("test_verlet_energy: OK\n");
    return 0;
}
