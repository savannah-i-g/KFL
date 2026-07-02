/* test_rk45_via_compute.c — verify RK4 + RK45 wrappers via Sun-Earth.
 *
 * Symplectic-free integrators don't conserve energy exactly, but
 * Dormand-Prince at rtol=1e-9 / atol=1e-12 holds Earth in orbit for
 * 10 days within ~1000 km (barycentric energy drift ~1e-8). RK4
 * fixed-step at dt=1hr is much looser (~1e-3); we check both at
 * their appropriate tolerances.
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void setup_(K26AstroBody bodies[2])
{
    memset(bodies, 0, sizeof(K26AstroBody) * 2);
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
}

int main(void)
{
    /* ---- RK45 adaptive ------------------------------------- */
    {
        K26AstroBody bodies[2];
        setup_(bodies);
        K26AstroGravState state;
        assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
        assert(k26astro_grav_set_integrator(&state,
                                              K26ASTRO_INTEGRATOR_RK45) == 0);

        /* 10 days of integration. */
        int rc = k26astro_grav_step(&state, 10.0 * 86400.0);
        if (rc != 0) {
            fprintf(stderr, "RK45 step rc=%d\n", rc);
        }
        assert(rc == 0);

        /* Earth's distance from Sun should still be ~1 AU
         * (circular orbit ≠ exactly 1 AU due to the slight
         * eccentricity introduced by the Sun-Earth coupling, but
         * within ~1000 km). */
        K26V3 r = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
        double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        double off = fabs(r_mag - K26A_AU_M);
        printf("test_rk45_via_compute: RK45 |r|-AU=%.3e m\n", off);
        assert(off < 1.0e5);

        k26astro_grav_state_destroy(&state);
    }

    /* ---- RK4 fixed-step ------------------------------------ */
    {
        K26AstroBody bodies[2];
        setup_(bodies);
        K26AstroGravState state;
        assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
        assert(k26astro_grav_set_integrator(&state,
                                              K26ASTRO_INTEGRATOR_RK4) == 0);

        /* RK4 fixed-step: 10 days in one go is too coarse; libk26compute's
         * RK4 driver runs n_steps=1 internally per our wrapper. So
         * step 240 times at 1 hour each. */
        for (int hour = 0; hour < 240; hour++) {
            int rc = k26astro_grav_step(&state, 3600.0);
            assert(rc == 0);
        }

        K26V3 r = k26astro_pos_sub(&bodies[1].pos, &bodies[0].pos);
        double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        double off = fabs(r_mag - K26A_AU_M);
        printf("test_rk45_via_compute: RK4  |r|-AU=%.3e m\n", off);
        /* RK4 at dt=1hr over 10 days: loose bound. */
        assert(off < 1.0e6);

        k26astro_grav_state_destroy(&state);
    }

    printf("test_rk45_via_compute: OK\n");
    return 0;
}
