/* test_perturb_gr_ppn1.c — Schwarzschild GR PPN-1 correction.
 *
 * Acceptance:
 *   - Returns non-zero acceleration for Mercury-around-Sun geometry
 *   - The correction is small (~10^-8 of Newtonian gravity for Mercury)
 *   - Sign: at periapsis (r̂·v = 0), the radial term dominates →
 *     accel is along -r̂ (perihelion shift forward over time)
 *   - Result is finite + non-NaN for typical inner-SS geometries
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_grav/perturb.h"
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

    /* Sun + Mercury at perihelion. */
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    /* Mercury: a = 0.387 AU, e = 0.206, periapsis = 0.307 AU. */
    double a_merc = 0.387 * K26A_AU_M;
    double e_merc = 0.206;
    double rp = a_merc * (1.0 - e_merc);
    double vp = sqrt(K26A_GM_SUN * (2.0/rp - 1.0/a_merc));

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = 2.203e13;            /* Mercury GM */
    bodies[1].pos  = k26astro_pos_from_m(rp, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, vp, 0.0);     /* purely tangential at periapsis */
    bodies[1].parent_body_idx = 0;

    K26AstroGravView view = k26astro_grav_view(bodies, 2);
    K26AstroGravState state;
    state.bodies = bodies;
    state.n_bodies = 2;

    K26V3 accel[2];
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_gr_ppn1(&state, &view, accel, NULL);

    /* Newtonian acceleration magnitude. */
    double a_newton = K26A_GM_SUN / (rp * rp);
    double a_gr = sqrt(accel[1].x*accel[1].x
                     + accel[1].y*accel[1].y
                     + accel[1].z*accel[1].z);
    double ratio = a_gr / a_newton;
    fprintf(stderr, "test_perturb_gr_ppn1: Newton=%.3e GR=%.3e ratio=%.3e\n",
            a_newton, a_gr, ratio);

    /* For Mercury, the GR correction is ~10⁻⁸ of Newtonian (this is
     * the famous 43"/century perihelion shift). */
    assert(ratio > 1e-9 && ratio < 1e-7);

    /* At periapsis (r̂·v = 0), the GR correction is purely radial:
     * the (r̂·v)·v term vanishes. So accel should be along ±x. */
    assert(fabs(accel[1].y) / a_gr < 0.01);  /* tangential negligible */
    /* Radial direction sign: prefactor = GM/(c²r²) > 0. radial term =
     * 4GM/r - v² > 0 for bound orbit at perihelion (v² < 2·GM/r since
     * bound). So accel.x should be positive (along +r̂). */
    assert(accel[1].x > 0.0);

    /* Finite + non-NaN. */
    assert(isfinite(accel[1].x));
    assert(isfinite(accel[1].y));
    assert(isfinite(accel[1].z));

    /* If we set GM_central to 0, the GR correction should vanish. */
    bodies[0].gm = 0.0;
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_gr_ppn1(&state, &view, accel, NULL);
    assert(accel[1].x == 0.0 && accel[1].y == 0.0 && accel[1].z == 0.0);

    printf("test_perturb_gr_ppn1: OK\n");
    return 0;
}
