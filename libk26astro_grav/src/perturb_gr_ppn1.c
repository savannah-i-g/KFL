/* perturb_gr_ppn1.c — General Relativity PPN-1 correction.
 *
 * The Schwarzschild post-Newtonian-1 correction to point-mass gravity
 * in the Einstein-Infeld-Hoffmann (EIH) form:
 *
 *   a_GR_on_i = (GM_c / (c² · r²)) · [(4·GM_c/r - v²) · r̂ + 4·(r̂·v)·v]
 *
 * where GM_c is the central body's GM, r is the body's heliocentric
 * distance, v is the heliocentric velocity, and c is the speed of
 * light.
 *
 * For Mercury at 0.4 AU, this term produces the famous 43"/century
 * perihelion advance (Einstein 1915, confirmed by observation).
 *
 * v0.1 limitation: only the central body's gravitational well
 * generates the GR correction. The full N-body EIH form would
 * include cross-pair GR terms; those matter only for tight binary
 * pulsars and aren't relevant at solar-system inner-planet scales.
 *
 * Reference: Will (1993) "Theory and Experiment in Gravitational
 * Physics", §8.3 Eq. 8.62. */
#include "k26astro_grav/perturb.h"
#include "k26astro_grav/grav.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <math.h>

void k26astro_perturb_gr_ppn1(const K26AstroGravState *state,
                              const K26AstroGravView  *view,
                              K26V3 *accel_out, void *ctx)
{
    (void)ctx;
    if (!state || !view || !accel_out) return;
    if (view->n < 2) return;

    const K26AstroBody *central = &view->bodies[0];
    double GM_c = central->gm;
    if (GM_c <= 0.0) return;

    const double c2 = K26A_C * K26A_C;

    for (int i = 1; i < view->n; i++) {
        const K26AstroBody *bi = &view->bodies[i];
        K26V3 r = k26astro_pos_sub(&bi->pos, &central->pos);
        K26V3 v_rel = {
            bi->vel.x - central->vel.x,
            bi->vel.y - central->vel.y,
            bi->vel.z - central->vel.z
        };

        double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
        if (r2 <= 0.0) continue;
        double rmag = sqrt(r2);
        K26V3 r_hat = { r.x/rmag, r.y/rmag, r.z/rmag };

        double v2 = v_rel.x*v_rel.x + v_rel.y*v_rel.y + v_rel.z*v_rel.z;
        double r_dot_v = r_hat.x*v_rel.x + r_hat.y*v_rel.y + r_hat.z*v_rel.z;

        /* a = (GM_c / (c²·r²)) · [(4·GM_c/r - v²)·r̂ + 4·(r̂·v)·v] */
        double prefactor = GM_c / (c2 * r2);
        double radial    = 4.0 * GM_c / rmag - v2;
        double velocity_term = 4.0 * r_dot_v;

        K26V3 a_gr = {
            prefactor * (radial * r_hat.x + velocity_term * v_rel.x),
            prefactor * (radial * r_hat.y + velocity_term * v_rel.y),
            prefactor * (radial * r_hat.z + velocity_term * v_rel.z)
        };
        accel_out[i].x += a_gr.x;
        accel_out[i].y += a_gr.y;
        accel_out[i].z += a_gr.z;
    }
}
