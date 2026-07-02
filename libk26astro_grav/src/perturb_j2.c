/* perturb_j2.c — J2 oblateness perturbation.
 *
 * For a body with non-zero J2 (second zonal harmonic), the
 * gravitational potential has an extra term that produces a small
 * non-radial force on orbiting bodies. The Sun's J2 is ~2e-7
 * (negligible at 1 AU but matters for Mercury); Earth's J2 is
 * 1.083e-3 (the dominant non-spherical effect on LEO satellites).
 *
 * Formula (potential expanded as spherical harmonics, gradient
 * gives acceleration; Vallado §9.4):
 *
 *   U_J2 = -(GM/r) · J2 · (R_eq/r)² · (3z²/r² - 1)/2
 *
 *   ∂U_J2/∂x = (3 GM J2 R_eq² / 2 r⁵) · x · (5 z²/r² - 1)
 *   ∂U_J2/∂y = (3 GM J2 R_eq² / 2 r⁵) · y · (5 z²/r² - 1)
 *   ∂U_J2/∂z = (3 GM J2 R_eq² / 2 r⁵) · z · (5 z²/r² - 3)
 *
 *   a_J2 = -∇U_J2
 *
 * For v0.1 we apply J2 from body[0] (the central body) to every
 * other body. The non-spherical body's rotation axis is assumed
 * aligned with +ẑ in the integration frame; future K26 versions
 * may transform via the body's attitude quaternion.
 *
 * Reference: Vallado (2013) Fundamentals of Astrodynamics §9.4. */
#include "k26astro_grav/perturb.h"
#include "k26astro_grav/grav.h"
#include "k26astro_core/pos.h"

#include <math.h>

void k26astro_perturb_j2(const K26AstroGravState *state,
                         const K26AstroGravView  *view,
                         K26V3 *accel_out, void *ctx)
{
    (void)ctx;
    if (!state || !view || !accel_out) return;
    if (view->n < 2) return;

    const K26AstroBody *central = &view->bodies[0];
    double J2 = central->j2;
    if (J2 == 0.0) return;

    double R_eq = central->radius;
    if (R_eq <= 0.0) return;
    double GM = central->gm;
    if (GM <= 0.0) return;

    double R_eq2 = R_eq * R_eq;

    for (int i = 1; i < view->n; i++) {
        const K26AstroBody *bi = &view->bodies[i];
        K26V3 r = k26astro_pos_sub(&bi->pos, &central->pos);
        double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
        if (r2 <= 0.0) continue;
        double rmag = sqrt(r2);
        double r5 = r2 * r2 * rmag;
        double inv_r5 = 1.0 / r5;
        double z2_over_r2 = (r.z * r.z) / r2;

        double common = (3.0 / 2.0) * GM * J2 * R_eq2 * inv_r5;
        K26V3 a_j2 = {
            common * r.x * (5.0 * z2_over_r2 - 1.0),
            common * r.y * (5.0 * z2_over_r2 - 1.0),
            common * r.z * (5.0 * z2_over_r2 - 3.0)
        };
        accel_out[i].x += a_j2.x;
        accel_out[i].y += a_j2.y;
        accel_out[i].z += a_j2.z;
        /* Newton's 3rd: equal-and-opposite on the central body,
         * scaled by mass ratio. */
        double mass_ratio = bi->gm / GM;
        accel_out[0].x -= a_j2.x * mass_ratio;
        accel_out[0].y -= a_j2.y * mass_ratio;
        accel_out[0].z -= a_j2.z * mass_ratio;
    }
}
