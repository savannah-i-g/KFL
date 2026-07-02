/* soi_crossing.c — SOI boundary crossing detection by linear
 * interpolation of distance.
 *
 * Given r₀ = |pos_before| and r₁ = |pos_after| relative to the
 * candidate parent, plus the parent's SOI radius R, return:
 *   1  if (r₀ - R) and (r₁ - R) have opposite signs (a crossing)
 *   0  if both stayed on the same side
 *  -1  if input is malformed
 *
 * The linearly-interpolated crossing fraction τ ∈ [0, 1] is:
 *   τ = (R - r₀) / (r₁ - r₀)
 * This is exact for a body on a radial trajectory; for tangential
 * motion the actual crossing happens between this estimate and the
 * point of closest approach. For the patched-conic handoff at typical
 * orbital scales (SOI is a few % of orbital radius), the linear-r
 * approximation is sufficient as the handoff trigger; the MERCURIUS
 * orchestration in libk26astro_rt refines via sub-step rewind.
 *
 * Tangential-only-crossing (the body grazes the SOI without entering)
 * is NOT detected by this function; it requires comparing minimum-r
 * to R, which involves the orbital geometry. The miss is intentional
 * at this layer; full motion-aware crossing detection lives in
 * libk26astro_grav's close_encounter.c. */
#include "k26astro_conics/soi.h"

#include <math.h>

int k26astro_soi_crossing_detect(K26V3 pos_before_rel_parent,
                                  K26V3 pos_after_rel_parent,
                                  double soi_radius,
                                  double *out_t_cross)
{
    if (soi_radius <= 0.0) return -1;

    double r0 = sqrt(pos_before_rel_parent.x * pos_before_rel_parent.x
                   + pos_before_rel_parent.y * pos_before_rel_parent.y
                   + pos_before_rel_parent.z * pos_before_rel_parent.z);
    double r1 = sqrt(pos_after_rel_parent.x * pos_after_rel_parent.x
                   + pos_after_rel_parent.y * pos_after_rel_parent.y
                   + pos_after_rel_parent.z * pos_after_rel_parent.z);

    double d0 = r0 - soi_radius;
    double d1 = r1 - soi_radius;

    if ((d0 >= 0.0 && d1 >= 0.0) || (d0 <= 0.0 && d1 <= 0.0)) {
        if (out_t_cross) *out_t_cross = 0.0;
        return 0;
    }

    if (out_t_cross) {
        /* Linear-r interpolation. */
        double denom = r1 - r0;
        if (fabs(denom) < 1e-30) {
            *out_t_cross = 0.5;
        } else {
            double tau = (soi_radius - r0) / denom;
            if (tau < 0.0) tau = 0.0;
            if (tau > 1.0) tau = 1.0;
            *out_t_cross = tau;
        }
    }
    return 1;
}
