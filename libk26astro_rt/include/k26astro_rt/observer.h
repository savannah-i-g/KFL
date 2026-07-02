/* libk26astro_rt — observer corrections.
 *
 * For in-world bodies the runtime implements its own light-time
 * iteration directly (using K26AstroBody.vel for linear-extrap
 * retarded position); for ephemeris-backed targets the iteration
 * delegates to libk26astro_ephem's k26astro_ephem_observe.
 * Adds stellar aberration (Aoki et al. 1983 first-order form) and
 * GR Shapiro delay (Will 1993 eq. 8.36) when the world has
 * gr_ppn1 enabled. */
#ifndef K26ASTRO_RT_OBSERVER_H
#define K26ASTRO_RT_OBSERVER_H

#include "k26astro_rt/world.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Compute an observation of `target_idx` as seen from `observer_idx`
 * at the world's current epoch. The world's current observer mode
 * (set via k26astro_world_set_observer_mode) selects the correction
 * tier:
 *   - GEOMETRIC: instantaneous straight-line vector
 *   - ASTROMETRIC: light-time-corrected emission point
 *   - APPARENT: + stellar aberration + Shapiro (if gr_ppn1 enabled)
 *   - TOPOCENTRIC: not implemented; returns E_NOT_IMPLEMENTED
 *
 * `out_target_pos` receives the position of the target in the same
 * frame as the observer's pos (sector-grid). `out_apparent_dir` is
 * the unit vector from observer to target as it actually appears
 * to the observer (post-correction). Either may be NULL. */
int k26astro_world_observe(const K26AstroWorld *world,
                            int target_idx, int observer_idx,
                            K26AstroPos *out_target_pos,
                            K26V3 *out_apparent_dir);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_OBSERVER_H */
