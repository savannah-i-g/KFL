/* k26astro_softkill/softkill.h — decoys + electronic warfare.
 *
 * libk26astro_softkill provides passive + active decoys, active
 * radar jamming, and chaff RCS statistics. The lib composes
 * against libk26astro_vehicle's PAYLOAD slot via the shared
 * K26AstroPayload base; each concrete (decoy, jammer) is its own
 * payload kind.
 *
 * Stated limitations:
 *   - Discrimination probabilities are tabulated values from
 *     open-literature missile-defence analyses (Garwin / Postol
 *     UCS reports). Application to the deep-space engagement-bubble
 *     regime is by analogy, not direct measurement.
 *   - Chaff statistics use the dipole-strip χ²(2) model; large
 *     clouds with significant multiple-scattering are not modelled.
 *   - Cyber-EW (intrusion against opponent avionics) is out of
 *     scope.
 *
 * Determinism: deterministic surfaces (J/S, burn-through, mean
 * RCS, CDF / quantile) are bit-identical across runs. Stochastic
 * surfaces (chaff sample draws) thread a K26CRng* explicitly;
 * passing NULL returns the deterministic mean.
 *
 * References:
 *   Skolnik 2008          Radar Handbook, 3rd ed.
 *   Knott et al. 2004     Radar Cross Section, 2nd ed.
 *   Garwin 2000 UCS       Countermeasures (PenAid evaluation)
 *   Postol 1991           Lessons of the Gulf War (Patriot)
 */
#ifndef K26ASTRO_SOFTKILL_SOFTKILL_H
#define K26ASTRO_SOFTKILL_SOFTKILL_H

#include "k26astro_softkill/softkill_consts.h"
#include "k26astro_softkill/decoy.h"
#include "k26astro_softkill/jammer.h"
#include "k26astro_softkill/chaff.h"

#endif /* K26ASTRO_SOFTKILL_SOFTKILL_H */
