/* decoy.h — passive + active decoys.
 *
 * A decoy is deployed with a separation Δv from the host and
 * advertised to the world as the host's signature surrogate.
 * Effectiveness depends on which discriminators the observer
 * runs:
 *
 *   - IR-only:   passive decoys fail (mass-temperature mismatch);
 *                active decoys can match if their heater + thrust
 *                budget supports it.
 *   - IR+RCS:    structural / radar-cross-section discrimination
 *                catches simple balloon decoys; active decoys with
 *                appropriately-shaped radar reflectors survive.
 *   - IR+RCS+accel-history: trajectory discrimination — under
 *                differential solar pressure or unmasked thrusting,
 *                decoys deviate from the host's accel profile.
 *
 * The probability that a given decoy successfully fools an
 * observer running a given discriminator regime is the
 * complement of `_probability_discriminated`.
 *
 * Sessler et al. 2000 (UCS Countermeasures) + Postol 1991 anchors.
 */
#ifndef K26ASTRO_SOFTKILL_DECOY_H
#define K26ASTRO_SOFTKILL_DECOY_H

#include "k26astro_softkill/softkill_consts.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct K26AstroDecoy K26AstroDecoy;

K26AstroDecoy *k26astro_decoy_new(
    K26AstroDecoyMode mode,
    double            dry_mass_kg,
    double            deploy_dv_mps,
    double            ir_match_quality,        /* in [0, 1] */
    double            rcs_match_quality,       /* in [0, 1] */
    double            accel_match_quality);    /* in [0, 1] */

void k26astro_decoy_destroy(K26AstroDecoy *d);
int  k26astro_decoy_bind(K26AstroDecoy *d, struct K26AstroVehicle *v);
struct K26AstroPayload *k26astro_decoy_as_payload(K26AstroDecoy *d);
K26AstroDecoyMode k26astro_decoy_mode(const K26AstroDecoy *d);

/* Probability that a decoy is correctly identified (discriminated)
 * by an observer running the given discriminator regime. Returns 0
 * if `d` is NULL. */
double k26astro_decoy_probability_discriminated(
    const K26AstroDecoy *d,
    K26AstroDiscriminatorRegime regime);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_SOFTKILL_DECOY_H */
