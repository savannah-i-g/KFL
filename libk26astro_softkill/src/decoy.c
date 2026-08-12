/* decoy.c — decoy opaque + discrimination probability dispatch.
 *
 * The lib's decoy payload deploys a single decoy with caller-
 * specified match qualities for each discriminator dimension
 * (IR / RCS / accel-history). The discrimination probability is
 * the table value (per mode, per regime) multiplied by the
 * caller-supplied 1 - match_quality factor — a higher match
 * quality reduces the chance the observer's discriminator
 * succeeds. */

#include "k26astro_softkill/decoy.h"

#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#include <stdlib.h>
#include <string.h>

struct K26AstroDecoy {
    K26AstroPayload    base;
    K26AstroDecoyMode mode;
    double            dry_mass_kg;
    double            deploy_dv_mps;
    double            ir_match_quality;
    double            rcs_match_quality;
    double            accel_match_quality;
};

K26AstroDecoy *k26astro_decoy_new(
    K26AstroDecoyMode mode,
    double            dry_mass_kg,
    double            deploy_dv_mps,
    double            ir_match_quality,
    double            rcs_match_quality,
    double            accel_match_quality)
{
    if (mode != K26ASTRO_DECOY_MODE_PASSIVE &&
        mode != K26ASTRO_DECOY_MODE_ACTIVE) return NULL;
    if (dry_mass_kg <= 0.0 || deploy_dv_mps < 0.0) return NULL;

    if (ir_match_quality    < 0.0) ir_match_quality    = 0.0;
    if (ir_match_quality    > 1.0) ir_match_quality    = 1.0;
    if (rcs_match_quality   < 0.0) rcs_match_quality   = 0.0;
    if (rcs_match_quality   > 1.0) rcs_match_quality   = 1.0;
    if (accel_match_quality < 0.0) accel_match_quality = 0.0;
    if (accel_match_quality > 1.0) accel_match_quality = 1.0;

    K26AstroDecoy *d = (K26AstroDecoy *)calloc(1, sizeof(*d));
    if (!d) return NULL;
    d->base.kind             = K26ASTRO_DEFENSE_KIND_DECOY;
    d->base.owner            = NULL;
    d->base.owner_slot       = -1;
    d->base.generation       = 1;
    d->base.on_owner_destroy = NULL;
    d->mode                  = mode;
    d->dry_mass_kg           = dry_mass_kg;
    d->deploy_dv_mps         = deploy_dv_mps;
    d->ir_match_quality      = ir_match_quality;
    d->rcs_match_quality     = rcs_match_quality;
    d->accel_match_quality   = accel_match_quality;
    return d;
}

void k26astro_decoy_destroy(K26AstroDecoy *d)
{
    if (!d) return;
    k26astro_payload_unlink_(&d->base);
    free(d);
}

int k26astro_decoy_bind(K26AstroDecoy *d, struct K26AstroVehicle *v)
{
    if (!d) return -1;
    return k26astro_payload_bind(&d->base, v);
}

struct K26AstroPayload *k26astro_decoy_as_payload(K26AstroDecoy *d)
{
    return d ? &d->base : NULL;
}

K26AstroDecoyMode k26astro_decoy_mode(const K26AstroDecoy *d)
{
    return d ? d->mode : K26ASTRO_DECOY_MODE_PASSIVE;
}

/* Discrimination probability dispatch.
 *
 * Canonical "any active discriminator catches the decoy"
 * composition (Garwin / Sessler et al. 2000 UCS, §Appendix A;
 * the standard form across the penetration-aids literature):
 *
 *     P_disc = 1 - ∏_i (1 - p_i × (1 - m_i))
 *
 * where p_i is the per-discriminator effectiveness from the
 * softkill_consts.h SINGLE-channel table, and m_i is this
 * decoy's match-quality for that discriminator. A perfectly-
 * matched signature on any one channel (m_i = 1) zeros that
 * channel's contribution; a complete mismatch (m_i = 0) leaves
 * the channel running at p_i.
 *
 * The product form (rather than averaging the match qualities)
 * keeps a single strong discriminator decisive: a decoy that
 * defeats two of three channels but fails the third badly is
 * still caught at that channel's single rate, where an average
 * derate would suppress it. */
double k26astro_decoy_probability_discriminated(
    const K26AstroDecoy *d,
    K26AstroDiscriminatorRegime regime)
{
    if (!d) return 0.0;

    /* Per-discriminator single-channel effectiveness for this mode. */
    const double p_ir =
        (d->mode == K26ASTRO_DECOY_MODE_PASSIVE)
            ? K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_SINGLE
            : K26ASTRO_SOFTKILL_P_DISC_ACTIVE_IR_SINGLE;
    const double p_rcs =
        (d->mode == K26ASTRO_DECOY_MODE_PASSIVE)
            ? K26ASTRO_SOFTKILL_P_DISC_PASSIVE_RCS_SINGLE
            : K26ASTRO_SOFTKILL_P_DISC_ACTIVE_RCS_SINGLE;
    const double p_accel =
        (d->mode == K26ASTRO_DECOY_MODE_PASSIVE)
            ? K26ASTRO_SOFTKILL_P_DISC_PASSIVE_ACCEL_SINGLE
            : K26ASTRO_SOFTKILL_P_DISC_ACTIVE_ACCEL_SINGLE;

    /* Product of channel-miss probabilities (probability that
     * NONE of the active discriminators catches the decoy). The
     * regime bit-field selects which channels are active:
     *
     *   0x01 — IR active (always, per existing API)
     *   0x02 — RCS active
     *   0x04 — accel-history active */
    double prod_miss = 1.0;
    /* IR always runs in every regime. */
    prod_miss *= (1.0 - p_ir * (1.0 - d->ir_match_quality));
    if (regime & 0x02) {
        prod_miss *= (1.0 - p_rcs * (1.0 - d->rcs_match_quality));
    }
    if (regime & 0x04) {
        prod_miss *= (1.0 - p_accel * (1.0 - d->accel_match_quality));
    }
    if (prod_miss < 0.0) prod_miss = 0.0;
    if (prod_miss > 1.0) prod_miss = 1.0;
    return 1.0 - prod_miss;
}
