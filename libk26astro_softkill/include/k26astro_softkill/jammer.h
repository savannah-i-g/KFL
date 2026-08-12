/* jammer.h — active radar jamming.
 *
 * The radar-EW equation governs the jamming-to-signal ratio at
 * the victim receiver:
 *
 *     J/S = (P_j · G_j · 4π · R_target²) / (P_t · G_t · σ)
 *
 * — Skolnik 2008 §15.3. The jammer is co-located with the host
 * or operates from a stand-off platform; the jamming power
 * arrives at the victim receiver at 1/R² (one-way) whereas the
 * radar return arrives at 1/R⁴ (two-way), so J/S grows as R².
 * The jammer wins at long range and is overcome at the
 * burn-through range R_BT given by setting J/S = threshold:
 *
 *     R_BT = √(SNR_th · P_t · G_t · σ / (P_j · G_j · 4π))
 *
 * Beyond R_BT the jammer dominates; within R_BT the radar burns
 * through.
 *
 * Counter-detection: the jammer announces its own position with
 * P_j watts of RF, detectable by passive RF sensors out to a
 * range far exceeding R_BT. Callers consume this via the
 * `_self_signature_W` helper alongside libk26astro_detect's
 * counter-detection range.
 *
 * Skolnik 2008 + Postol 1991 anchors.
 */
#ifndef K26ASTRO_SOFTKILL_JAMMER_H
#define K26ASTRO_SOFTKILL_JAMMER_H

#include "k26astro_softkill/softkill_consts.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    K26ASTRO_JAMMER_MODE_NOISE      = 1,   /* broadband noise */
    K26ASTRO_JAMMER_MODE_COVER_PULSE = 2,  /* narrowband cover pulse */
    K26ASTRO_JAMMER_MODE_DECEPTION  = 3    /* false-target deception */
} K26AstroJammerMode;

typedef struct K26AstroJammer K26AstroJammer;

K26AstroJammer *k26astro_jammer_new(
    K26AstroJammerMode mode,
    double             p_j_W,         /* jamming transmit power */
    double             g_j_db,        /* antenna gain */
    double             freq_Hz,       /* operating frequency */
    double             bandwidth_Hz,
    double             snr_threshold);

void k26astro_jammer_destroy(K26AstroJammer *j);
int  k26astro_jammer_bind(K26AstroJammer *j, struct K26AstroVehicle *v);
struct K26AstroPayload *k26astro_jammer_as_payload(K26AstroJammer *j);

K26AstroJammerMode k26astro_jammer_mode(const K26AstroJammer *j);

/* Compute J/S ratio at the victim radar receiver:
 *
 *   J/S = (P_j · G_j · 4π · R²) / (P_t · G_t · σ)
 *
 *   p_r_W            victim radar transmit power
 *   g_r_t_db         victim radar transmit antenna gain
 *   rcs_target_m2   victim's RCS (the jammer host's RCS,
 *                    typically — self-protection jamming)
 *   range_m          jammer-to-victim range */
double k26astro_jammer_js_ratio(const K26AstroJammer *j,
                                 double p_r_W,
                                 double g_r_t_db,
                                 double rcs_target_m2,
                                 double range_m);

/* Burn-through range (m) — the range at which J/S = snr_threshold.
 * Beyond this range the jammer dominates; within it, the radar
 * burns through and detects the target. */
double k26astro_jammer_burn_through_range(const K26AstroJammer *j,
                                           double p_r_W,
                                           double g_r_t_db,
                                           double rcs_target_m2);

/* Counter-detection: jammer's own emitted power (the value that
 * any passive RF observer can detect). Equal to p_j_W. */
double k26astro_jammer_self_signature_W(const K26AstroJammer *j);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_SOFTKILL_JAMMER_H */
