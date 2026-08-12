/* softkill_consts.h — decoy + electronic-warfare defaults.
 *
 * Anchored to published missile-defence countermeasure analyses
 * (Garwin / Bethe / Postol UCS reports) and the Knott et al.
 * Radar Cross Section textbook for chaff statistics.
 *
 * References:
 *   Sessler et al. 2000   UCS Countermeasures (penetration aids)
 *   Postol 1991           Lessons of the Gulf War (Patriot)
 *   Knott et al. 2004     Radar Cross Section, 2nd ed. (chaff)
 *   Skolnik 2008          Radar Handbook, 3rd ed. (jamming)
 */
#ifndef K26ASTRO_SOFTKILL_CONSTS_H
#define K26ASTRO_SOFTKILL_CONSTS_H

/* Decoy modality. */
typedef enum {
    K26ASTRO_DECOY_MODE_PASSIVE = 1,   /* mylar balloon, cold */
    K26ASTRO_DECOY_MODE_ACTIVE  = 2    /* heated + thrusting */
} K26AstroDecoyMode;

/* Observer-discriminator regime. A decoy that fools one regime
 * may not fool another. Bit-field so callers can combine. */
typedef enum {
    K26ASTRO_DISC_IR_ONLY         = 0x01,
    K26ASTRO_DISC_IR_PLUS_RCS     = 0x03,
    K26ASTRO_DISC_IR_RCS_ACCEL    = 0x07
} K26AstroDiscriminatorRegime;

/* Default RCS for a single dipole-strip chaff element at X-band,
 * m^2. Knott et al. 2004 §6.3. Order-of-magnitude. */
#define K26ASTRO_SOFTKILL_DIPOLE_RCS_DEFAULT_M2  1.0e-3

/* Default jammer noise figure (linear, dimensionless). */
#define K26ASTRO_SOFTKILL_JAMMER_NOISE_FIGURE_DEFAULT  3.0

/* Discrimination success rates: probability that the indicated
 * discriminator regime correctly identifies a decoy as a decoy
 * (i.e. rejects it). Deep-space-regime values are order-of-magnitude
 * consistent with the missile-defence penetration-aids literature
 * (Sessler, Cornwall, Garwin, Postol et al., "Countermeasures: A
 * Technical Evaluation of the Operational Effectiveness of the
 * Planned US National Missile Defense System," UCS/MIT Security
 * Studies Program, April 2000, Appendix A). Direct measurement
 * in the deep-space regime is not available; analogy stretch is
 * documented at the decoy lib's stated-limits section.
 *
 * Aggregate per-regime values (retained for reporting + driver
 * defaults — these are the "if every discriminator runs at full
 * effectiveness against a zero-match decoy, P_disc =" anchors). */
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_ONLY            0.15
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_PLUS_RCS        0.60
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_RCS_ACCEL       0.95

#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_IR_ONLY             0.05
#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_IR_PLUS_RCS         0.20
#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_IR_RCS_ACCEL        0.70

/* Per-discriminator single-channel effectiveness. Derived from the
 * aggregate values above by inverting the canonical product-of-
 * complements composition:
 *
 *     P_aggregate = 1 - ∏_i (1 - p_i)
 *
 * starting from p_ir = P_IR_ONLY and chaining:
 *
 *     p_rcs   = 1 - (1 - P_IR_PLUS_RCS)    / (1 - p_ir)
 *     p_accel = 1 - (1 - P_IR_RCS_ACCEL)   / ((1 - p_ir)(1 - p_rcs))
 *
 * This preserves the aggregate-table anchors while exposing
 * per-discriminator effectiveness for the canonical composition
 * used by the discrimination dispatch (which couples each p_i
 * with the decoy's per-discriminator match-quality m_i):
 *
 *     P_disc = 1 - ∏_i (1 - p_i × (1 - m_i)) */
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_SINGLE          0.15
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_RCS_SINGLE         0.5294117647058824 /* 1 - 0.40/0.85 */
#define K26ASTRO_SOFTKILL_P_DISC_PASSIVE_ACCEL_SINGLE       0.8750000000000000 /* 1 - 0.05/(0.85*0.4705882352941176) */

#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_IR_SINGLE           0.05
#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_RCS_SINGLE          0.1578947368421053 /* 1 - 0.80/0.95 */
#define K26ASTRO_SOFTKILL_P_DISC_ACTIVE_ACCEL_SINGLE        0.6249999999999999 /* 1 - 0.30/(0.95*0.8421052631578947) */

#endif /* K26ASTRO_SOFTKILL_CONSTS_H */
