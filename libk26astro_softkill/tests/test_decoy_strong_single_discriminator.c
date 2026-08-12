/* test_decoy_strong_single_discriminator.c — strong-single-
 * discriminator gate.
 *
 * An averaging composition (mean of the three per-channel match
 * qualities applied to the aggregate regime probability as
 * (1 − m_avg)) under-weights a single strong discriminator
 * — a decoy whose IR signature is badly mismatched but whose RCS
 * and accel-history are perfectly matched should be caught at
 * approximately the IR single-channel rate (≈ 15% for passive
 * decoys), not derated by the average match (which would bring it
 * down to ~5%).
 *
 * This test asserts the canonical product-of-complements
 * composition:
 *
 *     P_disc = 1 - ∏_i (1 - p_i × (1 - m_i))
 *
 * gives the right behaviour in the strong-single-discriminator
 * case, while showing that the averaging form would produce a
 * materially different (wrong) answer for the same inputs. */
#include "k26astro_softkill/decoy.h"
#include "k26astro_softkill/softkill_consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx_(double got, double want, double rel_tol)
{
    const double denom = (fabs(want) > 1.0) ? fabs(want) : 1.0;
    return fabs(got - want) / denom <= rel_tol;
}

int main(void)
{
    /* Strong-single-discriminator scenario: passive decoy with
     * m_ir = 0.0 (terrible IR match — observer's IR channel
     * sees through it), m_rcs = 1.0 (perfect RCS match), and
     * m_accel = 1.0 (perfect acceleration-history match). */
    K26AstroDecoy *d_strong_single = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE,
        /* dry_mass_kg   */ 10.0,
        /* deploy_dv_mps */ 1.0,
        /* m_ir          */ 0.0,
        /* m_rcs         */ 1.0,
        /* m_accel       */ 1.0);
    assert(d_strong_single != NULL);

    /* Canonical composition for IR+RCS+ACCEL regime:
     *
     *   p_ir = 0.15, p_rcs = 0.529, p_accel = 0.875 (passive)
     *
     *   P_disc = 1 - (1 - 0.15×(1-0)) × (1 - 0.529×(1-1)) × (1 - 0.875×(1-1))
     *          = 1 - (1 - 0.15) × (1 - 0) × (1 - 0)
     *          = 1 - 0.85
     *          = 0.15
     *
     * which is exactly the IR single-channel effectiveness. The
     * RCS and accel channels contribute zero (perfect match means
     * decoy looks indistinguishable on those channels). */
    const double p_canonical = k26astro_decoy_probability_discriminated(
        d_strong_single, K26ASTRO_DISC_IR_RCS_ACCEL);
    assert(approx_(p_canonical,
                   K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_SINGLE,
                   1.0e-12));
    fprintf(stderr,
            "test_decoy_strong_single_discriminator: "
            "canonical P_disc = %.6f (= p_ir single-channel %.6f)\n",
            p_canonical,
            (double)K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_SINGLE);

    /* An averaging composition would compute:
     *
     *   m_avg = (0.0 + 1.0 + 1.0) / 3 = 0.667
     *   P_disc_avg = P_PASSIVE_IR_RCS_ACCEL × (1 - 0.667)
     *              = 0.95 × 0.333
     *              ≈ 0.317
     *
     * The canonical answer is 0.15, not 0.317 — the averaging
     * form double-counts the RCS/accel perfect matches and
     * effectively suppresses the IR channel's correct discrimination. */
    const double p_averaging = 0.95 * (1.0 - (0.0 + 1.0 + 1.0) / 3.0);
    assert(fabs(p_averaging - p_canonical) > 0.10);
    fprintf(stderr,
            "test_decoy_strong_single_discriminator: "
            "averaging composition would give %.6f — "
            "canonical is %.6f\n",
            p_averaging, p_canonical);

    /* ---- Symmetric case: strong RCS, weak others ------------- *
     *
     * Decoy that defeats IR + accel but exposes itself on RCS
     * (m_ir=1, m_rcs=0, m_accel=1). The IR-only regime should
     * give 0 (decoy fools IR perfectly), the IR+RCS regime
     * should reveal it at the RCS single-channel rate.
     *
     *   p_ir + (1-p_ir) × p_rcs × (1 - 0) − cross-term
     *   = 1 - (1 - 0.15×0)(1 - 0.529×1)
     *   = 1 - 1.0 × 0.471
     *   = 0.529 */
    K26AstroDecoy *d_rcs_failure = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE,
        10.0, 1.0,
        /* m_ir */    1.0,
        /* m_rcs */   0.0,
        /* m_accel */ 1.0);
    const double p_ir_only =
        k26astro_decoy_probability_discriminated(
            d_rcs_failure, K26ASTRO_DISC_IR_ONLY);
    assert(approx_(p_ir_only, 0.0, 1.0e-12));
    const double p_ir_rcs =
        k26astro_decoy_probability_discriminated(
            d_rcs_failure, K26ASTRO_DISC_IR_PLUS_RCS);
    assert(approx_(p_ir_rcs,
                   K26ASTRO_SOFTKILL_P_DISC_PASSIVE_RCS_SINGLE,
                   1.0e-12));

    /* ---- Backward-compat: aggregate values at m=0 ------------ *
     *
     * Zero-match-on-all-channels (worst-case decoy) must give the
     * documented aggregate-regime values, since the per-channel
     * single constants were derived to invert exactly that
     * composition. */
    K26AstroDecoy *d_zero = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE,
        10.0, 1.0,
        0.0, 0.0, 0.0);
    assert(approx_(
        k26astro_decoy_probability_discriminated(
            d_zero, K26ASTRO_DISC_IR_ONLY),
        K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_ONLY,
        1.0e-12));
    assert(approx_(
        k26astro_decoy_probability_discriminated(
            d_zero, K26ASTRO_DISC_IR_PLUS_RCS),
        K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_PLUS_RCS,
        1.0e-12));
    assert(approx_(
        k26astro_decoy_probability_discriminated(
            d_zero, K26ASTRO_DISC_IR_RCS_ACCEL),
        K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_RCS_ACCEL,
        1.0e-12));

    /* ---- Perfect-match (m_all = 1) drives P_disc to 0 -------- */
    K26AstroDecoy *d_perfect = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_ACTIVE,
        50.0, 1.0,
        1.0, 1.0, 1.0);
    assert(k26astro_decoy_probability_discriminated(
        d_perfect, K26ASTRO_DISC_IR_RCS_ACCEL) == 0.0);

    k26astro_decoy_destroy(d_strong_single);
    k26astro_decoy_destroy(d_rcs_failure);
    k26astro_decoy_destroy(d_zero);
    k26astro_decoy_destroy(d_perfect);

    fprintf(stderr,
            "test_decoy_strong_single_discriminator: OK\n");
    return 0;
}
