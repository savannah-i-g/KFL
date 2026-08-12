/* test_softkill_physics.c — closed-form gates for decoy
 * discrimination + jammer EW equation + chaff χ² statistics. */

#include "k26astro_softkill/softkill.h"
#include "k26astro_core/consts.h"
#include "k26compute.h"

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
    /* ---- Decoy discrimination probability -------------------- *
     *
     * Passive decoy with all match qualities = 0 (poor match)
     * gets caught at the rates from the constant table. */
    K26AstroDecoy *d_perf_bad = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE,
        10.0, 1.0,
        /* ir/rcs/accel match */ 0.0, 0.0, 0.0);
    const double p_ir = k26astro_decoy_probability_discriminated(
        d_perf_bad, K26ASTRO_DISC_IR_ONLY);
    const double p_rcs = k26astro_decoy_probability_discriminated(
        d_perf_bad, K26ASTRO_DISC_IR_PLUS_RCS);
    const double p_all = k26astro_decoy_probability_discriminated(
        d_perf_bad, K26ASTRO_DISC_IR_RCS_ACCEL);
    assert(approx_(p_ir,  K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_ONLY, 1.0e-12));
    assert(approx_(p_rcs, K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_PLUS_RCS, 1.0e-12));
    assert(approx_(p_all, K26ASTRO_SOFTKILL_P_DISC_PASSIVE_IR_RCS_ACCEL, 1.0e-12));
    /* Discrimination escalates: IR-only < IR+RCS < full. */
    assert(p_ir < p_rcs);
    assert(p_rcs < p_all);

    /* Active decoy of equal-quality match is harder to catch. */
    K26AstroDecoy *d_act_bad = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_ACTIVE,
        50.0, 1.0,
        0.0, 0.0, 0.0);
    const double p_act_ir = k26astro_decoy_probability_discriminated(
        d_act_bad, K26ASTRO_DISC_IR_ONLY);
    assert(p_act_ir < p_ir);

    /* Perfect match quality drives discrimination to zero. */
    K26AstroDecoy *d_perfect = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_ACTIVE,
        50.0, 1.0,
        1.0, 1.0, 1.0);
    assert(k26astro_decoy_probability_discriminated(
        d_perfect, K26ASTRO_DISC_IR_RCS_ACCEL) == 0.0);

    k26astro_decoy_destroy(d_perf_bad);
    k26astro_decoy_destroy(d_act_bad);
    k26astro_decoy_destroy(d_perfect);

    /* ---- Jammer J/S equation -------------------------------- *
     *
     *   J/S = (P_j G_j 4π R²) / (P_t G_t σ)             Skolnik §15.3
     *
     * Closed-form at:
     *   P_j = 1 kW, G_j = 25 dB (10^2.5 ≈ 316.23)
     *   P_t = 10 kW, G_t = 35 dB (10^3.5 ≈ 3162.28)
     *   σ   = 1 m²
     *   R   = 1000 km = 1e6 m
     * J/S = (1e3 × 316.23 × 12.566 × 1e12) / (1e4 × 3162.28 × 1)
     *      = 3.974e18 / 3.162e7 ≈ 1.2566e+11.
     *
     * Anchor: the primary-source closed form predicts ≈1.2566e+11
     * at these inputs. Assert this independent of the recoded
     * comparison so the test catches any drift from Skolnik. */
    K26AstroJammer *j = k26astro_jammer_new(
        K26ASTRO_JAMMER_MODE_NOISE,
        1.0e3, 25.0, 10.0e9, 1.0e6, 1.0);
    const double js = k26astro_jammer_js_ratio(
        j, 1.0e4, 35.0, 1.0, 1.0e6);
    const double G_j = pow(10.0, 2.5);
    const double G_t = pow(10.0, 3.5);
    const double js_expected = (1.0e3 * G_j * 4.0 * K26A_PI * 1.0e12) /
                                (1.0e4 * G_t * 1.0);
    assert(approx_(js, js_expected, 1.0e-12));
    assert(approx_(js, 1.2566e+11, 1.0e-3));  /* Skolnik §15.3 anchor */

    /* Burn-through range: J/S(R_BT) = snr_threshold = 1.
     *   R_BT = √(snr_threshold · P_t · G_t · σ / (P_j · G_j · 4π))
     * Closed form at these inputs:
     *   R_BT² = (1 × 1e4 × 3162.28 × 1) / (1e3 × 316.23 × 12.566)
     *         = 3.16228e7 / 3.974e6
     *         ≈ 7.957
     *   R_BT  ≈ 2.821 m.
     *
     * Anchor: Skolnik §15.3 burn-through closed form at these inputs. */
    const double R_BT = k26astro_jammer_burn_through_range(
        j, 1.0e4, 35.0, 1.0);
    const double R_BT_expected = sqrt(
        (1.0 * 1.0e4 * G_t * 1.0) /
        (1.0e3 * G_j * 4.0 * K26A_PI));
    assert(approx_(R_BT, R_BT_expected, 1.0e-12));
    assert(approx_(R_BT, 2.821, 1.0e-3));  /* Skolnik §15.3 anchor */

    /* Self-signature equals transmitted power. */
    assert(k26astro_jammer_self_signature_W(j) == 1.0e3);

    /* Counter-detection consistency: J/S ratio at R_BT equals
     * exactly snr_threshold. */
    const double js_at_BT = k26astro_jammer_js_ratio(
        j, 1.0e4, 35.0, 1.0, R_BT);
    assert(approx_(js_at_BT, 1.0, 1.0e-9));

    k26astro_jammer_destroy(j);

    /* ---- Chaff χ²(2) statistics ----------------------------- *
     *
     * Mean RCS = N σ_dipole. CDF F(x) = 1 - exp(-x/2).
     * Quantile at p = 0.5 → -2 ln(0.5) ≈ 1.386. */
    const double mean = k26astro_chaff_mean_rcs(1000, 1e-3);
    assert(mean == 1.0);

    assert(approx_(k26astro_chaff_chi2_2_cdf(2.0),
                   1.0 - exp(-1.0), 1.0e-12));
    assert(approx_(k26astro_chaff_chi2_2_quantile(0.5),
                   2.0 * log(2.0), 1.0e-12));
    /* CDF and quantile are inverses. */
    const double x_test = 3.0;
    const double p_test = k26astro_chaff_chi2_2_cdf(x_test);
    const double x_round = k26astro_chaff_chi2_2_quantile(p_test);
    assert(approx_(x_round, x_test, 1.0e-12));

    /* Sample mean over 10 000 draws should be close to the mean RCS. */
    K26CRng rng;
    k26c_rng_init(&rng, 42);
    double sum = 0.0;
    const int N = 10000;
    for (int i = 0; i < N; i++) {
        sum += k26astro_chaff_sample_rcs(1000, 1e-3, &rng);
    }
    const double sample_mean = sum / N;
    /* Sample-mean of N draws from Exp distribution has σ ≈ mean / √N.
     * Allow 5% relative tolerance at N = 10 000. */
    assert(approx_(sample_mean, mean, 0.05));

    /* NULL RNG returns the deterministic mean. */
    assert(k26astro_chaff_sample_rcs(1000, 1e-3, NULL) == mean);

    /* Bit-identical replay: same seed → same sample. */
    K26CRng rng_a, rng_b;
    k26c_rng_init(&rng_a, 12345);
    k26c_rng_init(&rng_b, 12345);
    for (int i = 0; i < 100; i++) {
        const double a = k26astro_chaff_sample_rcs(500, 2e-3, &rng_a);
        const double b = k26astro_chaff_sample_rcs(500, 2e-3, &rng_b);
        assert(a == b);
    }

    printf("test_softkill_physics: OK\n");
    return 0;
}
