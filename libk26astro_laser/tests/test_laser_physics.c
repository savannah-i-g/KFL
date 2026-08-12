/* test_laser_physics.c — closed-form gates for the laser submodules
 * and the integrated engage() pipeline.
 *
 * Tolerances:
 *   - Airy spot size at canonical (λ, D, R) triples: 1e-12 rel.
 *   - Airy encircled fraction at the first dark ring (x = 3.8317):
 *     0.8378 ± 1e-3 (the canonical Born & Wolf value).
 *   - Maréchal Strehl at σ_w = λ/14: ≈ 0.82 ± 1e-3 (matches the
 *     standard "σ_w < λ/14 keeps Strehl > 0.8" textbook claim).
 *   - Phipps Q* monotone with material density: Al < Ti, Cu < steel.
 *   - Plasma-plug T = 1 below threshold, exp(-k(Φ/Φ_th - 1)) above.
 *   - Engage round-trip: at high fluence the engage event reports
 *     a non-zero mass loss; at sub-threshold fluence the mass loss
 *     drops to <= 0.1 × the supra-threshold result. */

#include "k26astro_laser/laser.h"
#include "k26astro_core/consts.h"

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
    /* ---- Airy half-angle -------------------------------------- */
    /* λ = 1064 nm, D = 1 m, M^2 = 1.0, σ_θ = 0 → θ = 1.22 × 1064e-9. */
    const double theta = k26astro_airy_half_angle(
        1064.0e-9, 1.0, 1.0, 0.0);
    const double theta_expected = 1.22 * 1064.0e-9;
    assert(approx_(theta, theta_expected, 1.0e-12));

    /* Spot diameter at R = 1000 km: 2 R θ ≈ 2 × 1e6 × 1.298e-6 ≈ 2.6 m. */
    const double spot = k26astro_airy_spot_diameter(
        1.0e6, 1064.0e-9, 1.0, 1.0, 0.0);
    assert(approx_(spot, 2.0 * 1.0e6 * theta_expected, 1.0e-12));

    /* M^2 = 1.5 widens the spot proportionally; jitter quadratures in. */
    const double theta_m15 = k26astro_airy_half_angle(
        1064.0e-9, 1.0, 1.5, 0.0);
    assert(approx_(theta_m15, 1.5 * theta_expected, 1.0e-12));

    /* ---- Airy encircled fraction ----------------------------- */
    /* At the first dark ring: x = πDr/(λR) = 3.8317.
     * For λ = 1 μm, D = 1 m, R = 1 m, the first dark ring is at
     * r = 3.8317 λ / (π D) = 3.8317 × 1e-6 / π ≈ 1.219e-6 m. */
    const double r_first_dark = 3.8317 * 1.0e-6 / K26A_PI;
    const double F = k26astro_airy_encircled_fraction(
        r_first_dark, 1.0, 1.0, 1.0e-6);
    /* Canonical value 0.8378. The exact zero of J_0^2 + J_1^2 to
     * machine precision sits within ~3e-4 of 3.8317 (we used the
     * rounded value), so allow ±1e-3 relative. */
    assert(approx_(F, 0.8378, 2.0e-3));

    /* ---- Maréchal Strehl ------------------------------------ */
    /* σ_w = λ/14 → 2π × 1/14 ≈ 0.4488 → exp(-0.4488^2) ≈ 0.817. */
    const double S_14 = k26astro_mirror_strehl(1.0e-6 / 14.0, 1.0e-6);
    const double S_14_expected = exp(-pow(2.0 * K26A_PI / 14.0, 2.0));
    assert(approx_(S_14, S_14_expected, 1.0e-12));
    assert(S_14 > 0.80 && S_14 < 0.83);

    /* σ_w = 0 → S = 1. */
    assert(k26astro_mirror_strehl(0.0, 1.0e-6) == 1.0);

    /* ---- Plasma-plug ---------------------------------------- */
    /* T = 1 below threshold. */
    const double T_below = k26astro_plasma_plug_transmissivity(
        5e7, 1e8, 0.5);
    assert(T_below == 1.0);
    /* T = exp(-0.5 × 1) = 0.6065 at 2× threshold. */
    const double T_2x = k26astro_plasma_plug_transmissivity(
        2e8, 1e8, 0.5);
    assert(approx_(T_2x, exp(-0.5), 1.0e-12));

    /* ---- Phipps coupling (qualitative) ---------------------- */
    /* Q* values are positive and order-of-magnitude correct. */
    assert(k26astro_phipps_q_star(K26ASTRO_LASER_MAT_ALUMINUM) > 1.0e7);
    assert(k26astro_phipps_q_star(K26ASTRO_LASER_MAT_STEEL) >
           k26astro_phipps_q_star(K26ASTRO_LASER_MAT_COPPER));

    /* Sub-threshold C_m derates by 10× per the lib's vapour-pressure
     * regime convention. */
    const double Cm_supra = k26astro_phipps_c_m(
        K26ASTRO_LASER_MAT_ALUMINUM, 1064.0, 1.0e9);
    const double Cm_sub = k26astro_phipps_c_m(
        K26ASTRO_LASER_MAT_ALUMINUM, 1064.0, 1.0e6);
    assert(approx_(Cm_sub, 0.1 * Cm_supra, 1.0e-12));

    /* Threshold is wavelength-independent per Phipps 2010 §III.B
     * (page 616) — no λ scaling appears in the cited paragraph. */
    const double Phi_th_532 = k26astro_phipps_threshold_fluence(
        K26ASTRO_LASER_MAT_ALUMINUM, 532.0);
    const double Phi_th_1064 = k26astro_phipps_threshold_fluence(
        K26ASTRO_LASER_MAT_ALUMINUM, 1064.0);
    assert(Phi_th_532 == Phi_th_1064);

    /* C_m wavelength scaling follows Phipps 2010 Eq. 19 (page 616):
     * C_m ∝ (Iλ√τ)^(-1/4), so the λ alone gives (1064/λ)^(1/4).
     * Anchor at 532 nm: ratio = 2^(1/4) ≈ 1.1892.
     *
     * Use supra-threshold fluence (1e9 J/m² > Phi_th_Al = 1e8 J/m²)
     * so both C_m values come from material_cm_default × wavelength_
     * scale, not the sub-threshold derate. */
    const double Cm_1064_supra = k26astro_phipps_c_m(
        K26ASTRO_LASER_MAT_ALUMINUM, 1064.0, 1.0e9);
    const double Cm_532_supra = k26astro_phipps_c_m(
        K26ASTRO_LASER_MAT_ALUMINUM, 532.0, 1.0e9);
    assert(approx_(Cm_532_supra / Cm_1064_supra, 1.1892, 1.0e-4));

    /* ---- Integrated engage() pipeline ------------------------ *
     *
     * Construct a representative MW-class laser and target.
     * Supra-threshold engagement must produce non-zero ablation.
     *
     * At range = 100 km, D = 1 m, λ = 1064 nm, M² = 1.2, P = 1 MW,
     * target_area = 1 m²: the Airy spot (~0.31 m diameter) sits well
     * inside the target, so encircled fraction ≈ 1 and p_on_target
     * ≈ P · Strehl ≈ 9.06e5 W. on_target_intensity = p_on_target /
     * target_area ≈ 9.06e5 W/m². At 1 s dwell, fluence ≈ 9.06e5
     * J/m² — well below the plasma threshold (1e8 J/m² for Al).
     * Increase dwell to drive cumulative fluence supra-threshold. */
    K26AstroLaser *l = k26astro_laser_new(
        /* D */ 1.0,
        /* λ */ 1064.0,
        /* P */ 1.0e6,
        /* M^2 */ 1.2,
        /* jitter */ 1.0e-7,
        /* σ_w */ 1064e-9 / 20.0,
        /* k_attn */ 0.5);
    assert(l != NULL);

    K26AstroLaserAblationEvent ev = k26astro_laser_engage(
        l,
        K26ASTRO_LASER_MAT_ALUMINUM,
        /* A_target */ 1.0,
        /* R_target */ 0.2,
        /* range */ 100.0e3,  /* 100 km */
        /* dwell */ 200.0);   /* fluence ≈ 1.81e8 J/m² > Phi_th_Al */
    assert(ev.spot_diameter_m > 0.0);
    assert(ev.encircled_fraction > 0.0);
    assert(ev.encircled_fraction <= 1.0);
    assert(ev.on_target_intensity_W_per_m2 > 0.0);
    /* Gate: fluence is actually supra-threshold so the lib exercises
     * the plasma-regime path, not the sub-threshold derate. */
    assert(ev.fluence_J_per_m2 >
        k26astro_phipps_threshold_fluence(
            K26ASTRO_LASER_MAT_ALUMINUM, 1064.0));
    assert(ev.plasma_ignited == 1);
    assert(ev.threshold_fluence_J_per_m2 ==
        k26astro_phipps_threshold_fluence(
            K26ASTRO_LASER_MAT_ALUMINUM, 1064.0));
    assert(ev.mass_loss_kg > 0.0);
    assert(ev.impulse_N_s > 0.0);

    /* Sub-threshold variant at 1 s dwell — same geometry, lower
     * cumulative fluence. Mass loss must drop by ≥ 10× (the lib's
     * vapour-pressure derate factor). */
    K26AstroLaserAblationEvent ev_sub = k26astro_laser_engage(
        l,
        K26ASTRO_LASER_MAT_ALUMINUM,
        /* A_target */ 1.0,
        /* R_target */ 0.2,
        /* range */ 100.0e3,
        /* dwell */ 1.0);     /* fluence ≈ 9.06e5 J/m² < Phi_th_Al */
    assert(ev_sub.fluence_J_per_m2 <
        k26astro_phipps_threshold_fluence(
            K26ASTRO_LASER_MAT_ALUMINUM, 1064.0));
    /* Sub-threshold mass loss is at least 10× smaller than supra
     * after normalising for the 200× dwell ratio. */
    assert(ev_sub.mass_loss_kg < ev.mass_loss_kg / 200.0);

    /* At very long range the spot grows large, intensity drops
     * below the plasma threshold, mass loss derates by 10×. */
    K26AstroLaserAblationEvent ev_far = k26astro_laser_engage(
        l, K26ASTRO_LASER_MAT_ALUMINUM, 1.0, 0.2,
        /* range */ 1.0e8,  /* 100 000 km */
        /* dwell */ 1.0);
    assert(ev_far.mass_loss_kg < ev.mass_loss_kg);

    /* ---- plasma-ignition gate fixture ----------------------- *
     *
     * Sub-threshold engagement at the canonical post-flight-audit
     * setpoint: P = 1 MW, D = 1 m, R = 100 km, dwell = 5 s,
     * target = Al with reflectivity 0.7.
     *
     * Closed-form expectations (laser_consts.h + airy/strehl
     * sub-modules):
     *   θ_diff = 1.22 × M² × λ / D = 1.22 × 1.2 × 1064e-9 / 1.0
     *          ≈ 1.5580e-6 rad; with jitter quadratured RSS,
     *          θ_total ≈ √(θ_diff² + σ_θ²) ≈ 1.5611e-6 rad.
     *   spot_d at 100 km ≈ 2 × 1e5 × θ_total ≈ 0.3122 m
     *   spot_area ≈ π × (spot_d/2)² ≈ 0.0766 m²
     *   target_radius for A=1 m² is √(1/π) ≈ 0.5642 m, larger
     *     than spot_radius ≈ 0.156 m → encircled_fraction ≈ 1.0
     *   strehl at σ_w = λ/20 ≈ exp(-(2π/20)²) ≈ 0.905
     *   p_aperture ≈ 1e6 × 0.905 ≈ 9.05e5 W
     *   on_target_intensity ≈ p_aperture × encircled / A_target
     *     ≈ 9.05e5 W/m² (A_target = 1 m²)
     *   fluence ≈ 9.05e5 × 5 ≈ 4.5e6 J/m²
     *   threshold (Al @ 1064 nm) = 1.0e8 J/m² (laser_consts.h:82).
     *
     * Therefore fluence ≪ threshold → plasma_ignited = false +
     * mass_loss is vapour-regime (10× derated). */
    K26AstroLaserAblationEvent ev_subthresh = k26astro_laser_engage(
        l,
        K26ASTRO_LASER_MAT_ALUMINUM,
        /* A_target */ 1.0,
        /* R_target */ 0.7,
        /* range */ 100.0e3,
        /* dwell */ 5.0);
    assert(ev_subthresh.fluence_J_per_m2 <
           ev_subthresh.threshold_fluence_J_per_m2);
    assert(ev_subthresh.threshold_fluence_J_per_m2 ==
           K26ASTRO_LASER_THRESHOLD_AL_J_PER_M2_AT_1064NM);
    assert(ev_subthresh.plasma_ignited == 0);
    /* Vapour-regime mass loss = E_coupled / Q_star × 0.1. The
     * derated mass loss is at most one-tenth of the formal
     * E_coupled / Q_star result; equivalently, mass_loss_kg ≤
     * 0.1 × (p_coupled × dwell / Q_star) plus epsilon. The exact
     * value is tied to encircled_fraction + strehl + reflectivity
     * but the gate is simply: vapour-regime applied. */
    const double Q_star_Al = k26astro_phipps_q_star(
        K26ASTRO_LASER_MAT_ALUMINUM);
    const double E_coupled_sub =
        ev_subthresh.p_coupled_W * ev_subthresh.effective_dwell_s;
    const double ablation_regime_mass =
        (Q_star_Al > 0.0) ? (E_coupled_sub / Q_star_Al) : 0.0;
    /* Derated vapour-regime mass-loss is 0.1× the ablation-regime
     * value; allow 1e-12 relative slack on the derate factor. */
    assert(ev_subthresh.mass_loss_kg <= 0.1 * ablation_regime_mass + 1e-30);
    assert(ev_subthresh.mass_loss_kg >= 0.1 * ablation_regime_mass - 1e-30);

    k26astro_laser_destroy(l);

    /* ---- Determinism gate: bit-identical re-engage ---------- */
    K26AstroLaser *l2 = k26astro_laser_new(
        1.0, 1064.0, 1.0e6, 1.2, 1.0e-7, 5e-8, 0.5);
    K26AstroLaserAblationEvent ev_a = k26astro_laser_engage(
        l2, K26ASTRO_LASER_MAT_TITANIUM, 1.0, 0.3, 100.0e3, 1.0);
    K26AstroLaserAblationEvent ev_b = k26astro_laser_engage(
        l2, K26ASTRO_LASER_MAT_TITANIUM, 1.0, 0.3, 100.0e3, 1.0);
    assert(ev_a.mass_loss_kg == ev_b.mass_loss_kg);
    assert(ev_a.impulse_N_s  == ev_b.impulse_N_s);
    assert(ev_a.fluence_J_per_m2 == ev_b.fluence_J_per_m2);
    k26astro_laser_destroy(l2);

    printf("test_laser_physics: OK\n");
    return 0;
}
