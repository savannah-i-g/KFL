/* test_detect_physics.c — closed-form physics gates.
 *
 * Tolerance choices:
 *   - Stefan-Boltzmann: ulp floor (1e-12 rel) — closed-form.
 *   - Planck inband at 5772 K, full pass-band: must agree with
 *     the bolometric Stefan-Boltzmann result to 1e-3 rel (limited
 *     by Simpson's rule + finite pass-band). The pass band covers
 *     the bulk of the solar spectrum.
 *   - Radar equation: ulp floor (1e-12 rel) — closed-form.
 *   - Lidar return: ulp floor — closed-form when atmospheric_tx=1.
 *   - Counter-detection: round-trip consistency — if observer
 *     parameters detect a passive target of given power at range R,
 *     and we ask "at what range would observer detect an emitter
 *     of the same effective power", the answer must equal R.
 *
 * Determinism: all calls pass NULL for the RNG so results are
 * the mean values, comparable across runs at ulp floor. */

#include "k26astro_detect/detect.h"
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
    /* ---- Stefan-Boltzmann bolometric --------------------------- */
    /* Sphere at 300 K, A = 1 m^2, ε = 1: P = σ T^4 ≈ 459.3 W. */
    const double T1 = 300.0;
    const double P1 = k26astro_signature_ir_bolometric_power(1.0, 1.0, T1);
    const double P1_expected = K26ASTRO_DETECT_SIGMA_SB * T1 * T1 * T1 * T1;
    assert(approx_(P1, P1_expected, 1.0e-12));
    /* Numerical value sanity check. */
    assert(P1 > 459.0 && P1 < 460.0);

    /* Black hole (T=0) emits zero. */
    assert(k26astro_signature_ir_bolometric_power(1.0, 1.0, 0.0) == 0.0);
    /* Greybody ε = 0.5 halves the power. */
    const double P_grey = k26astro_signature_ir_bolometric_power(1.0, 0.5, T1);
    assert(approx_(P_grey, 0.5 * P1, 1.0e-12));

    /* ---- Lambertian intensity --------------------------------- */
    /* Edge-on view: cos θ = 0 → I = 0. */
    assert(k26astro_signature_ir_lambertian(1.0, 1.0, T1, 0.0) == 0.0);
    /* Normal view: I = P / π. */
    const double I_normal = k26astro_signature_ir_lambertian(1.0, 1.0, T1, 1.0);
    assert(approx_(I_normal, P1 / K26A_PI, 1.0e-12));

    /* ---- Planck integration ----------------------------------- */
    /* Integrate over a wide pass-band that captures essentially
     * the entire blackbody spectrum (1 nm to 1 mm). Result must
     * approach the Stefan-Boltzmann bolometric value within
     * Simpson's quadrature error at 256 subdivisions. */
    const double P_full = k26astro_signature_ir_planck_inband(
        1.0, 1.0, T1,
        1.0e-9,    /* λ_lo = 1 nm */
        1.0e-3,    /* λ_hi = 1 mm */
        256);
    assert(approx_(P_full, P1, 1.0e-3));

    /* Narrow pass-band around the Wien peak (b/T = 9.66 μm at
     * 300 K) should be a substantial fraction of total power. */
    const double P_3_12 = k26astro_signature_ir_planck_inband(
        1.0, 1.0, T1,
        3.0e-6, 12.0e-6, 256);
    /* At 300 K the 3-12 μm band carries roughly 30-50% of total
     * radiance; the precise fraction depends on the integrator's
     * subdivision count. Assert lower bound 0.2, upper bound 1.0
     * (the band cannot exceed total). */
    assert(P_3_12 > 0.2 * P1 && P_3_12 < P1);

    /* ---- Faceted RCS ------------------------------------------ */
    /* Single 1 m^2 facet, normal along +z, look down -z (observer
     * above target). cos θ = 1, contribution = A^2 = 1.
     *
     * In K26 the look_unit_body vector points from observer
     * toward target. A target with its emitting / scattering
     * normal in +z and an observer above (+z relative to the
     * target) has look_unit_body = -z (observer-to-target points
     * downward in target frame). cos θ in detect_rcs is computed
     * via -look · normal = -(-z) · z = +1. */
    K26V3 normal_z   = { 0.0, 0.0, 1.0 };
    double area_1m2  = 1.0;
    K26V3 look_minus_z = { 0.0, 0.0, -1.0 };
    const double geom = k26astro_signature_rcs_facet_geometry(
        1, &normal_z, &area_1m2, look_minus_z);
    assert(approx_(geom, 1.0, 1.0e-12));

    /* Same facet but look from the side (look = +y): cos = 0 →
     * no contribution. */
    K26V3 look_y = { 0.0, 1.0, 0.0 };
    const double geom_y = k26astro_signature_rcs_facet_geometry(
        1, &normal_z, &area_1m2, look_y);
    assert(geom_y == 0.0);

    /* Monostatic RCS at λ = 0.03 m (X-band, 10 GHz):
     *   σ = (4π / λ^2) · A^2 cos^4 θ = (4π / 9e-4) · 1 ≈ 13963 m^2. */
    const double rcs = k26astro_signature_rcs_monostatic(
        1, &normal_z, &area_1m2, look_minus_z, 0.03);
    const double rcs_expected = (4.0 * K26A_PI / (0.03 * 0.03));
    assert(approx_(rcs, rcs_expected, 1.0e-12));

    /* ---- Radar equation -------------------------------------- */
    /* Reference values: P_t = 1 kW, G_t = G_r = 35 dB (linear =
     * 10^3.5 ≈ 3162), λ = 0.03 m (X-band), σ = 1 m^2, R = 1000 km,
     * L_sys = 4 dB, B = 1 MHz, T_sys = 290 K, F = 3 (≈4.77 dB).
     *
     * Closed-form P_r computed inline; the lib must agree to ulp. */
    const double p_tx = 1.0e3;
    const double g_db = 35.0;
    const double freq = 10.0e9;
    const double rcs_in = 1.0;
    const double R = 1.0e6;
    const double loss_db = 4.0;
    const double bw = 1.0e6;
    const double Tsys = 290.0;
    const double nf = 3.0;

    const double lambda = K26A_C / freq;
    const double G = pow(10.0, g_db / 10.0);
    const double L = pow(10.0, loss_db / 10.0);
    const double four_pi = 4.0 * K26A_PI;
    const double Pr_expected =
        (p_tx * G * G * lambda * lambda * rcs_in) /
        (four_pi * four_pi * four_pi * R * R * R * R * L);
    const double N_expected = K26A_K_BOLTZMANN * Tsys * bw * nf;

    K26AstroDetectRadarEvent re = k26astro_detect_radar_active(
        p_tx, g_db, g_db, freq, rcs_in, R, loss_db, bw, Tsys, nf,
        0.0, NULL);
    assert(approx_(re.power_received_W, Pr_expected, 1.0e-12));
    assert(approx_(re.noise_floor_W,    N_expected,  1.0e-12));
    assert(approx_(re.snr,              Pr_expected / N_expected, 1.0e-12));

    /* ---- Lidar return ---------------------------------------- */
    /* Reference: 1 mJ pulse at 1064 nm, 0.2 m aperture, 0.5 albedo,
     * 1 m^2 target, normal view, 100 km, vacuum, QE = 0.6, transmit
     * gain 130 dB (representative large-laser figure).
     *
     *   N = (E · G_t · η · ρ · A_tgt · cos θ · A_rx · T_atm)
     *       / (4π² · R⁴ · E_photon)
     *
     * E_ph = hc/λ ≈ 1.866e-19 J. */
    const double E_pulse = 1.0e-3;
    const double lambda_nm = 1064.0;
    const double A_rx_dia = 0.2;
    const double g_tx_db_test = 130.0;
    const double rho = 0.5;
    const double A_tgt = 1.0;
    const double Rl = 1.0e5;
    const double eta = 0.6;

    const double lambda_l = lambda_nm * 1.0e-9;
    const double E_ph = K26A_H_PLANCK * K26A_C / lambda_l;
    const double A_rx_area = K26A_PI * 0.25 * A_rx_dia * A_rx_dia;
    const double G_t_lidar = pow(10.0, g_tx_db_test / 10.0);
    const double N_expected_lidar =
        (E_pulse / E_ph) * G_t_lidar * eta * rho * A_tgt * 1.0 * A_rx_area /
        (4.0 * K26A_PI * K26A_PI * Rl * Rl * Rl * Rl);

    K26AstroDetectLidarEvent le = k26astro_detect_lidar_active(
        E_pulse, lambda_nm, A_rx_dia, g_tx_db_test, rho, A_tgt,
        1.0, Rl, 1.0, eta, 5.0, NULL);
    assert(approx_(le.mean_photons, N_expected_lidar, 1.0e-12));

    /* Anchor: diffraction-limited transmit gain for a uniform
     * aperture, G_t = (πD/λ)² (Born & Wolf §8.5 + Skolnik §1.4).
     * For a 1 m aperture at 1064 nm: G_t ≈ 8.72e+12, g_tx ≈ 129.4 dB.
     * Verifies the test's chosen 130 dB is close to a realistic
     * diffraction limit. */
    const double g_tx_diffraction_lim_db =
        20.0 * log10(K26A_PI * 1.0 / (1064.0e-9));
    assert(approx_(g_tx_diffraction_lim_db, 129.4, 0.6));

    /* ---- Counter-detection round-trip ------------------------- *
     *
     * For an observer with given aperture/integration/passband, the
     * counter-detection range against an emitter radiating P_eff
     * thermal watts must be the same range at which the observer
     * would just barely detect a target radiating P_eff watts in
     * the same band. Both derivations solve the same SNR equation,
     * so the round-trip must hold to ulp. */
    const double P_eff_W = 1000.0;
    const double Aper   = 0.3;
    const double Tint   = 0.1;
    const double Plo    = 3.0;
    const double Phi    = 12.0;
    const double Tput   = 0.5;
    const double SNRt   = 5.0;

    /* Counter-detection range. */
    const double R_cd = k26astro_counter_detect_ir_range(
        /* emitter_p_tx */ P_eff_W,
        /* emitter_radiator_T_K */ 300.0,
        /* observer_aperture_m */ Aper,
        /* observer_integration_s */ Tint,
        /* observer_passband_lo_um */ Plo,
        /* observer_passband_hi_um */ Phi,
        /* observer_throughput */ Tput,
        /* snr_threshold */ SNRt);
    assert(R_cd > 0.0);

    /* Confirm by computing the observer's SNR against a target
     * of P_eff at this range — must equal snr_threshold. */
    K26AstroDetectIrEvent ev_cd = k26astro_detect_ir_passive(
        /* emitter_power_W */ P_eff_W,
        /* emitter_T_K */ 300.0,
        /* range_m */ R_cd,
        /* aperture_m */ Aper,
        /* integration_s */ Tint,
        /* passband_lo_um */ Plo,
        /* passband_hi_um */ Phi,
        /* throughput */ Tput,
        /* snr_threshold */ SNRt,
        /* rng */ NULL);
    /* Round-trip tolerance loosens to 5% relative because the
     * passive-IR path adds a small CMB-background floor + a +1
     * dark-current term to the noise that the counter-detection
     * closed form omits. */
    assert(approx_(ev_cd.snr, SNRt, 0.05));

    /* ---- In-band T_K dependence anchor -------------------- *
     *
     * Counter-detection range scales as √(P_eff · f_in(T)), where
     * f_in is the fraction of bolometric Planck radiation in the
     * observer's pass-band. For (3-12) μm, standard radiation-
     * function tables give:
     *   f_in(300 K)  = F(3600 μm·K) - F(900 μm·K)  ≈ 0.418
     *   f_in(1000 K) = F(12000 μm·K) - F(3000 μm·K) ≈ 0.672
     * so the ratio R(1000K) / R(300K) ≈ √(0.672/0.418) ≈ 1.268.
     *
     * Anchor: primary-physics radiation function tables, not the
     * lib's Planck integrator. Tolerance loosens to 5% to absorb
     * Simpson-quadrature rounding at the 64-subdivision default. */
    const double R_cd_300 = k26astro_counter_detect_ir_range(
        P_eff_W, 300.0, Aper, Tint, Plo, Phi, Tput, SNRt);
    const double R_cd_1000 = k26astro_counter_detect_ir_range(
        P_eff_W, 1000.0, Aper, Tint, Plo, Phi, Tput, SNRt);
    assert(R_cd_300 > 0.0 && R_cd_1000 > 0.0);
    const double R_ratio = R_cd_1000 / R_cd_300;
    assert(approx_(R_ratio, 1.268, 0.05));

    /* ---- "No stealth in space" qualitative gate -------------- *
     *
     * A 300 K spacecraft with 100 m^2 radiator area emits roughly
     * 460 W per m^2 × 100 m^2 ≈ 46 kW bolometric. With a 0.3 m
     * observer aperture and 100 ms integration, the counter-
     * detection range should be enormous — many millions of km.
     * This is the popular "there's no stealth in space" result
     * ("There Ain't No Stealth In Space", Atomic Rockets). The gate asserts R > 1e6 m. */
    const double R_stealth = k26astro_counter_detect_ir_range(
        /* emitter_p_tx */ 46000.0,
        /* emitter_radiator_T_K */ 300.0,
        /* observer_aperture_m */ 0.3,
        /* observer_integration_s */ 0.1,
        /* observer_passband_lo_um */ 3.0,
        /* observer_passband_hi_um */ 12.0,
        /* observer_throughput */ 0.5,
        /* snr_threshold */ 5.0);
    assert(R_stealth > 1.0e6);

    printf("test_detect_physics: OK\n");
    return 0;
}
