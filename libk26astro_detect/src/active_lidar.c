/* active_lidar.c — photon-counting lidar detection.
 *
 * Returned photon count per pulse from a Lambertian-scattering
 * target at range R, with transmitter gain G_t (Measures 1992,
 * Chapter 7):
 *
 *     N_signal = (E_pulse · G_t · η · ρ · A_tgt · cos θ · A_rx · T_atm)
 *                / (4π² · R⁴ · E_photon)
 *
 * Derivation:
 *   - Outbound: irradiance at target = E_pulse · G_t / (4π R²).
 *   - Lambertian reflection of (A_tgt · cos θ · ρ) of that irradiance,
 *     with intensity per steradian = (P_re-emitted / π).
 *   - Return-leg: irradiance at receiver = (P/π) / R²; capture A_rx.
 *   - Multiply by detector efficiency and divide by E_photon.
 *
 * Shot noise is Poisson: SNR = √N for mean-rate N. The optional
 * RNG returns Poisson-sampled photon counts in REFERENCED-mode
 * runs; passing NULL returns the deterministic mean. */

#include "k26astro_detect/detect_active_lidar.h"

#include "k26astro_core/consts.h"
#include "k26compute.h"

#include <math.h>
#include <stddef.h>

K26AstroDetectLidarEvent k26astro_detect_lidar_active(
    double pulse_energy_J,
    double wavelength_nm,
    double aperture_rx_m,
    double g_tx_db,
    double target_albedo,
    double target_area_m2,
    double cos_view_angle,
    double range_m,
    double atmospheric_tx,
    double detector_efficiency,
    double snr_threshold,
    K26CRng *rng)
{
    /* Legacy entry: diffraction-limited (M²=1.0). Equivalent to
     * _with_beam_quality(... m_squared=1.0). */
    return k26astro_detect_lidar_active_with_beam_quality(
        pulse_energy_J, wavelength_nm, aperture_rx_m, g_tx_db,
        /* m_squared = */ 1.0,
        target_albedo, target_area_m2, cos_view_angle, range_m,
        atmospheric_tx, detector_efficiency, snr_threshold, rng);
}

K26AstroDetectLidarEvent k26astro_detect_lidar_active_with_beam_quality(
    double pulse_energy_J,
    double wavelength_nm,
    double aperture_rx_m,
    double g_tx_db,
    double m_squared,
    double target_albedo,
    double target_area_m2,
    double cos_view_angle,
    double range_m,
    double atmospheric_tx,
    double detector_efficiency,
    double snr_threshold,
    K26CRng *rng)
{
    K26AstroDetectLidarEvent ev;
    ev.range_m      = range_m;
    ev.mean_photons = 0.0;
    ev.snr          = 0.0;
    ev.detected     = 0;

    if (pulse_energy_J <= 0.0 || wavelength_nm <= 0.0) return ev;
    if (aperture_rx_m <= 0.0 || target_area_m2 <= 0.0) return ev;
    if (range_m <= 0.0) return ev;
    if (target_albedo <= 0.0) return ev;
    if (target_albedo > 1.0) target_albedo = 1.0;
    if (atmospheric_tx <= 0.0) return ev;
    if (atmospheric_tx > 1.0) atmospheric_tx = 1.0;
    if (detector_efficiency <= 0.0) return ev;
    if (detector_efficiency > 1.0) detector_efficiency = 1.0;
    if (cos_view_angle <= 0.0) return ev;
    if (cos_view_angle > 1.0) cos_view_angle = 1.0;
    /* M² < 1 violates the diffraction limit; clamp up. */
    if (m_squared < 1.0) m_squared = 1.0;

    const double lambda_m = wavelength_nm * 1.0e-9;
    const double E_photon = K26A_H_PLANCK * K26A_C / lambda_m;
    const double A_rx = K26A_PI * 0.25 * aperture_rx_m * aperture_rx_m;
    const double R4 = range_m * range_m * range_m * range_m;
    const double G_t_dl = pow(10.0, g_tx_db / 10.0);
    /* M² degradation: transmitter peak intensity per steradian
     * scales as 1/M² (beam divergence half-angle scales as √M²,
     * solid-angle as M², on-axis intensity as 1/M²). */
    const double G_t_effective = G_t_dl / m_squared;

    /* Round-trip Lambertian return:
     *   (E_pulse · G_t · η · ρ · A_tgt · cos θ · A_rx · T_atm)
     *   ÷ (4π² · R⁴ · E_photon) */
    double N = (pulse_energy_J / E_photon) *
               G_t_effective *
               detector_efficiency *
               target_albedo *
               target_area_m2 *
               cos_view_angle *
               A_rx *
               atmospheric_tx /
               (4.0 * K26A_PI * K26A_PI * R4);
    if (N < 0.0) N = 0.0;

    if (rng) {
        /* Poisson approximated by Gaussian with σ = √N for N >= 1. */
        const double sigma = (N > 1.0) ? sqrt(N) : 1.0;
        N += k26c_rng_normal(rng, 0.0, sigma);
        if (N < 0.0) N = 0.0;
    }

    ev.mean_photons = N;
    ev.snr      = (N > 0.0) ? sqrt(N) : 0.0;
    ev.detected = (ev.snr >= snr_threshold) ? 1 : 0;
    return ev;
}
