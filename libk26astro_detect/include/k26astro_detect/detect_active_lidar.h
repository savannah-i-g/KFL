/* detect_active_lidar.h — photon-counting lidar detection.
 *
 * The round-trip lidar equation in the photon-counting regime:
 *
 *     N_signal = (E_pulse · G_t · η · ρ · A_rx · A_tgt · cos θ · T_atm)
 *                ÷ (4π² · R⁴ · E_photon)
 *
 * where:
 *     E_pulse    transmitted pulse energy (J)
 *     G_t        transmitter antenna gain (linear; from g_tx_db)
 *     λ          wavelength (m)
 *     ρ          target albedo (Lambertian; in [0, 1])
 *     A_rx       receiver aperture area (m²)
 *     A_tgt      target projected area (m²)
 *     cos θ      cos of target normal vs look vector (Lambertian)
 *     T_atm      two-way atmospheric transmissivity (1.0 vacuum)
 *     R          range (m)
 *     η          detector quantum efficiency
 *     E_photon   per-photon energy (h c / λ)
 *
 * Derivation:
 *   - Outbound: irradiance at target = E_pulse · G_t / (4π R²).
 *   - Lambertian reflection: re-emits A_tgt · cos θ · ρ of incident
 *     power, distributed as Lambertian (intensity = P/π per sr).
 *   - Return-leg: irradiance at receiver = (re-emitted/π) / R²,
 *     times A_rx captures the photon count.
 *
 * For a diffraction-limited transmit aperture, G_t = (πD/λ)² ≈
 * (π · 1 m / 1.064 μm)² ≈ 8.72e+12 ≈ 129.4 dB at 1064 nm. Real
 * systems take a few-dB hit from beam quality and pointing.
 *
 * Photon shot noise is Poisson; detection threshold compares the
 * expected photon count to the SNR threshold under Poisson
 * statistics: SNR = √N for sufficiently large N.
 *
 * The two-way Beer-Lambert atmospheric loss factor is exposed as
 * an opt-in parameter; vacuum-bubble engagements set it to 1.0.
 *
 * Reference: Measures, R.M. 1992. Laser Remote Sensing:
 * Fundamentals and Applications. Krieger. */
#ifndef K26ASTRO_DETECT_ACTIVE_LIDAR_H
#define K26ASTRO_DETECT_ACTIVE_LIDAR_H

#include "k26compute.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double range_m;
    double mean_photons;
    double snr;
    int    detected;
} K26AstroDetectLidarEvent;

/* Compute a lidar return + detection decision.
 *
 * Inputs:
 *   pulse_energy_J     transmitted pulse energy
 *   wavelength_nm      operating wavelength (nm)
 *   aperture_rx_m      receiver aperture diameter (full diameter,
 *                      not area)
 *   g_tx_db            transmitter antenna gain in dB. For a
 *                      diffraction-limited uniform aperture
 *                      G_t = (πD/λ)²; typical near-130 dB for a
 *                      1-metre Nd:YAG aperture at 1064 nm.
 *   target_albedo      Lambertian albedo, [0, 1]
 *   target_area_m2     target projected area facing the lidar
 *   cos_view_angle     cos of target normal vs return path
 *   range_m            range to target
 *   atmospheric_tx     two-way atmospheric transmissivity in [0,1];
 *                      1.0 for vacuum
 *   detector_efficiency  detector quantum efficiency, [0, 1]
 *   snr_threshold      detection-trigger SNR threshold
 *   rng                optional RNG for shot-noise jitter (NULL for
 *                      mean-only deterministic output). */
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
    K26CRng *rng);

/* Compute a lidar return with explicit beam-quality degradation.
 *
 * Identical to k26astro_detect_lidar_active with one extra
 * parameter: m_squared, the laser beam-quality factor (M²).
 *
 * A diffraction-limited single-transverse-mode Gaussian beam has
 * M²=1.0. Real laser transmitters degrade by a factor that
 * captures multi-mode content, aberrations, and pointing jitter
 * not separately modelled:
 *
 *   M² ≈ 1.0      ideal diffraction-limited oscillator
 *   M² = 1.2-1.5  good-quality solid-state TEM₀₀ laser
 *   M² = 1.5-3.0  high-power slab/disk lasers; fibre lasers
 *   M² > 3.0      multi-mode, poor beam quality
 *
 * The transmitter gain is scaled as G_t_effective = G_t / M²
 * (power per unit solid angle scales inversely with the
 * beam-divergence squared, which scales with M²). The returned
 * photon count is therefore reduced by the same factor.
 *
 * For M²=1.0 the result is bit-identical to the diffraction-
 * limited legacy entry. The legacy
 * k26astro_detect_lidar_active is now a thin wrapper that
 * passes m_squared=1.0. */
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
    K26CRng *rng);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_ACTIVE_LIDAR_H */
