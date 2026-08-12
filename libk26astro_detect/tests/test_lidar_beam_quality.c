/* test_lidar_beam_quality.c — M² beam-quality factor gate.
 *
 * Asserts the new k26astro_detect_lidar_active_with_beam_quality
 * entry point applies M² beam-quality degradation correctly:
 *
 *   - M² = 1.0 is bit-identical to the legacy diffraction-limited
 *     entry point.
 *   - M² = 1.5 reduces the returned photon count by exactly 1.5×
 *     (transmitter on-axis intensity scales as 1/M²).
 *   - M² < 1.0 is clamped to 1.0 (the diffraction limit).
 *
 * The plan-time motivation: real laser transmitters degrade by
 * M² ∈ [1.2, 2] for beam quality; without this parameter the
 * substrate reports the diffraction-limited upper bound, which
 * makes the writeup's quoted lidar SNR an optimistic anchor
 * rather than a realistic operating point. */
#include "k26astro_detect/detect_active_lidar.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx_rel(double got, double want, double rel_tol)
{
    const double denom = (fabs(want) > 1.0) ? fabs(want) : 1.0;
    return fabs(got - want) / denom <= rel_tol;
}

int main(void)
{
    /* Reference shot: 10 mJ pulse at 1064 nm, 0.3 m receiver
     * aperture, 130 dB transmitter gain (diffraction-limited
     * 1 m primary), Lambertian 0.5 albedo, 2 m² target at
     * cos θ = 1.0 normal incidence, 5×10⁵ m range, vacuum
     * transmissivity, 50% detector efficiency, SNR threshold 5. */
    const double E_pulse_J        = 1.0e-2;
    const double lambda_nm        = 1064.0;
    const double aperture_rx_m    = 0.3;
    const double g_tx_db          = 130.0;
    const double target_albedo    = 0.5;
    const double target_area_m2   = 2.0;
    const double cos_view_angle   = 1.0;
    const double range_m          = 5.0e5;
    const double atmospheric_tx   = 1.0;
    const double detector_eff     = 0.5;
    const double snr_threshold    = 5.0;

    /* ---- Legacy entry: diffraction-limited (M²=1.0) ----------- */
    K26AstroDetectLidarEvent ev_dl =
        k26astro_detect_lidar_active(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(ev_dl.mean_photons > 0.0);
    fprintf(stderr,
            "test_lidar_beam_quality: legacy DL — "
            "N=%.3e SNR=%.3f\n",
            ev_dl.mean_photons, ev_dl.snr);

    /* ---- _with_beam_quality at M²=1.0 should be bit-identical */
    K26AstroDetectLidarEvent ev_m1 =
        k26astro_detect_lidar_active_with_beam_quality(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db,
            /* m_squared = */ 1.0,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(ev_m1.mean_photons == ev_dl.mean_photons);
    assert(ev_m1.snr == ev_dl.snr);

    /* ---- M²=1.5 should reduce photon count by exactly 1.5× ---- */
    K26AstroDetectLidarEvent ev_m15 =
        k26astro_detect_lidar_active_with_beam_quality(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db,
            /* m_squared = */ 1.5,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(approx_rel(ev_m15.mean_photons,
                      ev_dl.mean_photons / 1.5, 1.0e-12));
    /* SNR = √N, so SNR_15 / SNR_dl = √(1/1.5). */
    assert(approx_rel(ev_m15.snr,
                      ev_dl.snr / sqrt(1.5), 1.0e-12));

    fprintf(stderr,
            "test_lidar_beam_quality: M²=1.5 — N=%.3e (= "
            "DL/1.5 = %.3e ✓), SNR=%.3f\n",
            ev_m15.mean_photons, ev_dl.mean_photons / 1.5,
            ev_m15.snr);

    /* ---- M²=2.0 reduces by 2×, M²=3.0 by 3× ------------------- */
    K26AstroDetectLidarEvent ev_m2 =
        k26astro_detect_lidar_active_with_beam_quality(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db, 2.0,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(approx_rel(ev_m2.mean_photons,
                      ev_dl.mean_photons / 2.0, 1.0e-12));

    K26AstroDetectLidarEvent ev_m3 =
        k26astro_detect_lidar_active_with_beam_quality(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db, 3.0,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(approx_rel(ev_m3.mean_photons,
                      ev_dl.mean_photons / 3.0, 1.0e-12));

    /* ---- M² < 1.0 is clamped to 1.0 --------------------------- */
    K26AstroDetectLidarEvent ev_m_below_dl =
        k26astro_detect_lidar_active_with_beam_quality(
            E_pulse_J, lambda_nm, aperture_rx_m, g_tx_db, 0.5,
            target_albedo, target_area_m2, cos_view_angle,
            range_m, atmospheric_tx, detector_eff,
            snr_threshold, NULL);
    assert(ev_m_below_dl.mean_photons == ev_dl.mean_photons);

    fprintf(stderr, "test_lidar_beam_quality: OK\n");
    return 0;
}
