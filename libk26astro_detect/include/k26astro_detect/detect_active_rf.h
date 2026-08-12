/* detect_active_rf.h — monostatic radar detection.
 *
 * Standard monostatic radar equation:
 *
 *     P_r = (P_t · G_t · G_r · λ^2 · σ) / ((4π)^3 · R^4 · L_sys)
 *
 * where:
 *     P_r        received power (W)
 *     P_t        transmitted power (W)
 *     G_t, G_r   transmit + receive antenna gains (linear, not dB)
 *     λ          wavelength (m)
 *     σ          target RCS (m^2)
 *     R          observer-target range (m)
 *     L_sys      system loss factor (>= 1 for actual loss)
 *
 * Detection condition: SNR = P_r / N_total >= snr_threshold.
 * Noise floor is kTBF + propagation noise; the lib uses the
 * standard radar receiver-noise model
 *
 *     N_total = k_B · T_sys · B · F
 *
 * with T_sys the receiver system noise temperature (K), B the
 * bandwidth (Hz), and F the noise figure (dimensionless, >= 1).
 *
 * Reference: Skolnik, M.I. 2008. Radar Handbook, 3rd ed., §1.4. */
#ifndef K26ASTRO_DETECT_ACTIVE_RF_H
#define K26ASTRO_DETECT_ACTIVE_RF_H

#include "k26compute.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double range_m;
    double power_received_W;
    double snr;
    double noise_floor_W;
    int    detected;
} K26AstroDetectRadarEvent;

/* Compute a monostatic radar return + detection decision.
 *
 * Inputs:
 *   p_tx_W          transmitted power
 *   g_tx_db         transmit antenna gain (decibels)
 *   g_rx_db         receive antenna gain (decibels)
 *   freq_Hz         operating frequency
 *   rcs_m2          target radar cross-section
 *   range_m         observer-target range
 *   loss_sys_db     system loss (dB, positive number)
 *   bandwidth_Hz    receiver bandwidth
 *   T_sys_K         receiver system noise temperature
 *   noise_figure    receiver noise figure (linear)
 *   snr_threshold   detection-trigger SNR threshold
 *   rng             optional RNG for shot-noise jitter (NULL for
 *                   mean-only deterministic output). */
K26AstroDetectRadarEvent k26astro_detect_radar_active(
    double p_tx_W,
    double g_tx_db,
    double g_rx_db,
    double freq_Hz,
    double rcs_m2,
    double range_m,
    double loss_sys_db,
    double bandwidth_Hz,
    double T_sys_K,
    double noise_figure,
    double snr_threshold,
    K26CRng *rng);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_ACTIVE_RF_H */
