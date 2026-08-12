/* active_rf.c — monostatic radar detection.
 *
 * Implements the standard radar equation (Skolnik 2008 §1.4):
 *
 *     P_r = (P_t · G_t · G_r · λ^2 · σ) / ((4π)^3 · R^4 · L_sys)
 *
 * with G_t, G_r supplied as decibel gains (linearised inside the
 * function) and L_sys supplied as positive dB system loss
 * (linearised similarly). Noise floor uses kTBF.
 *
 * Optional RNG jitter is Gaussian on the received power with
 * standard deviation derived from the noise floor; passing NULL
 * returns the deterministic mean SNR. */

#include "k26astro_detect/detect_active_rf.h"

#include "k26astro_core/consts.h"
#include "k26compute.h"

#include <math.h>
#include <stddef.h>

static double db_to_linear_(double db) { return pow(10.0, db / 10.0); }

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
    K26CRng *rng)
{
    K26AstroDetectRadarEvent ev;
    ev.range_m          = range_m;
    ev.power_received_W = 0.0;
    ev.snr              = 0.0;
    ev.noise_floor_W    = 0.0;
    ev.detected         = 0;

    if (p_tx_W <= 0.0 || freq_Hz <= 0.0 || rcs_m2 <= 0.0 ||
        range_m <= 0.0 || bandwidth_Hz <= 0.0 || T_sys_K <= 0.0)
        return ev;
    if (noise_figure < 1.0) noise_figure = 1.0;

    const double lambda_m = K26A_C / freq_Hz;
    const double G_t = db_to_linear_(g_tx_db);
    const double G_r = db_to_linear_(g_rx_db);
    const double L   = db_to_linear_(loss_sys_db);
    const double four_pi = 4.0 * K26A_PI;
    const double four_pi_cubed = four_pi * four_pi * four_pi;
    const double R4 = range_m * range_m * range_m * range_m;

    double P_r = (p_tx_W * G_t * G_r * lambda_m * lambda_m * rcs_m2) /
                 (four_pi_cubed * R4 * L);
    if (P_r < 0.0) P_r = 0.0;

    /* Receiver noise floor: N = k_B · T_sys · B · F. */
    const double N = K26A_K_BOLTZMANN * T_sys_K * bandwidth_Hz * noise_figure;
    ev.noise_floor_W = N;

    if (rng) {
        const double sigma = N;
        P_r += k26c_rng_normal(rng, 0.0, sigma);
        if (P_r < 0.0) P_r = 0.0;
    }

    ev.power_received_W = P_r;
    ev.snr      = (N > 0.0) ? (P_r / N) : 0.0;
    ev.detected = (ev.snr >= snr_threshold) ? 1 : 0;
    return ev;
}
