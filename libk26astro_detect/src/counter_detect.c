/* counter_detect.c — active-emission counter-detection range.
 *
 * Solves the SNR-equals-threshold condition of detect_ir_passive
 * for R, given the emitter's self-emitted thermal power. The
 * emitter's transmitted power P_tx_W is multiplied by the
 * dissipation factor to obtain the effective thermal load on the
 * platform's radiator; the radiator's effective temperature T_K
 * is the caller's input (set by the platform thermal model).
 *
 * For the photon-shot-noise dominated regime where N_signal >>
 * N_background:
 *     SNR^2 = N_signal = (P_eff · η · A_obs · τ) /
 *                        (4π R^2 · E_photon)
 *
 * Setting SNR = snr_threshold and solving for R:
 *     R = √( P_eff · η · A_obs · τ / (4π · snr_threshold^2 · E_photon) )
 *
 * The dissipation factor is the detect_consts.h compile-time
 * default; v0.1 exposes no per-call override.
 */

#include "k26astro_detect/counter_detect.h"
#include "k26astro_detect/detect_consts.h"
#include "k26astro_detect/signature.h"

#include "k26astro_core/consts.h"

#include <math.h>

double k26astro_counter_detect_ir_range(
    double emitter_p_tx_W,
    double emitter_radiator_T_K,
    double observer_aperture_m,
    double observer_integration_s,
    double observer_passband_lo_um,
    double observer_passband_hi_um,
    double observer_throughput,
    double snr_threshold)
{
    if (emitter_p_tx_W <= 0.0 || observer_aperture_m <= 0.0)
        return 0.0;
    if (emitter_radiator_T_K <= 0.0) return 0.0;
    if (observer_integration_s <= 0.0) return 0.0;
    if (observer_passband_hi_um <= observer_passband_lo_um) return 0.0;
    if (observer_throughput <= 0.0) return 0.0;
    if (observer_throughput > 1.0) observer_throughput = 1.0;
    if (snr_threshold <= 0.0) return 0.0;

    /* Fraction of the radiator's bolometric emission that falls in
     * the observer's pass-band. With ε=1, area=1 the Planck helper
     * returns π·∫B_λ dλ (W/m²); the bolometric per-area Stefan-
     * Boltzmann emission is σT⁴, so f_in = π·radiance_inband / σT⁴
     * is a dimensionless fraction in [0, 1]. */
    const double lo_m = 1.0e-6 * observer_passband_lo_um;
    const double hi_m = 1.0e-6 * observer_passband_hi_um;
    const double inband =
        k26astro_signature_ir_planck_inband(
            1.0, 1.0, emitter_radiator_T_K, lo_m, hi_m, 64);
    const double bolo = K26ASTRO_DETECT_SIGMA_SB *
        emitter_radiator_T_K * emitter_radiator_T_K *
        emitter_radiator_T_K * emitter_radiator_T_K;
    double f_in = (bolo > 0.0) ? (inband / bolo) : 0.0;
    if (f_in < 0.0) f_in = 0.0;
    if (f_in > 1.0) f_in = 1.0;

    const double P_eff = emitter_p_tx_W *
        K26ASTRO_DETECT_COUNTER_DETECT_DISSIPATION_FACTOR *
        f_in;

    const double A_obs = K26A_PI * 0.25 *
        observer_aperture_m * observer_aperture_m;

    /* Pass-band central wavelength for the per-photon energy. */
    const double lambda_mean_m = 0.5e-6 *
        (observer_passband_lo_um + observer_passband_hi_um);
    const double E_photon =
        K26ASTRO_DETECT_PHOTON_ENERGY_AT(lambda_mean_m);

    const double numer = P_eff * observer_throughput *
                         A_obs * observer_integration_s;
    const double denom = 4.0 * K26A_PI *
                         snr_threshold * snr_threshold * E_photon;
    if (denom <= 0.0) return 0.0;
    return sqrt(numer / denom);
}
