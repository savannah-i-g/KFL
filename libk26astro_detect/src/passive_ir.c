/* passive_ir.c — passive infrared detection geometry.
 *
 * Photon-counting detection model under the assumption that the
 * received in-band irradiance dominates over the CMB background
 * for the engagement-bubble regime of interest (target T >> 2.725 K).
 *
 * Signal photon count:
 *     N_signal = (E_in_band × η × A_observer × τ) / E_photon_mean
 *
 * where E_in_band is the in-band irradiance at the observer
 * (W / m^2 ), η is the optical throughput × QE, A_observer is the
 * receiving aperture area (π D^2 / 4), τ is the integration time,
 * and E_photon_mean is the mean per-photon energy across the pass
 * band (h c / λ_mean for centre wavelength λ_mean).
 *
 * Background photon count is computed identically with the
 * emitter quantity replaced by the CMB-background radiance at
 * 2.725 K. For all practical IR pass-bands and observer apertures
 * the CMB contribution is many orders of magnitude below any
 * thermal target's signature; the background term is retained so
 * the noise model degrades gracefully at very long range or very
 * cold targets.
 *
 * SNR = N_signal / √(N_signal + N_background + 1). Shot-noise
 * dominated; the +1 photon term keeps the denominator finite at
 * zero counts (dark floor).
 *
 * The optional K26CRng* applies Poisson jitter to N_signal in
 * REFERENCED-mode runs; passing NULL returns the deterministic
 * mean. */

#include "k26astro_detect/detect_ir_passive.h"
#include "k26astro_detect/detect_consts.h"
#include "k26astro_detect/signature.h"

#include "k26astro_core/consts.h"
#include "k26compute.h"

#include <math.h>
#include <stddef.h>

K26AstroDetectIrEvent k26astro_detect_ir_passive(
    double emitter_power_W,
    double emitter_T_K,
    double range_m,
    double aperture_m,
    double integration_s,
    double passband_lo_um,
    double passband_hi_um,
    double throughput,
    double snr_threshold,
    K26CRng *rng)
{
    /* Legacy entry: CMB-only background. Equivalent to
     * _with_optics(... t_optics_k=0, optics_emissivity=0). */
    return k26astro_detect_ir_passive_with_optics(
        emitter_power_W, emitter_T_K, range_m, aperture_m,
        integration_s, passband_lo_um, passband_hi_um,
        throughput, snr_threshold,
        /* t_optics_k       */ 0.0,
        /* optics_emissivity*/ 0.0,
        rng);
}

K26AstroDetectIrEvent k26astro_detect_ir_passive_with_optics(
    double emitter_power_W,
    double emitter_T_K,
    double range_m,
    double aperture_m,
    double integration_s,
    double passband_lo_um,
    double passband_hi_um,
    double throughput,
    double snr_threshold,
    double t_optics_k,
    double optics_emissivity,
    K26CRng *rng)
{
    K26AstroDetectIrEvent ev;
    ev.range_m                = range_m;
    ev.flux_received_W_per_m2 = 0.0;
    ev.snr                    = 0.0;
    ev.n_signal_photons       = 0.0;
    ev.n_background_photons   = 0.0;
    ev.detected               = 0;

    if (range_m <= 0.0 || aperture_m <= 0.0 || integration_s <= 0.0)
        return ev;
    if (passband_hi_um <= passband_lo_um) return ev;
    if (throughput <= 0.0) return ev;
    if (throughput > 1.0) throughput = 1.0;

    /* Irradiance at observer: the emitter radiates isotropically
     * into 4π sr (assumed isotropic Lambertian-averaged); the
     * observer sees E_observer = emitter_power_W / (4πR²) of
     * bolometric power per unit area. The in-band fraction is
     * applied below using the emitter temperature; for a directional
     * source the caller should supply the in-band power already
     * integrated over the observer's viewing cone (with f_in = 1). */
    if (emitter_T_K <= 0.0) return ev;
    const double E_observer =
        emitter_power_W / (4.0 * K26A_PI * range_m * range_m);
    ev.flux_received_W_per_m2 = E_observer;

    const double A_observer = K26A_PI * 0.25 * aperture_m * aperture_m;

    /* Mean photon energy at pass-band centre. */
    const double lambda_mean_m =
        0.5e-6 * (passband_lo_um + passband_hi_um);
    if (lambda_mean_m <= 0.0) return ev;
    const double E_photon = K26ASTRO_DETECT_PHOTON_ENERGY_AT(lambda_mean_m);

    /* Pass-band edges in metres; reused for both the emitter in-band
     * fraction and the CMB background term below. */
    const double passband_lo_m = 1.0e-6 * passband_lo_um;
    const double passband_hi_m = 1.0e-6 * passband_hi_um;

    /* Fraction of emitter_power_W that lands in the observer's pass-
     * band, derived from the emitter's Planck spectrum at emitter_T_K. */
    const double inband_emitter =
        k26astro_signature_ir_planck_inband(
            1.0, 1.0, emitter_T_K,
            passband_lo_m, passband_hi_m, 64);
    const double bolo_emitter = K26ASTRO_DETECT_SIGMA_SB *
        emitter_T_K * emitter_T_K * emitter_T_K * emitter_T_K;
    double f_in = (bolo_emitter > 0.0)
        ? (inband_emitter / bolo_emitter) : 0.0;
    if (f_in < 0.0) f_in = 0.0;
    if (f_in > 1.0) f_in = 1.0;

    /* Expected signal photon count, weighted by the in-band fraction
     * of the emitter's bolometric power. */
    double n_sig = (E_observer * f_in * throughput *
                    A_observer * integration_s)
                   / E_photon;
    if (n_sig < 0.0) n_sig = 0.0;

    /* CMB background photon count. We integrate the Planck CMB
     * radiance over the same pass band and apply the observer's
     * solid-angle acceptance (per-pixel solid angle = (λ/D)^2 ≈
     * diffraction-limited point-source FOV). For tactical-engagement
     * regimes this term is negligible; we retain it so that the
     * shot-noise denominator does not become singular at very
     * faint signals. */
    const double cmb_radiance =
        k26astro_signature_ir_planck_inband(
            1.0, 1.0, K26ASTRO_DETECT_T_CMB,
            passband_lo_m, passband_hi_m, 64) /
        K26A_PI;  /* convert back from Lambertian-Π to W/m^2/sr */
    /* Diffraction-limited solid angle ≈ (λ/D)^2 sr. */
    const double omega_pix =
        (lambda_mean_m / aperture_m) * (lambda_mean_m / aperture_m);
    const double E_bg = cmb_radiance * omega_pix;
    double n_bg = (E_bg * throughput * A_observer * integration_s)
                  / E_photon;
    if (n_bg < 0.0) n_bg = 0.0;

    /* Optics self-emission background. A telescope at temperature
     * T_optics with in-band emissivity ε_optics radiates a Planck
     * spectrum that fills the detector's diffraction-limited solid
     * angle. The collected photon rate at the focal plane is
     *
     *     N_optics = ε × L_planck_inband(T_optics) × Ω_pix
     *                × η × A_obs × τ / E_photon
     *
     * For a warm room-temperature aperture in the 3-12 μm band this
     * term dwarfs the CMB contribution by ~190 orders of magnitude
     * and sets the actual noise floor. Disabled when t_optics_k=0
     * (legacy CMB-only behaviour). */
    double n_optics = 0.0;
    if (t_optics_k > 0.0 && optics_emissivity > 0.0) {
        double eps_clamped = optics_emissivity;
        if (eps_clamped > 1.0) eps_clamped = 1.0;
        const double optics_radiance =
            k26astro_signature_ir_planck_inband(
                1.0, 1.0, t_optics_k,
                passband_lo_m, passband_hi_m, 64) /
            K26A_PI;
        const double E_optics =
            eps_clamped * optics_radiance * omega_pix;
        n_optics = (E_optics * throughput * A_observer * integration_s)
                   / E_photon;
        if (n_optics < 0.0) n_optics = 0.0;
    }

    /* Aggregate background for the noise denominator. The detection
     * event's n_background_photons field reports the total
     * background (CMB + optics) so consumers don't have to sum two
     * fields. */
    const double n_bg_total = n_bg + n_optics;

    /* Optional Poisson jitter for the signal term (background is
     * deterministically the mean — its variance is captured in
     * the shot-noise denominator). */
    if (rng) {
        const double sigma = sqrt(n_sig + n_bg_total);
        n_sig += k26c_rng_normal(rng, 0.0, sigma);
        if (n_sig < 0.0) n_sig = 0.0;
    }

    ev.n_signal_photons     = n_sig;
    ev.n_background_photons = n_bg_total;
    const double noise = sqrt(n_sig + n_bg_total + 1.0);  /* +1 dark floor */
    ev.snr      = (noise > 0.0) ? (n_sig / noise) : 0.0;
    ev.detected = (ev.snr >= snr_threshold) ? 1 : 0;
    return ev;
}
