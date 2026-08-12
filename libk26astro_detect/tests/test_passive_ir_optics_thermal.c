/* test_passive_ir_optics_thermal.c — telescope-self-emission gate.
 *
 * Demonstrates that adding an optics thermal self-emission term to
 * the passive IR detection model establishes the realistic noise
 * floor for a warm-aperture instrument looking in its own emission
 * band — the binding constraint that was silently set by the
 * detector noise-figure +1 dark floor in the CMB-only model.
 *
 * Three regimes:
 *   1. CMB-only (legacy): SNR is effectively unbounded for any
 *      non-trivial in-band signal because n_bg ≈ 0.
 *   2. Cold optics (T_optics=80 K, ε=0.05, cryo-baffled cold stop):
 *      SNR remains very high — cold telescopes are essentially
 *      background-free against warm targets.
 *   3. Warm optics (T_optics=280 K, ε=0.05): the optics in-band
 *      emission dominates noise; SNR drops dramatically and the
 *      detection regime becomes "warm vs warm" with realistic
 *      contrast limits.
 *
 * Scenario: a 300 K target with ~100 m² effective radiating area
 * observed at 5×10⁵ m through a 1 m aperture with 0.5 s
 * integration in the 3-12 μm pass-band. */
#include "k26astro_detect/detect_ir_passive.h"
#include "k26astro_detect/signature.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* Target signature: 300 K greybody, ε_target=0.9, A=100 m². */
    const double T_target_K   = 300.0;
    const double area_m2      = 100.0;
    const double eps_target   = 0.9;
    const double bolo_total_W = eps_target * 5.670374419e-8 * area_m2
                              * T_target_K * T_target_K
                              * T_target_K * T_target_K;

    const double range_m        = 5.0e5;     /* 500 km — engagement bubble */
    const double aperture_m     = 1.0;
    const double integration_s  = 0.5;
    const double passband_lo_um = 3.0;
    const double passband_hi_um = 12.0;
    const double throughput     = 0.7;
    const double snr_threshold  = 5.0;

    /* ---- Regime 1: CMB-only (legacy entry point) -------------- */
    K26AstroDetectIrEvent ev_cmb_only =
        k26astro_detect_ir_passive(
            bolo_total_W, T_target_K,
            range_m, aperture_m, integration_s,
            passband_lo_um, passband_hi_um,
            throughput, snr_threshold,
            /* rng */ NULL);

    /* Sanity: emitter is detected (signal far exceeds CMB shot
     * noise). */
    assert(ev_cmb_only.detected == 1);
    fprintf(stderr,
            "test_passive_ir_optics_thermal: CMB-only — "
            "n_sig=%.3e n_bg=%.3e SNR=%.3e\n",
            ev_cmb_only.n_signal_photons,
            ev_cmb_only.n_background_photons,
            ev_cmb_only.snr);

    /* ---- Regime 2: Cold cryo-baffled optics (T_optics=80 K) -- */
    K26AstroDetectIrEvent ev_cold =
        k26astro_detect_ir_passive_with_optics(
            bolo_total_W, T_target_K,
            range_m, aperture_m, integration_s,
            passband_lo_um, passband_hi_um,
            throughput, snr_threshold,
            /* t_optics_k       */ 80.0,
            /* optics_emissivity*/ 0.05,
            NULL);

    /* Cold optics still detects the 300 K target — cryo-baffling
     * is *the* technique that makes IR astronomy work. */
    assert(ev_cold.detected == 1);
    /* The added optics noise should be small but non-zero,
     * making n_bg_total strictly larger than the CMB-only case. */
    assert(ev_cold.n_background_photons >
           ev_cmb_only.n_background_photons);
    fprintf(stderr,
            "test_passive_ir_optics_thermal: 80 K cryo — "
            "n_sig=%.3e n_bg=%.3e SNR=%.3e\n",
            ev_cold.n_signal_photons,
            ev_cold.n_background_photons,
            ev_cold.snr);

    /* ---- Regime 3: Warm room-temperature optics (T_optics=280 K) */
    K26AstroDetectIrEvent ev_warm =
        k26astro_detect_ir_passive_with_optics(
            bolo_total_W, T_target_K,
            range_m, aperture_m, integration_s,
            passband_lo_um, passband_hi_um,
            throughput, snr_threshold,
            /* t_optics_k       */ 280.0,
            /* optics_emissivity*/ 0.05,
            NULL);

    /* Warm optics inject a substantial background photon count.
     * The 280 K optics radiate Planck-spectrum 3-12 μm flux into
     * the diffraction-limited pixel, and at the small (λ/D)²
     * solid angle this collects to a finite-but-significant rate. */
    assert(ev_warm.n_background_photons >
           1.0e6 * ev_cmb_only.n_background_photons);
    fprintf(stderr,
            "test_passive_ir_optics_thermal: 280 K warm — "
            "n_sig=%.3e n_bg=%.3e SNR=%.3e detected=%d\n",
            ev_warm.n_signal_photons,
            ev_warm.n_background_photons,
            ev_warm.snr,
            ev_warm.detected);

    /* SNR ordering: cold > warm, monotone in optics temperature. */
    assert(ev_cold.snr > ev_warm.snr);

    /* ---- Disabling optics yields legacy behaviour --------------*
     *
     * t_optics_k=0 or optics_emissivity=0 should recover the
     * CMB-only result bit-identically. */
    K26AstroDetectIrEvent ev_disabled =
        k26astro_detect_ir_passive_with_optics(
            bolo_total_W, T_target_K,
            range_m, aperture_m, integration_s,
            passband_lo_um, passband_hi_um,
            throughput, snr_threshold,
            /* t_optics_k       */ 0.0,
            /* optics_emissivity*/ 0.05,
            NULL);
    assert(ev_disabled.n_background_photons ==
           ev_cmb_only.n_background_photons);
    assert(ev_disabled.snr == ev_cmb_only.snr);

    K26AstroDetectIrEvent ev_zero_emiss =
        k26astro_detect_ir_passive_with_optics(
            bolo_total_W, T_target_K,
            range_m, aperture_m, integration_s,
            passband_lo_um, passband_hi_um,
            throughput, snr_threshold,
            /* t_optics_k       */ 280.0,
            /* optics_emissivity*/ 0.0,
            NULL);
    assert(ev_zero_emiss.n_background_photons ==
           ev_cmb_only.n_background_photons);

    fprintf(stderr, "test_passive_ir_optics_thermal: OK\n");
    return 0;
}
