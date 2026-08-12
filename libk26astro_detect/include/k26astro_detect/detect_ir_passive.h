/* detect_ir_passive.h — passive infrared detection geometry.
 *
 * An observer with optical aperture D and dwell time τ collects
 * photons from a target emitter against the cosmic-microwave-
 * background floor (T_CMB ≈ 2.725 K). Detection range scales as
 *
 *     R_d ∝ √A · T^2
 *
 * for an emitter with projected radiating area A and effective
 * temperature T (Hudson 1969; "There Ain't No Stealth In Space", Atomic Rockets). The
 * exact range depends on observer aperture, integration time,
 * pass-band, SNR threshold, and optics throughput.
 *
 * The detection model is photon-shot-noise dominated:
 *
 *     SNR = N_signal / √(N_signal + N_background + N_dark)
 *
 * with N_signal and N_background derived from the in-band Planck
 * integral over the observer's pass-band. The lib assumes:
 *   - Greybody emission (caller supplies effective temperature
 *     and emissivity-weighted area).
 *   - Lambertian surface emission (the projected area in the
 *     observer's line of sight equals the geometric area times
 *     cos θ; caller may supply an effective area that already
 *     includes the projection).
 *   - Zero atmospheric attenuation (vacuum engagement bubble).
 *
 * Stated limitations:
 *   - Spectrally-structured plumes are not modelled — caller
 *     supplies a greybody equivalent.
 *   - Dark-current noise is per-system; the default value is a
 *     conservative cryo-cooled HgCdTe baseline.
 *   - Aperture diffraction limits the angular resolution but
 *     does not bound the integrated flux measurement directly;
 *     consumers wanting to model point-source-vs-extended-target
 *     discrimination should track angular size separately. */
#ifndef K26ASTRO_DETECT_IR_PASSIVE_H
#define K26ASTRO_DETECT_IR_PASSIVE_H

#include "k26compute.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Detection event record. */
typedef struct {
    double range_m;                 /* observer-target distance */
    double flux_received_W_per_m2;  /* in-band irradiance at observer */
    double snr;                     /* signal-to-noise */
    double n_signal_photons;        /* expected in-band photon count */
    double n_background_photons;    /* expected CMB-background photons */
    int    detected;                /* 1 if snr >= snr_threshold */
} K26AstroDetectIrEvent;

/* Compute a passive IR detection event.
 *
 * Inputs:
 *   emitter_power_W    Total in-band radiated power (W) — caller
 *                      derives via k26astro_signature_ir_planck_inband
 *                      from the target's surface state.
 *   emitter_T_K        Effective emitter temperature (used for the
 *                      noise-floor model, not the signal model).
 *   range_m            Observer-target geometric distance.
 *   aperture_m         Observer optical aperture diameter.
 *   integration_s      Dwell / integration time.
 *   passband_lo_um     Lower pass-band edge (μm).
 *   passband_hi_um     Upper pass-band edge (μm).
 *   throughput         Optical throughput (transmissivity × QE), [0, 1].
 *   snr_threshold      Detection-trigger SNR threshold.
 *   rng                Optional RNG for noise jitter (NULL for
 *                      deterministic mean-only output).
 *
 * Background model: CMB Planck radiance over the pass-band, summed
 * over the diffraction-limited pixel solid angle. For warm-aperture
 * instruments the optics' own thermal emission dominates the CMB
 * floor by many orders of magnitude in the 3-12 μm band — prefer
 * the _with_optics variant for any room-temperature instrument
 * looking in its own emission band. */
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
    K26CRng *rng);

/* Compute a passive IR detection event with telescope self-emission
 * folded into the background.
 *
 * Adds two parameters to k26astro_detect_ir_passive:
 *
 *   t_optics_k         Telescope (optics + baffle + dewar window)
 *                      effective temperature (K). Pass 0 to disable
 *                      optics self-emission (recovers the
 *                      CMB-only behaviour of the legacy function).
 *                      Typical values:
 *                        290 K — uncooled room-temperature telescope
 *                        280 K — minimally-shaded warm baffle
 *                        80-120 K — cryo-cooled cold stop
 *                        4-30 K — space-cryogenic (JWST-class)
 *   optics_emissivity  In-band emissivity of the optics path (0..1).
 *                      Typical:
 *                        0.05 — well-baffled cryo, low-emissivity coatings
 *                        0.30 — bare aluminium primary
 *                        1.0 — black-painted absorber surface
 *
 * Background photon count is augmented by the optics self-emission:
 *
 *     N_optics = (B_optics × η × A_obs × Ω × τ) / E_photon
 *
 * where B_optics = ε × L_Planck_inband(T_optics) is the optics
 * in-band radiance and Ω is the diffraction-limited pixel solid
 * angle. The shot-noise denominator becomes
 * √(N_signal + N_CMB + N_optics + 1).
 *
 * For a 280 K telescope with ε=0.05 looking in 3-12 μm with a 1 m
 * aperture, N_optics dominates N_CMB by ~190 orders of magnitude
 * and sets the actual noise floor — without it the model is
 * background-free and effectively any non-zero signal triggers
 * detection at any SNR threshold. This is why the legacy
 * CMB-only function should not be relied on for warm-aperture
 * instruments. */
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
    K26CRng *rng);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_IR_PASSIVE_H */
