/* counter_detect.h — active-emission counter-detection.
 *
 * Closes the loop on §detect_active_rf and §detect_active_lidar:
 * a vehicle transmitting at P_t announces its own position to
 * any passive IR observer because the eventual thermalisation of
 * transmitted power radiates as heat. The counter-detection
 * range is the distance at which a passive observer with given
 * aperture / integration / SNR threshold detects the emitter's
 * self-emitted IR signature.
 *
 * Without this surface, driver decision logic could
 * silently treat radar/lidar/jammer use as "free" — observable
 * counter-detection ranges typically exceed the active emitter's
 * own detection range by orders of magnitude.
 *
 * Source anchor: same as detect_ir_passive — Hudson 1969 + the
 * "There Ain't No Stealth In Space" synthesis (Atomic Rockets,
 * projectrho.com). */
#ifndef K26ASTRO_DETECT_COUNTER_DETECT_H
#define K26ASTRO_DETECT_COUNTER_DETECT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Compute the range at which a passive IR observer detects the
 * emitter's self-radiated heat. The emitter's transmitted power
 * P_tx_W is multiplied by a dissipation factor (see
 * K26ASTRO_DETECT_COUNTER_DETECT_DISSIPATION_FACTOR) to produce
 * an equivalent thermal load on the emitter's radiator structure;
 * the radiator's effective temperature T_K is the caller's input
 * (typically set by the platform thermal model).
 *
 * Returns the range (m) at which SNR = snr_threshold for the
 * observer's aperture + integration setup. Returns 0.0 on invalid
 * inputs.
 *
 * The relation derives from the SNR equation by setting SNR =
 * snr_threshold and solving for R; the solution scales as
 *
 *     R_counter ∝ √(P_eff · A_observer · τ)
 *
 * with P_eff = P_tx · dissipation_factor the thermalised power. */
double k26astro_counter_detect_ir_range(
    double emitter_p_tx_W,
    double emitter_radiator_T_K,
    double observer_aperture_m,
    double observer_integration_s,
    double observer_passband_lo_um,
    double observer_passband_hi_um,
    double observer_throughput,
    double snr_threshold);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_COUNTER_DETECT_H */
