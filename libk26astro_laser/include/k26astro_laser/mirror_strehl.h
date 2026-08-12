/* mirror_strehl.h — wavefront-error Strehl degradation.
 *
 * Maréchal approximation: for a mirror system with RMS wavefront
 * error σ_w (in waves at wavelength λ), the on-axis Strehl ratio is
 *
 *     S = exp( -(2π σ_w_waves)^2 )
 *
 * The approximation holds for σ_w < ~λ/10; beyond that the actual
 * Strehl drops faster than the formula predicts and consumers
 * should consult diffraction-integral computations.
 *
 * In thermal-loaded mirrors the wavefront error scales with
 * internal flux (W / m^2) above a per-system threshold; the lib
 * does not model the thermal-distortion path directly — the
 * caller's fit supplies the current σ_w state. */
#ifndef K26ASTRO_LASER_MIRROR_STREHL_H
#define K26ASTRO_LASER_MIRROR_STREHL_H

#ifdef __cplusplus
extern "C" {
#endif

/* On-axis Strehl ratio for RMS wavefront error σ_w (m) at
 * wavelength λ (m). σ_w / λ < 0.1 is the regime where Maréchal
 * is accurate; outside that range the formula remains
 * monotonically decreasing but should be treated as an
 * order-of-magnitude estimate. */
double k26astro_mirror_strehl(double rms_wavefront_m,
                              double wavelength_m);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_LASER_MIRROR_STREHL_H */
