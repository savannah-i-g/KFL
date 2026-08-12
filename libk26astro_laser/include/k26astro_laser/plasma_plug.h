/* plasma_plug.h — ablation-plume self-attenuation.
 *
 * Once the target's ablation plume reaches inverse-bremsstrahlung
 * absorption density, subsequent photons are partially absorbed
 * in transit, attenuating the beam before it reaches the surface.
 * The effect is well-documented in long-pulse laser-target
 * experiments and is the dominant degradation mechanism at high
 * fluence (Phipps 2010 §IV.C).
 *
 * The empirical anchor for the attenuation factor is weak —
 * published values span a factor of ~3 depending on plume
 * geometry, ambient gas, and pulse-shape details. The lib exposes
 * the attenuation as a first-order configurable factor; the
 * caller's fit table supplies the per-material constant.
 *
 * Model: T_plasma = exp(-k · max(0, Φ/Φ_th - 1))
 *
 * where T_plasma is the plume transmissivity (in [0, 1]) and k is
 * the per-material attenuation strength. T_plasma → 1 at Φ_th,
 * → e^-k at Φ = 2 Φ_th; the default k = 0.5 reflects the
 * published 30%+ attenuation observed at 2-3x threshold. */
#ifndef K26ASTRO_LASER_PLASMA_PLUG_H
#define K26ASTRO_LASER_PLASMA_PLUG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Default attenuation strength factor (dimensionless). Caller's
 * fit table may override per-material. */
#define K26ASTRO_LASER_PLASMA_PLUG_DEFAULT_K  0.5

/* Plume transmissivity for incident fluence Φ relative to the
 * threshold Φ_th, with attenuation strength k. Returns 1.0
 * below threshold; otherwise exp(-k * (Φ/Φ_th - 1)). */
double k26astro_plasma_plug_transmissivity(double fluence_J_per_m2,
                                            double threshold_J_per_m2,
                                            double k_attenuation);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_LASER_PLASMA_PLUG_H */
