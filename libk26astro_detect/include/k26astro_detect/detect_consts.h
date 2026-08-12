/* detect_consts.h — radiometric and detection-physics constants.
 *
 * The Stefan-Boltzmann constant is derived from K26A_H_PLANCK,
 * K26A_K_BOLTZMANN, and K26A_C per the SI definition:
 *
 *     σ = (2 π^5 k_B^4) / (15 h^3 c^2)
 *
 * with CODATA-2019 values for h, k_B, c (defined-exact in SI).
 *
 * The cosmic microwave background temperature is the canonical
 * Planck 2018 result (Fixsen 2009 + Planck 2018 cosmological
 * parameters; 2.7255 ± 0.0006 K — we use the rounded 2.725 K
 * value standard in space-object IR detection literature).
 *
 * Detection-side defaults (SNR threshold, optics throughput) are
 * conservative values used by the AIAA / SPIE space-object IR
 * detection literature (Hudson 1969 IR system engineering; the
 * "There Ain't No Stealth In Space" synthesis on the Atomic
 * Rockets site, projectrho.com). Consumers should override
 * these for specific sensor specifications.
 */
#ifndef K26ASTRO_DETECT_CONSTS_H
#define K26ASTRO_DETECT_CONSTS_H

/* Stefan-Boltzmann constant: W · m^-2 · K^-4. Derived; provided
 * here because libk26astro_core does not currently export σ
 * directly. */
#define K26ASTRO_DETECT_SIGMA_SB        5.670374419e-8

/* Cosmic microwave background temperature, K. Fixsen 2009 /
 * Planck 2018. */
#define K26ASTRO_DETECT_T_CMB           2.725

/* Default optical throughput for an unspecified observer aperture
 * (transmissivity × quantum efficiency); replaced by per-sensor
 * spec if known. Conservative 0.5 covers most practical IR
 * systems (Hudson 1969, Table 4-3). */
#define K26ASTRO_DETECT_DEFAULT_THROUGHPUT  0.5

/* Default IR pass-band edges, μm. The 3-12 μm window covers the
 * mid-wave / long-wave IR regimes where blackbody emission from
 * ~100-1000 K bodies dominates and atmospheric / radiator
 * emission peaks. */
#define K26ASTRO_DETECT_DEFAULT_PASSBAND_LO_UM   3.0
#define K26ASTRO_DETECT_DEFAULT_PASSBAND_HI_UM  12.0

/* Default detection SNR threshold (dimensionless). 5σ over noise
 * floor is the standard space-object detection convention. */
#define K26ASTRO_DETECT_DEFAULT_SNR_THRESH       5.0

/* Active-emission self-dissipation factor for counter-detection.
 * A radar emitter's transmitted power is eventually thermalised
 * by the platform's radiator + structure; the default 1.0
 * assumes every watt transmitted radiates as heat. Calibrated
 * platforms can override with a smaller factor (some power
 * leaves as radio-frequency that does NOT thermalise). */
#define K26ASTRO_DETECT_COUNTER_DETECT_DISSIPATION_FACTOR  1.0

/* Photon energy at wavelength λ (m): E = h c / λ. */
#define K26ASTRO_DETECT_PHOTON_ENERGY_AT(lambda_m)  \
    (K26A_H_PLANCK * K26A_C / (lambda_m))

#endif /* K26ASTRO_DETECT_CONSTS_H */
