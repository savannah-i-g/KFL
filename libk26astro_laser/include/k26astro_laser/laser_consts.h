/* laser_consts.h — diffraction, coupling, and material defaults.
 *
 * Phipps momentum-coupling coefficients are tabulated for common
 * target materials at canonical wavelengths (Phipps 1989 + 2010).
 * The table is small by design — callers configuring a specific
 * vehicle should override the per-material values from their
 * mission fit tables. Defaults are representative; the
 * actual coupling is a function of wavelength, fluence,
 * intensity, and surface state, so per-engagement overrides are
 * the norm in published analyses.
 *
 * References:
 *   - Phipps, C.R. et al. 1988. Impulse coupling to targets in
 *     vacuum by KrF, HF, and CO2 single-pulse lasers.
 *     J. Appl. Phys. 64:1083.
 *   - Phipps, C.R. et al. 2010. Review: Laser-ablation propulsion.
 *     J. Propulsion and Power 26:609.
 *   - Born, M. and Wolf, E. 2002. Principles of Optics, 7th ed.
 *     §8.5 (Airy diffraction).
 */
#ifndef K26ASTRO_LASER_CONSTS_H
#define K26ASTRO_LASER_CONSTS_H

/* Airy first-null half-angle factor: 1.22. */
#define K26ASTRO_LASER_AIRY_FACTOR              1.22

/* Encircled-energy fraction inside the Airy first dark ring
 * (Born & Wolf 2002 §8.5.2). The remaining 16.2% spills into
 * the surrounding rings. */
#define K26ASTRO_LASER_AIRY_ENCIRCLED_CENTRAL   0.8378

/* Default beam quality M^2 for an unspecified emitter (diffraction-
 * limited single-mode). High-quality beam directors run 1.1–1.5;
 * commercial industrial fibre lasers 2–10. */
#define K26ASTRO_LASER_DEFAULT_M_SQUARED        1.2

/* Default RMS pointing jitter for an unspecified mount, rad. */
#define K26ASTRO_LASER_DEFAULT_JITTER_RAD       1.0e-6

/* Default mirror wavefront RMS error in waves; Maréchal Strehl
 * stays > 0.8 for σ_w < λ/14. */
#define K26ASTRO_LASER_DEFAULT_WAVEFRONT_WAVES  (1.0 / 20.0)

/* Phipps coupling material identifiers. Mission fit tables are
 * keyed by these enum values; new materials extend the enum but
 * never renumber existing entries. */
typedef enum {
    K26ASTRO_LASER_MAT_NONE        = 0,
    K26ASTRO_LASER_MAT_ALUMINUM    = 1,
    K26ASTRO_LASER_MAT_STEEL       = 2,
    K26ASTRO_LASER_MAT_TITANIUM    = 3,
    K26ASTRO_LASER_MAT_COPPER      = 4,
    K26ASTRO_LASER_MAT_COMPOSITE   = 5,   /* generic CFRP composite */
    K26ASTRO_LASER_MAT_FUSED_SILICA = 6
} K26AstroLaserMaterial;

/* Default Phipps momentum-coupling coefficient C_m (N/W) at
 * representative IR / near-IR wavelengths. Tabulated values are
 * order-of-magnitude anchors at fluences above the plasma
 * ignition threshold; below threshold, vapour-pressure coupling
 * applies and C_m is roughly 10x smaller. Mission fits override
 * per-engagement.
 *
 * Source: Phipps 1989 Table II + Phipps 2010 §III.A. Values are
 * for normal incidence in vacuum; off-axis incidence reduces C_m
 * by cos θ.
 *
 * Units: N · W^-1 (one newton of impulse per watt of optical
 * power coupled into the target). */
#define K26ASTRO_LASER_CM_DEFAULT_AL_N_PER_W           5.0e-5
#define K26ASTRO_LASER_CM_DEFAULT_STEEL_N_PER_W        4.0e-5
#define K26ASTRO_LASER_CM_DEFAULT_TITANIUM_N_PER_W     3.5e-5
#define K26ASTRO_LASER_CM_DEFAULT_COPPER_N_PER_W       2.5e-5
#define K26ASTRO_LASER_CM_DEFAULT_COMPOSITE_N_PER_W    8.0e-5
#define K26ASTRO_LASER_CM_DEFAULT_FUSED_SILICA_N_PER_W 1.5e-5

/* Default plasma-ignition fluence threshold (J · m^-2) at 1064 nm.
 * Phipps 2010 §III.B (page 616) gives I·τ^(1/2) = constant for the
 * plasma-ignition intensity, i.e. pulse-width (not wavelength)
 * scaling. The lib returns these 1064 nm values for all
 * wavelengths; per-wavelength + per-pulse-width tables override. */
#define K26ASTRO_LASER_THRESHOLD_AL_J_PER_M2_AT_1064NM       1.0e8
#define K26ASTRO_LASER_THRESHOLD_STEEL_J_PER_M2_AT_1064NM    1.5e8
#define K26ASTRO_LASER_THRESHOLD_TITANIUM_J_PER_M2_AT_1064NM 1.2e8
#define K26ASTRO_LASER_THRESHOLD_COPPER_J_PER_M2_AT_1064NM   2.0e8
#define K26ASTRO_LASER_THRESHOLD_COMPOSITE_J_PER_M2_AT_1064NM 5.0e7
#define K26ASTRO_LASER_THRESHOLD_FUSED_SILICA_J_PER_M2_AT_1064NM 3.0e8

/* Default specific-ablation energy Q* (J / kg) — energy required
 * to ablate unit target mass once the plasma plug is established.
 * Phipps 2010 §III.A gives Q* ~ 1e7-1e8 J/kg for metals; the
 * coupling efficiency η = C_m c / 2 implies a self-consistent
 * Q* derivable from C_m + jet velocity; the lib uses the
 * tabulated Q* directly to avoid double-counting. */
#define K26ASTRO_LASER_Q_STAR_AL_J_PER_KG            2.5e7
#define K26ASTRO_LASER_Q_STAR_STEEL_J_PER_KG         3.5e7
#define K26ASTRO_LASER_Q_STAR_TITANIUM_J_PER_KG      4.0e7
#define K26ASTRO_LASER_Q_STAR_COPPER_J_PER_KG        2.0e7
#define K26ASTRO_LASER_Q_STAR_COMPOSITE_J_PER_KG     1.5e7
#define K26ASTRO_LASER_Q_STAR_FUSED_SILICA_J_PER_KG  3.0e7

#endif /* K26ASTRO_LASER_CONSTS_H */
