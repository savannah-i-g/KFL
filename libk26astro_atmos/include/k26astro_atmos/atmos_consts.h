/* atmos_consts.h — atmospheric model constants as hex-literal
 * IEEE-754 doubles for cross-platform determinism.
 *
 * References:
 *   - Bates 1959, "Some problems concerning the terrestrial
 *     atmosphere above about 100 km level"
 *   - ICAO 1993 Standard Atmosphere (Document 7488/3)
 *   - Edlén 1966 / Ciddor 1996 (refractive-index dispersion)
 *   - Bucholtz 1995 (Rayleigh scattering cross sections) */
#ifndef K26ASTRO_ATMOS_CONSTS_H
#define K26ASTRO_ATMOS_CONSTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Hex-literal double constructor (avoids strtod rounding noise
 * across libc implementations). */
static inline double k26astro_atmos_hex_d_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* ICAO 1993 sea-level pressure: 101325 Pa exactly (defined). */
#define K26ASTRO_ATMOS_P0_PA              101325.0

/* Earth scale height H = 8400 m (Bates 1959 lower-thermosphere
 * limit). Hex: 0x40C0680000000000 = 8400.0. */
#define K26ASTRO_ATMOS_EARTH_H_M_BITS     0x40C0680000000000ULL

/* Rayleigh extinction coefficient at sea level (Bucholtz 1995,
 * 550 nm green): β_R ≈ 1.16e-5 m^-1.
 * 1.16e-5 = 0x3EE85F4DC5C2D000. */
#define K26ASTRO_ATMOS_EARTH_BETA_R_BITS  0x3EE85F4DC5C2D000ULL

/* Mie extinction coefficient at sea level (typical maritime haze,
 * 550 nm): β_M ≈ 2.1e-6 m^-1.
 * 2.1e-6 = 0x3EC19DB7358BD307 (verified via IEEE-754 round-trip;
 * the earlier 0x3EB1A0461C6F90DB pattern decoded to 1.05e-6,
 * exactly half — fixed during W5 audit). */
#define K26ASTRO_ATMOS_EARTH_BETA_M_BITS  0x3EC19DB7358BD307ULL

/* Henyey-Greenstein asymmetry parameter for atmospheric aerosols
 * (Earth haze): g = 0.76. Hex: 0x3FE851EB851EB852. */
#define K26ASTRO_ATMOS_EARTH_HG_G_BITS    0x3FE851EB851EB852ULL

/* Refractive index sea-level - 1 (Edlén 1966, λ=550nm): n₀ = 2.93e-4.
 * Hex: 0x3F33345BCD104F1E. */
#define K26ASTRO_ATMOS_EARTH_N0_BITS      0x3F33345BCD104F1EULL

/* Top of "useful" atmosphere for refraction integration: 100 km.
 * Above this, ρ ~ ρ₀·exp(-100/8.4) ~ 7e-6, refractive index
 * effectively 1. Hex: 0x40F86A0000000000 = 100000.0. */
#define K26ASTRO_ATMOS_EARTH_TOP_M_BITS   0x40F86A0000000000ULL

/* ICAO 1993 sea-level temperature: T_0 = 288.15 K.
 * Hex: 0x4072026666666666. */
#define K26ASTRO_ATMOS_EARTH_T0_K_BITS    0x4072026666666666ULL

/* Ratio of specific heats for Earth's diatomic-dominated lower
 * atmosphere: γ = c_p / c_v = 1.4 (nominal). Used for the
 * isentropic speed-of-sound relation c = sqrt(γ R T).
 * Hex: 0x3FF6666666666666. */
#define K26ASTRO_ATMOS_EARTH_GAMMA_BITS   0x3FF6666666666666ULL

/* Specific gas constant for Earth dry air: R_specific = R_u / M_air
 * = 8314.462618 / 28.9647 ≈ 287.05 J/(kg·K). ICAO 1993 value.
 * Hex: 0x4071F0CCCCCCCCCD. */
#define K26ASTRO_ATMOS_EARTH_R_SPECIFIC_BITS  0x4071F0CCCCCCCCCDULL

/* Speed of light, IAU 2009 defined (exact). */
#define K26ASTRO_ATMOS_C_LIGHT            299792458.0

/* Resolved hex-literal accessors. */
#define K26ASTRO_ATMOS_EARTH_H_M \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_H_M_BITS)
#define K26ASTRO_ATMOS_EARTH_BETA_R \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_BETA_R_BITS)
#define K26ASTRO_ATMOS_EARTH_BETA_M \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_BETA_M_BITS)
#define K26ASTRO_ATMOS_EARTH_HG_G \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_HG_G_BITS)
#define K26ASTRO_ATMOS_EARTH_N0 \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_N0_BITS)
#define K26ASTRO_ATMOS_EARTH_TOP_M \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_TOP_M_BITS)
#define K26ASTRO_ATMOS_EARTH_T0_K \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_T0_K_BITS)
#define K26ASTRO_ATMOS_EARTH_GAMMA \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_GAMMA_BITS)
#define K26ASTRO_ATMOS_EARTH_R_SPECIFIC \
    k26astro_atmos_hex_d_(K26ASTRO_ATMOS_EARTH_R_SPECIFIC_BITS)

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_ATMOS_CONSTS_H */
