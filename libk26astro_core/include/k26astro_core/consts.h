/* k26astro_core/consts.h — physical + astronomical constants.
 *
 * Sources:
 *   - IAU 2015 Resolution B3 nominal solar + planetary conversion
 *     constants (defined as exact SI conversion factors so all
 *     downstream computation uses a single canonical reference).
 *   - CODATA 2018 fundamental constants (G, c, h, k_B).
 *
 * Naming: SI base unit assumed (metres, kilograms, seconds, kelvin)
 * unless suffixed. All values are double-precision; long-double is
 * not used so the determinism contract holds across platforms whose
 * `long double` mantissa width varies (x86_64 has 80-bit extended,
 * aarch64 has 128-bit IEEE-754 quad — both diverge from the binary64
 * the rest of the library uses).
 *
 * The IAU Resolution B3 nominal values are the canonical source for
 * GM_sun / GM_earth etc. Earlier "best-known" values from CODATA
 * have larger uncertainties and shift between releases; the IAU
 * nominal values are *defined* exactly so they don't drift. */
#ifndef K26ASTRO_CORE_CONSTS_H
#define K26ASTRO_CORE_CONSTS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Mathematical constants ------------------------------------- */
#define K26A_PI               3.141592653589793238462643383279502884
#define K26A_TWO_PI           6.283185307179586476925286766559005768
#define K26A_HALF_PI          1.570796326794896619231321691639751442
#define K26A_DEG_PER_RAD      (180.0 / K26A_PI)
#define K26A_RAD_PER_DEG      (K26A_PI / 180.0)
#define K26A_ARCSEC_PER_RAD   (3600.0 * K26A_DEG_PER_RAD)
#define K26A_RAD_PER_ARCSEC   (K26A_RAD_PER_DEG / 3600.0)

/* ---- CODATA 2018 fundamentals (SI, exact or "best") ------------- */
/* Speed of light in vacuum, defining constant (exact). */
#define K26A_C                299792458.0
/* Newtonian constant of gravitation. CODATA 2018 recommended; ~10⁻⁵
 * relative uncertainty — the largest of the major fundamentals. */
#define K26A_G                6.67430e-11
/* Planck constant, defining constant (exact). */
#define K26A_H_PLANCK         6.62607015e-34
/* Boltzmann constant, defining constant (exact). */
#define K26A_K_BOLTZMANN      1.380649e-23

/* ---- IAU 2015 Resolution B3: nominal solar conversion constants - */
/* Defined exactly. Use these in GM products to avoid drift across
 * CODATA updates. */
#define K26A_GM_SUN           1.3271244e20       /* m^3 / s^2 — heliocentric (IAU 2015 B3 nominal, exact) */
#define K26A_R_SUN_EQU        6.957e8            /* m — solar equatorial radius */
#define K26A_L_SUN            3.828e26           /* W — total radiative output  */
#define K26A_T_SUN_EFF        5772.0             /* K — effective temperature   */

/* IAU 2015 nominal terrestrial conversion constants (geocentric). */
#define K26A_GM_EARTH         3.986004e14        /* m^3 / s^2 */
#define K26A_R_EARTH_EQU      6.3781e6           /* m — equatorial radius */
#define K26A_R_EARTH_POL      6.3568e6           /* m — polar radius */
#define K26A_EARTH_FLATTENING (1.0 / 298.257223563) /* WGS84 */
/* Sidereal day (relative to vernal equinox). */
#define K26A_EARTH_SIDEREAL_DAY_S 86164.0905    /* s */
/* Solar day (relative to Sun, mean). */
#define K26A_EARTH_SOLAR_DAY_S    86400.0       /* s — definition */

/* IAU 2015 nominal Jovicentric conversion constant. */
#define K26A_GM_JUPITER       1.2668653e17       /* m^3 / s^2 */

/* ---- Astronomical units of distance ----------------------------- */
/* IAU 2012 redefinition: 1 AU = 149_597_870_700 m exactly. */
#define K26A_AU_M             149597870700.0
/* Parsec — defined from AU + arcsec definition; rounded to the bit. */
#define K26A_PC_M             (K26A_AU_M / K26A_RAD_PER_ARCSEC)
#define K26A_LY_M             9.4607304725808e15  /* IAU 2015: light-year exact */
#define K26A_LIGHT_SEC_M      K26A_C              /* trivially */

/* ---- Time constants --------------------------------------------- */
/* Julian year (365.25 d × 86400 s). The astronomical year. */
#define K26A_YEAR_JULIAN_S    (365.25 * 86400.0)
/* TT - TAI offset (Resolution B1.9, IAU 2000). Exact. */
#define K26A_TT_TAI_OFFSET_S  32.184
/* GMST sidereal-day length used by simple-frame ECI→ECEF mappings. */
#define K26A_SIDEREAL_DAY_SECONDS  86164.0905

/* ---- Selected planetary masses (informational) ----------------- */
/* M_planet / M_sun. From the GM ratios published with DE441; useful
 * when the user wants quick mass-ratio math without pulling in the
 * full ephemeris. */
#define K26A_M_MERCURY_RATIO  (1.66012082e-7)
#define K26A_M_VENUS_RATIO    (2.44783824e-6)
#define K26A_M_EARTH_RATIO    (3.00348962e-6)
#define K26A_M_MARS_RATIO     (3.22715608e-7)
#define K26A_M_JUPITER_RATIO  (9.54791915e-4)
#define K26A_M_SATURN_RATIO   (2.85885670e-4)
#define K26A_M_URANUS_RATIO   (4.36624374e-5)
#define K26A_M_NEPTUNE_RATIO  (5.15138902e-5)

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CORE_CONSTS_H */
