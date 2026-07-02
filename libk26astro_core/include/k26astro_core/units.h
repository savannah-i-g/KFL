/* k26astro_core/units.h — SI unit conversions.
 *
 * All internal computation is SI (metres, seconds, kilograms, kelvin).
 * The conversion helpers here exist for the API surface — KFL source
 * authors writing `let r: double = au(1.5)` get a cleaner read than
 * `let r: double = 1.5 * K26A_AU_M`, and downstream lib boundaries
 * (KFL builtin manifests, JSON exports of K26AstroPos) document
 * the unit explicitly.
 *
 * Inline `static` so callers compile against constant-folded
 * literals. The .c file exists for symmetry (the linker pulls
 * nothing from it today) and for future non-inlineable helpers
 * (e.g. catalogue-magnitude conversions with bandpass tables). */
#ifndef K26ASTRO_CORE_UNITS_H
#define K26ASTRO_CORE_UNITS_H

#include "consts.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Length ----------------------------------------------------- */
static inline double k26a_au(double v)            { return v * K26A_AU_M; }
static inline double k26a_pc(double v)            { return v * K26A_PC_M; }
static inline double k26a_ly(double v)            { return v * K26A_LY_M; }
static inline double k26a_light_seconds(double v) { return v * K26A_LIGHT_SEC_M; }
static inline double k26a_km(double v)            { return v * 1.0e3; }

/* Inverse conversions: metres back to nominal unit. */
static inline double k26a_to_au(double m)            { return m / K26A_AU_M; }
static inline double k26a_to_pc(double m)            { return m / K26A_PC_M; }
static inline double k26a_to_ly(double m)            { return m / K26A_LY_M; }
static inline double k26a_to_light_seconds(double m) { return m / K26A_LIGHT_SEC_M; }
static inline double k26a_to_km(double m)            { return m * 1.0e-3; }

/* ---- Time ------------------------------------------------------- */
static inline double k26a_year(double v) { return v * K26A_YEAR_JULIAN_S; }
static inline double k26a_day (double v) { return v * 86400.0; }
static inline double k26a_hr  (double v) { return v * 3600.0; }
static inline double k26a_min (double v) { return v * 60.0; }

static inline double k26a_to_year(double s) { return s / K26A_YEAR_JULIAN_S; }
static inline double k26a_to_day (double s) { return s / 86400.0; }
static inline double k26a_to_hr  (double s) { return s / 3600.0; }
static inline double k26a_to_min (double s) { return s / 60.0; }

/* ---- Angle ------------------------------------------------------ */
static inline double k26a_deg     (double v) { return v * K26A_RAD_PER_DEG; }
static inline double k26a_arcmin  (double v) { return v * K26A_RAD_PER_DEG / 60.0; }
static inline double k26a_arcsec  (double v) { return v * K26A_RAD_PER_ARCSEC; }
static inline double k26a_mas     (double v) { return v * K26A_RAD_PER_ARCSEC * 1.0e-3; }

static inline double k26a_to_deg    (double r) { return r * K26A_DEG_PER_RAD; }
static inline double k26a_to_arcmin (double r) { return r * K26A_DEG_PER_RAD * 60.0; }
static inline double k26a_to_arcsec (double r) { return r * K26A_ARCSEC_PER_RAD; }
static inline double k26a_to_mas    (double r) { return r * K26A_ARCSEC_PER_RAD * 1.0e3; }

/* ---- Mass ------------------------------------------------------- *
 * Note: solar mass is `M_sun = GM_sun / G`. The IAU nominal GM_sun
 * is exact, but G has CODATA uncertainty — so M_sun in kg is *not*
 * exact. Use GM products directly whenever possible (gravity
 * computations don't need M, they need GM). */
static inline double k26a_msun(double v) {
    return v * K26A_GM_SUN / K26A_G;
}
static inline double k26a_to_msun(double kg) {
    return kg * K26A_G / K26A_GM_SUN;
}
/* GM in solar units — useful when KFL source wants
 * `let gm: double = msun_gm(1.0)` for the Sun. */
static inline double k26a_msun_gm(double m_in_sun_units) {
    return m_in_sun_units * K26A_GM_SUN;
}

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CORE_UNITS_H */
