/* libk26astro_atmos — planetary atmospheric model.
 *
 * Single library covering BOTH:
 *   1. Refraction (deterministic) — used by libk26astro_rt's
 *      Topocentric observer mode to bend the apparent direction
 *      of a target body relative to the geometric direction.
 *   2. Single-scatter Rayleigh + Henyey-Greenstein Mie scattering
 *      (rendering) — consumed by libk26astro_render's scatter
 *      dispatch to produce qualitative sky-glow / horizon
 *      gradients.
 *
 * One model object (K26AstroAtmos) holds the parameters shared by
 * both code paths so a Topocentric observation and its rendered
 * sky stay consistent about where the horizon is.
 *
 * Atmospheric model: barometric density profile
 *   ρ(h) = ρ₀ · exp(-h / H)
 * with H = scale height (~8400 m for Earth). For thermospheric
 * structure and drag-grade density, see the NRLMSISE-00 wrapper
 * below.
 *
 * Refraction: Snell-bent ray, 64-step trapezoidal integration of
 *   dθ/ds = -(1/n) · dn/dh · sin(zenith)
 * along the ray path. Halts when the ray exits the atmosphere
 * (h > 100 km). Bennett 1982 closed-form is the gate reference.
 *
 * Single-scatter: integrates phase × extinction along view ray.
 * No shadow integration, no multiple scattering. LUT-based
 * multi-scattering (Hillaire 2020) is not currently provided. */
#ifndef K26ASTRO_ATMOS_H
#define K26ASTRO_ATMOS_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_ATMOS_LIB_VERSION "0.1.0"

typedef struct K26AstroAtmosParams {
    double scale_height_m;       /* H in barometric ρ(h) = ρ₀·exp(-h/H) */
    double rayleigh_coeff;       /* β_R at sea level (m^-1) */
    double mie_coeff;            /* β_M at sea level (m^-1) */
    double ground_pressure_pa;   /* P₀; informational + drag-future */
    double mie_asymmetry_g;      /* HG g, in (-1, +1); 0.76 for Earth haze */
    double n0;                   /* refractive index at sea level - 1
                                   * (~2.93e-4 for Earth STP). */
    double atmos_top_m;          /* upper bound for refraction
                                   * integration (~100 km for Earth) */
    /* Thermodynamic state for the isothermal speed-of-sound relation
     * c = sqrt(γ · R_specific · T). The barometric density profile is
     * isothermal (T independent of altitude); these fields supply the
     * single representative T₀ and the gas-property pair. Per-planet
     * defaults are populated by k26astro_atmos_earth_standard. Custom
     * AtmosParams that leave these zero will have
     * k26astro_atmos_sound_speed_at return 0. */
    double t_0_k;                        /* sea-level temperature, K */
    double gamma_ratio_specific_heats;   /* γ = c_p / c_v */
    double r_specific_j_per_kg_k;        /* specific gas constant R/M */
} K26AstroAtmosParams;

typedef struct K26AstroAtmos K26AstroAtmos;

/* ---- Lifecycle -------------------------------------------------- */

/* Returns NULL on bad params (any required field non-positive). */
K26AstroAtmos *k26astro_atmos_init(const K26AstroAtmosParams *p);
void           k26astro_atmos_destroy(K26AstroAtmos *a);

/* Standard-atmosphere convenience: Earth's ICAO 1993 standard
 * with Bates 1959 scale height. */
K26AstroAtmos *k26astro_atmos_earth_standard(void);

/* ---- Refraction (deterministic; consumed by libk26astro_rt) ---- */

/* Returns the apparent direction with elevation refraction
 * applied. Inputs:
 *   geometric_dir : unit vector from observer to target before
 *                   atmospheric refraction (computed by the
 *                   light-time / aberration pipeline upstream)
 *   zenith_dir    : unit vector pointing from observer
 *                   straight up (away from planet centre)
 *
 * Returns a unit vector along the apparent (refraction-bent)
 * direction of the target — geometric_dir rotated toward zenith
 * by the elevation-dependent refraction angle. Atmospheric
 * refraction LIFTS the apparent direction: a target on the
 * geometric horizon appears about 35 arcmin above on Earth.
 *
 * For targets already at or above zenith, returns geometric_dir
 * unchanged. For targets below the geometric horizon, applies
 * the refraction until the apparent elevation reaches 0
 * (apparent horizon), then clamps. */
K26V3 k26astro_atmos_apparent(const K26AstroAtmos *a,
                               K26V3 geometric_dir,
                               K26V3 zenith_dir);

/* Diagnostic: refraction magnitude in radians for a given true
 * elevation. Bennett 1982 closed-form scaled by n₀/n₀_earth. */
double k26astro_atmos_refraction_rad(const K26AstroAtmos *a,
                                      double true_elevation_rad);

/* ---- Single-scatter rendering (consumed by libk26astro_render) - */

/* In-scattered radiance along a view ray. view_origin_m is the
 * ray origin in metres relative to the planet centre; view_dir
 * is the unit ray direction; sun_dir is the unit vector toward
 * the Sun. Returns spectral radiance approximation as (R, G, B)
 * triple in arbitrary linear units (caller renormalises to
 * display range). */
K26V3 k26astro_atmos_inscatter(const K26AstroAtmos *a,
                                K26V3 view_origin_m,
                                K26V3 view_dir,
                                K26V3 sun_dir);

/* ---- Diagnostics ------------------------------------------------ */

/* Density at altitude h above sea level. */
double k26astro_atmos_density_at(const K26AstroAtmos *a, double h_m);

/* Speed of sound at altitude h. The barometric profile is isothermal
 * (T(h) = T_0), so the returned value is altitude-independent within
 * the barometric region; it nonetheless takes `h_m` so callers can
 * substitute a non-isothermal model later without changing the
 * call site.
 *
 *   c(h) = sqrt(γ · R_specific · T_0)
 *
 * Returns 0 if any of T_0, γ, R_specific is non-positive (e.g. on a
 * zero-initialised AtmosParams that didn't populate the
 * thermodynamic fields) or if a is NULL. */
double k26astro_atmos_sound_speed_at(const K26AstroAtmos *a,
                                     double h_m);

/* Read-back of params (NULL atmos → all zeros). */
K26AstroAtmosParams k26astro_atmos_params(const K26AstroAtmos *a);

/* ---- NRLMSISE-00 thermosphere ---------------------------------- */

/* Drag-bearing density inputs. F10.7 (solar 10.7 cm flux) and ap
 * (geomagnetic indices, daily + 3-hourly) drive thermospheric
 * heating + density. Day-of-year + UT seconds set the diurnal
 * model phase.
 *
 * `ap` must point at 7 doubles (daily ap + six 3-hourly samples
 * spanning the previous 24 hours). Sources: NOAA SWPC daily
 * archives, ESA Space Weather, or any IAGA Bulletin distribution. */
typedef struct K26AstroAtmosNrlmsiseInputs {
    int    year;           /* e.g. 2026 */
    int    day_of_year;    /* 1..366 */
    double ut_seconds;     /* seconds past 00:00 UT */
    double alt_km;
    double latitude_deg;
    double longitude_deg;
    double local_solar_time_hr;
    double f107_daily;     /* yesterday's F10.7 flux, SFU */
    double f107_81day_avg; /* 81-day running mean */
    double ap[7];          /* daily ap + 3-hourly samples */
} K26AstroAtmosNrlmsiseInputs;

typedef struct K26AstroAtmosNrlmsiseDensities {
    double total_kg_m3;
    double n_he_m3;        /* number density, per m³ */
    double n_o_m3;
    double n_n2_m3;
    double n_o2_m3;
    double n_ar_m3;
    double n_h_m3;
    double n_n_m3;
    double t_exo_k;        /* exospheric temperature */
    double t_alt_k;        /* temperature at altitude */
} K26AstroAtmosNrlmsiseDensities;

/* Query NRLMSISE-00 thermospheric density + composition. Returns
 * 0 on success and populates `out`. The barometric profile
 * (k26astro_atmos_density_at) remains the default for simple
 * altitude-density queries; this entry point is the high-fidelity
 * path consumed by orbital-decay / drag work.
 *
 * Backed by the NRL FORTRAN distribution of NRLMSISE-00 (Picone,
 * Hedin, Drob, Aikin 2002; US-government public domain), vendored
 * at libk26astro_atmos/src/upstream/nrlmsise00/. */
#define K26ASTRO_ATMOS_OK                    0
#define K26ASTRO_ATMOS_E_NOT_IMPLEMENTED    -7
int k26astro_atmos_density_nrlmsise(
    const K26AstroAtmosNrlmsiseInputs *in,
    K26AstroAtmosNrlmsiseDensities    *out);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_ATMOS_H */
