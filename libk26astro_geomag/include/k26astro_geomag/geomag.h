/* libk26astro_geomag - Earth magnetic field via IGRF-14.
 *
 * IGRF-14 (NOAA NCEI / IAGA Working Group V-MOD, released 2024-11)
 * is the canonical reference Earth-magnetic-field model. K26 uses
 * the synthesis routine `igrf14syn` to evaluate the field at any
 * (lat, lon, altitude, epoch) within the coefficient validity range
 * [1900.0, 2030.0]. Fortran-link library, alongside libk26astro_quad
 * and libk26astro_ode.
 *
 * Consumer use cases:
 *   - Lorentz-force perturbation in LEO/charged-particle dynamics
 *   - Ionospheric coupling for libk26astro_render atmospheric extensions
 *   - Magnetic-coordinate transforms (apex / quasi-dipole / IGRF
 *     -based MLT)
 *   - Compass-deviation correction for navigation tools
 *
 * Determinism contract: pure function — `field_at(lat, lon, alt, t)`
 * is bit-stable given the same inputs. -O2 -ffp-contract=off, no
 * thread state, no allocations in the synthesis path. */
#ifndef K26ASTRO_GEOMAG_H
#define K26ASTRO_GEOMAG_H

#include "k26m3d.h"  /* K26V3 */

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_GEOMAG_LIB_VERSION "0.1.0"

/* ----- Return codes ------------------------------------------------- */

#define K26ASTRO_GEOMAG_OK              0
#define K26ASTRO_GEOMAG_E_BAD_INPUT    -1
#define K26ASTRO_GEOMAG_E_OUT_OF_RANGE -2   /* date outside [1900,2035] */

/* Full magnetic-field state at a point. Components are in the local
 * geographic ENU frame (X=north, Y=east, Z=down — IGRF convention).
 * F_total_T is the total field magnitude (Euclidean norm of N/E/D).
 * Units: tesla (1 T = 1e4 gauss). IGRF natively returns nT; the
 * wrapper converts to T for SI consistency with the rest of K26. */
typedef struct K26AstroGeomagField {
    double B_north_T;
    double B_east_T;
    double B_down_T;
    double F_total_T;
} K26AstroGeomagField;

/* Coefficient validity range (Julian years AD). Outside this range
 * IGRF-14 prints a warning and returns f=1e8 + zero components. */
#define K26ASTRO_GEOMAG_EPOCH_MIN_Y   1900.0
#define K26ASTRO_GEOMAG_EPOCH_MAX_Y   2035.0

/* Epoch convention: K26 expresses time as years past J2000.0
 * (2000-01-01 12:00:00 TT). Convert to IGRF's AD-year via
 * `igrf_date = 2000.0 + epoch_y2000_yr`. */
static inline double k26astro_geomag_y2000_to_ad_(double epoch_y2000_yr)
{
    return 2000.0 + epoch_y2000_yr;
}

/* ----- C-direct API (struct-out, full info) ------------------------- */

/* Evaluate the IGRF-14 main field at (lat, lon, altitude, epoch).
 *
 * Parameters:
 *   lat_rad        : geodetic latitude in radians, range [-π/2, π/2]
 *   lon_rad        : geodetic longitude in radians (east positive),
 *                    will be wrapped into [0, 2π)
 *   altitude_m     : height above WGS84 ellipsoid in metres
 *                    (must be > -6378137 to avoid going inside Earth)
 *   epoch_y2000_yr : time in years past J2000.0 (so 2025 ≈ +25)
 *   out            : output struct (must be non-NULL)
 *
 * Returns K26ASTRO_GEOMAG_OK on success;
 *         K26ASTRO_GEOMAG_E_OUT_OF_RANGE if the date is outside
 *         the [1900, 2035] coefficient range;
 *         K26ASTRO_GEOMAG_E_BAD_INPUT for NULL out or invalid range. */
int k26astro_geomag_field_at(double lat_rad, double lon_rad,
                             double altitude_m,
                             double epoch_y2000_yr,
                             K26AstroGeomagField *out);

/* Secular variation (dB/dt at the requested epoch). Same arg
 * convention as field_at; output components are tesla per year. */
int k26astro_geomag_secular_variation_at(double lat_rad, double lon_rad,
                                         double altitude_m,
                                         double epoch_y2000_yr,
                                         K26AstroGeomagField *out);

/* ----- KFL-callable surface (K26V3 return) -------------------------- */

/* IGRF main field at (lat, lon, alt, epoch) as a K26V3 in the local
 * ENU frame (x=north, y=east, z=down), units of tesla. If the date
 * is out of range, returns the zero vector.
 *
 * This is the primary KFL-callable entry point — KFL handles K26V3
 * returns directly (same convention as libk26astro_atmos_apparent). */
K26V3 k26astro_geomag_field_v3(double lat_rad, double lon_rad,
                               double altitude_m,
                               double epoch_y2000_yr);

/* Total field magnitude at (lat, lon, alt, epoch). Tesla.
 * Convenience for KFL programs that only need |B|. */
double k26astro_geomag_field_magnitude(double lat_rad, double lon_rad,
                                       double altitude_m,
                                       double epoch_y2000_yr);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_GEOMAG_H */
