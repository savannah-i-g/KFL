/* k26astro_geomag.c — C-side wrapper for libk26astro_geomag.
 *
 * Converts K26's natural coordinate system (geodetic lat/lon in
 * radians, altitude in metres, epoch in years past J2000.0) into
 * IGRF-14's native units (colatitude/east-longitude in degrees,
 * altitude in km, AD-year date), invokes the Fortran synthesis
 * routine, then converts the result from nT to tesla. */

#include "k26astro_geomag/geomag.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ----- Fortran bind(C) entries (defined in k26astro_geomag_iface.f90) */

extern void k26astro_geomag_field_call(double date,
                                       double colat_deg, double elong_deg,
                                       double alt_km,
                                       double *out_x_nT, double *out_y_nT,
                                       double *out_z_nT, double *out_f_nT);

extern void k26astro_geomag_secvar_call(double date,
                                        double colat_deg, double elong_deg,
                                        double alt_km,
                                        double *out_x_nT, double *out_y_nT,
                                        double *out_z_nT, double *out_f_nT);

/* ----- Coordinate transforms --------------------------------------- */

static double k_wrap_lon_deg_(double lon_deg)
{
    /* Wrap into [0, 360). */
    double w = fmod(lon_deg, 360.0);
    if (w < 0.0) w += 360.0;
    return w;
}

/* Convert K26's altitude-above-WGS84 (metres) → IGRF's geodetic
 * altitude (km). Identity scaling, just unit conversion. */
static double k_altitude_m_to_km_(double altitude_m)
{
    return altitude_m * 1.0e-3;
}

/* Validate inputs common to both field_at and secular_variation_at.
 * Returns 0 if OK, non-zero K26ASTRO_GEOMAG_E_* otherwise. */
static int k_validate_(double lat_rad, double lon_rad, double altitude_m,
                       double epoch_y2000_yr, void *out)
{
    (void)lon_rad;
    if (!out) return K26ASTRO_GEOMAG_E_BAD_INPUT;
    if (lat_rad < -M_PI * 0.5 - 1.0e-9 ||
        lat_rad >  M_PI * 0.5 + 1.0e-9) {
        return K26ASTRO_GEOMAG_E_BAD_INPUT;
    }
    /* Roughly Earth-radius sanity: don't allow altitude inside the
     * core (-6378 km is the surface, anything below ~-3500 km would
     * be inside the Earth). IGRF itself only validates altitude>3485
     * km in geocentric mode (itype=2); we use itype=1 (geodetic)
     * which is altitude above WGS84 surface, so negative-but-shallow
     * values (e.g., depth below sea level) are physically meaningful
     * for ocean / underground use cases. Clamp at -10 km below sea
     * level. */
    if (altitude_m < -10000.0) {
        return K26ASTRO_GEOMAG_E_BAD_INPUT;
    }
    /* Epoch range check (IGRF-14 valid 1900.0–2030.0, definitive
     * 1945.0–2020.0, predictive thereafter; warning above 2030.0,
     * hard refuse above 2035.0). */
    const double ad_year = k26astro_geomag_y2000_to_ad_(epoch_y2000_yr);
    if (ad_year < K26ASTRO_GEOMAG_EPOCH_MIN_Y ||
        ad_year > K26ASTRO_GEOMAG_EPOCH_MAX_Y) {
        return K26ASTRO_GEOMAG_E_OUT_OF_RANGE;
    }
    return K26ASTRO_GEOMAG_OK;
}

/* ----- C-direct struct-out API ------------------------------------- */

int k26astro_geomag_field_at(double lat_rad, double lon_rad,
                             double altitude_m,
                             double epoch_y2000_yr,
                             K26AstroGeomagField *out)
{
    const int rc = k_validate_(lat_rad, lon_rad, altitude_m,
                                epoch_y2000_yr, out);
    if (rc != K26ASTRO_GEOMAG_OK) {
        if (out) {
            out->B_north_T = 0.0;
            out->B_east_T  = 0.0;
            out->B_down_T  = 0.0;
            out->F_total_T = 0.0;
        }
        return rc;
    }

    const double date      = k26astro_geomag_y2000_to_ad_(epoch_y2000_yr);
    const double colat_deg = 90.0 - (lat_rad * 180.0 / M_PI);
    const double elong_deg = k_wrap_lon_deg_(lon_rad * 180.0 / M_PI);
    const double alt_km    = k_altitude_m_to_km_(altitude_m);

    double x_nT, y_nT, z_nT, f_nT;
    k26astro_geomag_field_call(date, colat_deg, elong_deg, alt_km,
                                &x_nT, &y_nT, &z_nT, &f_nT);

    /* nT → T (1 nT = 1e-9 T). */
    out->B_north_T = x_nT * 1.0e-9;
    out->B_east_T  = y_nT * 1.0e-9;
    out->B_down_T  = z_nT * 1.0e-9;
    out->F_total_T = f_nT * 1.0e-9;
    return K26ASTRO_GEOMAG_OK;
}

int k26astro_geomag_secular_variation_at(double lat_rad, double lon_rad,
                                         double altitude_m,
                                         double epoch_y2000_yr,
                                         K26AstroGeomagField *out)
{
    const int rc = k_validate_(lat_rad, lon_rad, altitude_m,
                                epoch_y2000_yr, out);
    if (rc != K26ASTRO_GEOMAG_OK) {
        if (out) {
            out->B_north_T = 0.0;
            out->B_east_T  = 0.0;
            out->B_down_T  = 0.0;
            out->F_total_T = 0.0;
        }
        return rc;
    }

    const double date      = k26astro_geomag_y2000_to_ad_(epoch_y2000_yr);
    const double colat_deg = 90.0 - (lat_rad * 180.0 / M_PI);
    const double elong_deg = k_wrap_lon_deg_(lon_rad * 180.0 / M_PI);
    const double alt_km    = k_altitude_m_to_km_(altitude_m);

    double x_nT, y_nT, z_nT, f_nT;
    k26astro_geomag_secvar_call(date, colat_deg, elong_deg, alt_km,
                                 &x_nT, &y_nT, &z_nT, &f_nT);

    /* nT/year → T/year. IGRF in isv=1 returns "rubbish" for f; the
     * caller recomputes total from the magnitudes. */
    out->B_north_T = x_nT * 1.0e-9;
    out->B_east_T  = y_nT * 1.0e-9;
    out->B_down_T  = z_nT * 1.0e-9;
    out->F_total_T = sqrt(out->B_north_T * out->B_north_T +
                          out->B_east_T  * out->B_east_T  +
                          out->B_down_T  * out->B_down_T);
    return K26ASTRO_GEOMAG_OK;
}

/* ----- KFL-callable surface (K26V3 return) ------------------------- */

K26V3 k26astro_geomag_field_v3(double lat_rad, double lon_rad,
                               double altitude_m,
                               double epoch_y2000_yr)
{
    K26AstroGeomagField field;
    const int rc = k26astro_geomag_field_at(lat_rad, lon_rad, altitude_m,
                                             epoch_y2000_yr, &field);
    K26V3 v = { 0.0, 0.0, 0.0 };
    if (rc == K26ASTRO_GEOMAG_OK) {
        v.x = field.B_north_T;
        v.y = field.B_east_T;
        v.z = field.B_down_T;
    }
    return v;
}

double k26astro_geomag_field_magnitude(double lat_rad, double lon_rad,
                                       double altitude_m,
                                       double epoch_y2000_yr)
{
    K26AstroGeomagField field;
    if (k26astro_geomag_field_at(lat_rad, lon_rad, altitude_m,
                                  epoch_y2000_yr, &field)
        != K26ASTRO_GEOMAG_OK) {
        return 0.0;
    }
    return field.F_total_T;
}
