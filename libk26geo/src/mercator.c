/* Web Mercator tangent-plane projection.
 *
 * We use spherical Mercator (R = 6378137 m, the WGS84 equatorial radius)
 * because the eccentricity correction is well below our 1 m world unit
 * at city scales — and we never claim sub-metre geodetic accuracy. The
 * extra cos(lat0) factor in `scale` linearises the projection at the
 * loaded data's origin latitude so a 1 km × 1 km neighbourhood spans
 * (approximately) 1000 world units in both X and Z.
 *
 * Forward:
 *     mx = (lon - lon0) * π/180 * R * cos(lat0)
 *     my = (mercY(lat) - mercY(lat0)) * R * cos(lat0)
 *
 *     mercY(φ) = log(tan(π/4 + φ·π/360))
 *
 * World axes (matches libk26m3d right-handed +X right / +Y up / -Z fwd):
 *     world.x =  +mx     (east)
 *     world.z =  -my     (so +north goes into the screen)
 *     world.y =   0      (altitude — caller sets)
 *
 * Inverse:
 *     lon = lon0 + (world.x / (R · cos(lat0))) · 180/π
 *     lat = (2·atan(exp(merc_y0 + (-world.z)/(R·cos(lat0)))) - π/2) · 180/π
 */

#include <math.h>

#include "k26geo.h"

double k26geo_merc_y(double lat_deg)
{
    return log(tan(M_PI / 4.0 + lat_deg * M_PI / 360.0));
}

K26GeoStatus k26geo_origin_init(K26GeoOrigin *out, double lon_deg, double lat_deg)
{
    if (!out)                                          return K26GEO_ERR_INVAL;
    if (!isfinite(lon_deg) || !isfinite(lat_deg))      return K26GEO_ERR_INVAL;
    if (fabs(lat_deg) > K26GEO_MAX_ABS_LAT_DEG)        return K26GEO_ERR_RANGE;

    out->lon_deg  = lon_deg;
    out->lat_deg  = lat_deg;
    out->merc_y0  = k26geo_merc_y(lat_deg);
    out->cos_lat0 = cos(lat_deg * M_PI / 180.0);
    out->scale    = K26GEO_EARTH_RADIUS_M * out->cos_lat0;
    return K26GEO_OK;
}

K26V3 k26geo_lonlat_to_world(const K26GeoOrigin *o, double lon_deg, double lat_deg)
{
    K26V3 v = { 0, 0, 0 };
    if (!o) return v;

    const double dlon_rad = (lon_deg - o->lon_deg) * (M_PI / 180.0);
    const double dy_merc  = k26geo_merc_y(lat_deg) - o->merc_y0;

    v.x =  dlon_rad * o->scale;
    v.y =  0.0;
    v.z = -dy_merc  * o->scale;
    return v;
}

void k26geo_world_to_lonlat(const K26GeoOrigin *o,
                            double world_x, double world_z,
                            double *out_lon_deg, double *out_lat_deg)
{
    if (!o) {
        if (out_lon_deg) *out_lon_deg = 0.0;
        if (out_lat_deg) *out_lat_deg = 0.0;
        return;
    }

    const double dlon_rad = world_x / o->scale;
    const double dy_merc  = -world_z / o->scale;

    double lon = o->lon_deg + dlon_rad * (180.0 / M_PI);
    if (lon >  180.0) lon =  180.0;
    if (lon < -180.0) lon = -180.0;

    const double merc_y = o->merc_y0 + dy_merc;
    double lat = (2.0 * atan(exp(merc_y)) - M_PI / 2.0) * (180.0 / M_PI);
    if (lat >  85.0) lat =  85.0;
    if (lat < -85.0) lat = -85.0;

    if (out_lon_deg) *out_lon_deg = lon;
    if (out_lat_deg) *out_lat_deg = lat;
}
