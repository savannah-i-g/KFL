/* Round-trip + spot-projection tests for k26geo Mercator. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "k26geo.h"

static int fail = 0;

static void expect_near(double got, double want, double tol, const char *name)
{
    if (!isfinite(got) || fabs(got - want) > tol) {
        fprintf(stderr, "FAIL %s: got %.12g, want %.12g (tol %.3g)\n",
                name, got, want, tol);
        fail++;
    }
}

static void roundtrip(double lon0, double lat0, double dlat_km, double dlon_km)
{
    K26GeoOrigin o;
    K26GeoStatus s = k26geo_origin_init(&o, lon0, lat0);
    if (s != K26GEO_OK) {
        fprintf(stderr, "FAIL roundtrip init at (%.3f,%.3f): %s\n",
                lon0, lat0, k26geo_status_str(s));
        fail++;
        return;
    }

    /* Probe a point ~dlon_km east, ~dlat_km north of origin. */
    const double dlat_deg = dlat_km / 111.0;
    const double dlon_deg = dlon_km / (111.0 * cos(lat0 * M_PI / 180.0));
    const double lon = lon0 + dlon_deg;
    const double lat = lat0 + dlat_deg;

    K26V3 w = k26geo_lonlat_to_world(&o, lon, lat);

    double rlon, rlat;
    k26geo_world_to_lonlat(&o, w.x, w.z, &rlon, &rlat);

    char nameLon[128], nameLat[128];
    snprintf(nameLon, sizeof nameLon,
             "roundtrip lon @(%.3f,%.3f) +%.0fkm", lon0, lat0, dlon_km);
    snprintf(nameLat, sizeof nameLat,
             "roundtrip lat @(%.3f,%.3f) +%.0fkm", lon0, lat0, dlat_km);
    expect_near(rlon, lon, 1e-9, nameLon);
    expect_near(rlat, lat, 1e-9, nameLat);
}

int main(void)
{
    /* Origin clamp */
    K26GeoOrigin o;
    if (k26geo_origin_init(&o,  0.0,  85.0) != K26GEO_ERR_RANGE) {
        fprintf(stderr, "FAIL polar lat 85 should be rejected\n"); fail++;
    }
    if (k26geo_origin_init(&o,  0.0, -85.0) != K26GEO_ERR_RANGE) {
        fprintf(stderr, "FAIL polar lat -85 should be rejected\n"); fail++;
    }
    if (k26geo_origin_init(&o,  0.0,  70.01) != K26GEO_ERR_RANGE) {
        fprintf(stderr, "FAIL lat 70.01 should be rejected\n"); fail++;
    }
    if (k26geo_origin_init(&o,  0.0,  69.99) != K26GEO_OK) {
        fprintf(stderr, "FAIL lat 69.99 should be accepted\n"); fail++;
    }

    /* Origin maps to (0,0). */
    if (k26geo_origin_init(&o, -2.5879, 51.4544) != K26GEO_OK) {
        fprintf(stderr, "FAIL Bristol init\n"); fail++;
    } else {
        K26V3 w = k26geo_lonlat_to_world(&o, -2.5879, 51.4544);
        expect_near(w.x, 0.0, 1e-9, "origin maps to x=0");
        expect_near(w.y, 0.0, 1e-9, "origin maps to y=0");
        expect_near(w.z, 0.0, 1e-9, "origin maps to z=0");
    }

    /* East/north sign convention. */
    {
        K26GeoOrigin eq;
        k26geo_origin_init(&eq, 0.0, 0.0);
        K26V3 east  = k26geo_lonlat_to_world(&eq,  0.01, 0.00);
        K26V3 north = k26geo_lonlat_to_world(&eq,  0.00, 0.01);
        if (!(east.x > 0.0))  { fprintf(stderr, "FAIL east.x must be +ve\n"); fail++; }
        if (!(north.z < 0.0)) { fprintf(stderr, "FAIL north.z must be -ve\n"); fail++; }
    }

    /* Round-trip at several latitudes + offsets. */
    roundtrip(  0.0000,  0.0000, 1.0, 1.0);   /* Equator             */
    roundtrip( -0.1276, 51.5074, 5.0, 5.0);   /* London              */
    roundtrip(139.6917, 35.6895, 5.0, 5.0);   /* Tokyo               */
    roundtrip(-21.9426, 64.1466, 5.0, 5.0);   /* Reykjavik (mid-lat) */
    roundtrip(-58.3816,-34.6037, 5.0, 5.0);   /* Buenos Aires        */
    roundtrip( -2.5879, 51.4544, 0.05, 0.05); /* Bristol fine-grain  */

    if (fail) {
        fprintf(stderr, "%d failure(s)\n", fail);
        return 1;
    }
    puts("test_mercator: OK");
    return 0;
}
