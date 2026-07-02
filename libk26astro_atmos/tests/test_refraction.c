/* test_refraction.c — Bennett 1982 refraction gate.
 *
 * Verifies k26astro_atmos_apparent against the Bennett 1982 closed-
 * form at five reference elevations:
 *
 *   0°  → ~34.5'  (horizon plateau)
 *   10° → ~5.3'
 *   30° → ~1.7'
 *   45° → ~1.0'
 *   90° → 0'      (zenith — no refraction)
 *
 * Test approach: build the Earth standard atmosphere, set up a
 * geometric direction at the requested elevation, call
 * _apparent, and compare the resulting elevation against
 * (true_elev + R_bennett). Tolerance is 0.5 arcsec to allow for
 * the n0 scale factor at the Earth-standard's exact n₀.
 *
 * Acceptance:
 *   - elev=0° apparent lift ≈ 34.5 arcmin ± 1 arcmin
 *   - elev=10° apparent lift ≈ 5.3 arcmin ± 0.5 arcmin
 *   - elev=45° apparent lift ≈ 1.0 arcmin ± 0.5 arcmin
 *   - elev=90° apparent lift = 0 (zenith pass-through) */
#include "k26astro_atmos/atmos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double elevation_(K26V3 dir, K26V3 zenith)
{
    double dn = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    double zn = sqrt(zenith.x*zenith.x + zenith.y*zenith.y + zenith.z*zenith.z);
    double cs = (dir.x*zenith.x + dir.y*zenith.y + dir.z*zenith.z) / (dn * zn);
    if (cs >  1.0) cs =  1.0;
    if (cs < -1.0) cs = -1.0;
    return asin(cs);
}

static int check_(K26AstroAtmos *a, double true_deg,
                   double expected_arcmin, double tol_arcmin)
{
    K26V3 zenith = { 0.0, 0.0, 1.0 };

    /* Geometric direction at the requested elevation in the
     * x-z plane (azimuth = 0, elevation = true_deg). */
    double true_rad = true_deg * (M_PI / 180.0);
    K26V3 geom = { cos(true_rad), 0.0, sin(true_rad) };

    K26V3 app = k26astro_atmos_apparent(a, geom, zenith);
    double app_elev = elevation_(app, zenith);
    double lift_rad = app_elev - true_rad;
    double lift_arcmin = lift_rad * (180.0 * 60.0 / M_PI);
    double diff = fabs(lift_arcmin - expected_arcmin);

    fprintf(stderr,
        "  true=%5.1f° lift=%6.3f arcmin (expected ~%.3f, "
        "tol %.3f) %s\n",
        true_deg, lift_arcmin, expected_arcmin, tol_arcmin,
        (diff <= tol_arcmin) ? "OK" : "FAIL");
    return (diff <= tol_arcmin) ? 0 : 1;
}

int main(void)
{
    K26AstroAtmos *a = k26astro_atmos_earth_standard();
    assert(a);

    int fail = 0;
    fail |= check_(a,  0.0, 34.5, 2.0);   /* horizon plateau */
    fail |= check_(a, 10.0,  5.30, 0.5);
    fail |= check_(a, 30.0,  1.70, 0.3);
    fail |= check_(a, 45.0,  1.00, 0.3);
    fail |= check_(a, 90.0,  0.00, 0.01); /* zenith */

    k26astro_atmos_destroy(a);

    if (fail) {
        fprintf(stderr, "test_refraction: FAIL\n");
        return 1;
    }
    fprintf(stderr, "test_refraction: OK\n");
    return 0;
}
