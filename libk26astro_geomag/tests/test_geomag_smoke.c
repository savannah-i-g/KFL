/* test_geomag_smoke.c — gate for libk26astro_geomag.
 *
 * IGRF-14 produces physically-meaningful Earth-magnetic-field values.
 * The smoke test exercises:
 *   - Boulder, CO (mid-latitude N hemisphere) at 2025.0 — expect
 *     total field 40-65 μT, downward component positive
 *   - Sydney (S hemisphere) at 2025.0 — expect downward component
 *     NEGATIVE (field exits Earth in S hemisphere)
 *   - Reykjavik (high-latitude N) at 2025.0 — expect total field
 *     > 45 μT
 *   - K26V3 form matches the struct form bit-for-bit
 *   - Out-of-range epoch (year 2050) returns E_OUT_OF_RANGE
 *   - Determinism: two same-input calls return bit-identical results
 *
 * Reference values use broad physical-range checks rather than
 * bit-exact pinning — IGRF predictions for 2025-2030 are by design
 * non-definitive (the 5-year extrapolation isn't bit-stable across
 * IGRF generations). For bit-stable regression, future work can
 * pull values from BGS's validation table once it's vendored
 * alongside the coefficients. */

#include "k26astro_geomag/geomag.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures_ = 0;

static double deg_(double d) { return d * M_PI / 180.0; }

static void test_boulder_2025_(void)
{
    /* Boulder, CO: 40.015 N, 105.27 W, altitude 1655 m. Epoch 2025
     * (years past J2000 = 25.0). */
    K26AstroGeomagField field;
    int rc = k26astro_geomag_field_at(deg_(40.015), deg_(-105.27),
                                       1655.0, 25.0, &field);
    if (rc != K26ASTRO_GEOMAG_OK) {
        fprintf(stderr, "FAIL boulder_2025: rc=%d\n", rc);
        failures_++;
        return;
    }
    /* Total field magnitude: Boulder is well-characterised at ~52 μT
     * (52,000 nT) circa 2025; allow a broad 40-65 μT range. */
    const double F = field.F_total_T * 1e9;   /* back to nT for log */
    if (F < 40000.0 || F > 65000.0) {
        fprintf(stderr,
                "FAIL boulder_2025: F=%.1f nT outside [40000, 65000]\n", F);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS boulder_2025: F=%.1f nT in expected range\n", F);
    }
    /* Z (down) should be positive in N hemisphere. */
    if (field.B_down_T <= 0.0) {
        fprintf(stderr,
                "FAIL boulder_2025: B_down=%.6f T should be > 0 (N hemisphere)\n",
                field.B_down_T);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS boulder_2025: B_down=%.6f T positive (N hemisphere)\n",
                field.B_down_T);
    }
    /* X (north) should be positive at Boulder (field points generally
     * north of east in N hemisphere). */
    if (field.B_north_T <= 0.0) {
        fprintf(stderr,
                "FAIL boulder_2025: B_north=%.6f T should be > 0\n",
                field.B_north_T);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS boulder_2025: B_north=%.6f T positive\n",
                field.B_north_T);
    }
}

static void test_sydney_2025_(void)
{
    /* Sydney: 33.8688 S, 151.2093 E. */
    K26AstroGeomagField field;
    int rc = k26astro_geomag_field_at(deg_(-33.8688), deg_(151.2093),
                                       58.0, 25.0, &field);
    if (rc != K26ASTRO_GEOMAG_OK) {
        fprintf(stderr, "FAIL sydney_2025: rc=%d\n", rc);
        failures_++;
        return;
    }
    /* S hemisphere: Z (down) should be NEGATIVE. */
    if (field.B_down_T >= 0.0) {
        fprintf(stderr,
                "FAIL sydney_2025: B_down=%.6f T should be < 0 (S hemisphere)\n",
                field.B_down_T);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS sydney_2025: B_down=%.6f T negative (S hemisphere)\n",
                field.B_down_T);
    }
}

static void test_reykjavik_2025_(void)
{
    /* Reykjavik: 64.1466 N, 21.9426 W. High-latitude N → strong field. */
    K26AstroGeomagField field;
    int rc = k26astro_geomag_field_at(deg_(64.1466), deg_(-21.9426),
                                       0.0, 25.0, &field);
    if (rc != K26ASTRO_GEOMAG_OK) {
        fprintf(stderr, "FAIL reykjavik_2025: rc=%d\n", rc);
        failures_++;
        return;
    }
    const double F = field.F_total_T * 1e9;
    /* High-latitude N: typically 50-55 μT. */
    if (F < 45000.0 || F > 65000.0) {
        fprintf(stderr,
                "FAIL reykjavik_2025: F=%.1f nT outside [45000, 65000]\n", F);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS reykjavik_2025: F=%.1f nT in expected range\n", F);
    }
}

static void test_v3_matches_struct_(void)
{
    K26AstroGeomagField field;
    k26astro_geomag_field_at(deg_(40.015), deg_(-105.27), 1655.0, 25.0,
                              &field);
    K26V3 v = k26astro_geomag_field_v3(deg_(40.015), deg_(-105.27),
                                        1655.0, 25.0);
    if (v.x != field.B_north_T || v.y != field.B_east_T ||
        v.z != field.B_down_T) {
        fprintf(stderr,
                "FAIL v3_matches_struct: v=[%.17g %.17g %.17g] vs struct=[%.17g %.17g %.17g]\n",
                v.x, v.y, v.z,
                field.B_north_T, field.B_east_T, field.B_down_T);
        failures_++;
    } else {
        fprintf(stderr, "PASS v3_matches_struct: bit-identical\n");
    }
}

static void test_out_of_range_epoch_(void)
{
    /* Year 2050 is far beyond IGRF-14's [1900, 2035] range. */
    K26AstroGeomagField field;
    int rc = k26astro_geomag_field_at(deg_(40.0), deg_(0.0), 0.0, 50.0,
                                       &field);
    if (rc != K26ASTRO_GEOMAG_E_OUT_OF_RANGE) {
        fprintf(stderr,
                "FAIL out_of_range_epoch: rc=%d (expected E_OUT_OF_RANGE=%d)\n",
                rc, K26ASTRO_GEOMAG_E_OUT_OF_RANGE);
        failures_++;
        return;
    }
    /* Out-of-range should produce zero output. */
    if (field.B_north_T != 0.0 || field.B_east_T != 0.0 ||
        field.B_down_T != 0.0 || field.F_total_T != 0.0) {
        fprintf(stderr,
                "FAIL out_of_range_epoch: nonzero output despite error\n");
        failures_++;
    } else {
        fprintf(stderr,
                "PASS out_of_range_epoch: E_OUT_OF_RANGE + zero output\n");
    }
}

static void test_magnitude_matches_components_(void)
{
    K26AstroGeomagField field;
    k26astro_geomag_field_at(deg_(40.015), deg_(-105.27), 1655.0, 25.0,
                              &field);
    const double computed = sqrt(field.B_north_T * field.B_north_T +
                                  field.B_east_T  * field.B_east_T  +
                                  field.B_down_T  * field.B_down_T);
    const double err = fabs(computed - field.F_total_T) / field.F_total_T;
    if (err > 1.0e-6) {
        fprintf(stderr,
                "FAIL magnitude_matches: computed %.6f vs F %.6f (rel err %.3g)\n",
                computed, field.F_total_T, err);
        failures_++;
    } else {
        fprintf(stderr,
                "PASS magnitude_matches: sqrt(N²+E²+D²) = F to %.3g rel\n", err);
    }
}

static void test_determinism_(void)
{
    K26AstroGeomagField a, b;
    k26astro_geomag_field_at(deg_(40.015), deg_(-105.27), 1655.0, 25.0,
                              &a);
    k26astro_geomag_field_at(deg_(40.015), deg_(-105.27), 1655.0, 25.0,
                              &b);
    if (a.B_north_T != b.B_north_T ||
        a.B_east_T  != b.B_east_T  ||
        a.B_down_T  != b.B_down_T  ||
        a.F_total_T != b.F_total_T) {
        fprintf(stderr,
                "FAIL determinism: not bit-identical across runs\n");
        failures_++;
    } else {
        fprintf(stderr,
                "PASS determinism: byte-identical across two runs\n");
    }
}

int main(void)
{
    fprintf(stderr, "libk26astro_geomag smoke test\n");

    test_boulder_2025_();
    test_sydney_2025_();
    test_reykjavik_2025_();
    test_v3_matches_struct_();
    test_out_of_range_epoch_();
    test_magnitude_matches_components_();
    test_determinism_();

    if (failures_) {
        fprintf(stderr, "TOTAL: %d FAIL\n", failures_);
        return 1;
    }
    fprintf(stderr, "TOTAL: all PASS\n");
    return 0;
}
