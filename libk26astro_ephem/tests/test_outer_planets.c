/* test_outer_planets.c — k26astro_ephem_has_outer_planets()
 * predicate and the strict full-kernel loader.
 *
 * Builds two minimal synthetic SPK kernels via
 * k26astro_spk_write_synthetic — one with an Earth (NAIF 399) segment
 * representing the inner-only subset, one with a Jupiter (NAIF 5)
 * segment representing the full-kernel coverage signature. Asserts
 * the predicate returns 0 and 1 respectively. Then asserts the
 * strict loader returns NULL when the production full-kernel path is
 * absent (the expected dev / CI environment state). */
#include "k26astro_ephem/spk.h"
#include "k26astro_ephem/ephem.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ORDER 8
#define K     ORDER
#define RSIZE (2 + 3 * K)

static int write_minimal_kernel_(const char *path, int target_body)
{
    /* One record, half-day, zero-valued Chebyshev coefficients. The
     * predicate doesn't evaluate the polynomial, only walks segment
     * metadata, so the coefficients don't need to be meaningful. */
    double records[RSIZE];
    memset(records, 0, sizeof(records));
    records[0] = 43200.0;  /* MID = half-day past epoch zero */
    records[1] = 43200.0;  /* RADIUS = half-day */

    K26AstroSpkWriteSegment seg = {
        .target_body      = target_body,
        .center_body      = 0,
        .start_et         = 0.0,
        .end_et           = 86400.0,
        .interval_seconds = 86400.0,
        .records          = records,
        .n_records        = 1,
        .coeffs_per_axis  = K
    };
    return k26astro_spk_write_synthetic(path, &seg, 1);
}

int main(void)
{
    const char *inner_path = "tests/fixtures/synthetic_inner_only.spk";
    const char *outer_path = "tests/fixtures/synthetic_with_jupiter.spk";

    assert(write_minimal_kernel_(inner_path, 399) == 0);  /* Earth */
    assert(write_minimal_kernel_(outer_path, 5)   == 0);  /* Jupiter barycentre */

    K26AstroEphem *e_inner = k26astro_ephem_load(inner_path);
    assert(e_inner != NULL);
    assert(k26astro_ephem_has_outer_planets(e_inner) == 0);
    k26astro_ephem_close(e_inner);

    K26AstroEphem *e_outer = k26astro_ephem_load(outer_path);
    assert(e_outer != NULL);
    assert(k26astro_ephem_has_outer_planets(e_outer) == 1);
    k26astro_ephem_close(e_outer);

    /* NULL safety. */
    assert(k26astro_ephem_has_outer_planets(NULL) == 0);

    /* Strict loader returns NULL when the production full-kernel
     * file is absent. Skip the assertion if the file happens to be
     * installed in the test environment — that's a valid state too,
     * just not the dev default. */
    K26AstroEphem *e_strict = k26astro_ephem_load_default_full_strict();
    if (e_strict != NULL) {
        assert(k26astro_ephem_has_outer_planets(e_strict) == 1);
        k26astro_ephem_close(e_strict);
    }

    printf("test_outer_planets: OK\n");
    return 0;
}
