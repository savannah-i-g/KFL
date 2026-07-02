/* test_lambert_tof_curve.c — Lambert TOF-curve diagnostic gate.
 *
 * Loads the reference CSV at tests/data/lambert_tof_reference.csv
 * (single-branch Lagrange form, ~216 rows over λ ∈ [-0.8, 0.8],
 * x ∈ [-0.9, 0.9], n_rev ∈ [0, 3]) and compares each row against
 * lambert_multi.c's tof_curve_ output (exposed via the
 * k26astro_lambert_tof_for_test diagnostic hook).
 *
 * The evaluator splits Battin (near x=1) and Lagrange (away
 * from x=1) branches. Both should be mathematically equivalent;
 * a suspected scaling mismatch at the branch boundary causes the
 * n_rev>=1 multi-branch solver to converge to the wrong x, and
 * this gate is the per-row diagnostic for closing that gap.
 *
 * Gate strategy:
 *
 *   |x| ≤ 0.7  →  K26 must agree with truth to 1e-12 relative
 *                 (well inside the Lagrange branch where K26
 *                 mathematically matches the truth derivation)
 *
 *   0.7 < |x| < 1.0 →  K26 must agree to 1e-8 relative (precision
 *                      loss near singularities allowed)
 *
 * Failure mode dump: any row that misses prints its (λ, x, n_rev,
 * K26_T, truth_T, abs_err, rel_err) so the bisection can begin from
 * a precise miss point.
 *
 * Build: standalone host C, links against libk26astro_conics.a and
 * libk26m3d.a. Run from the libk26astro_conics directory so the
 * relative CSV path resolves correctly. */
#define _GNU_SOURCE
#include "k26astro_conics/lambert_multi.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CSV_PATH "tests/data/lambert_tof_reference.csv"

#define TOL_LAGRANGE_BAND 1.0e-12   /* |x| ≤ 0.7 — must match exactly */
#define TOL_NEAR_BOUNDARY 1.0e-8    /* 0.7 < |x| < 1.0 — loose for precision */

static int run_row(double lam, double x, int n_rev, double truth,
                   int *misses_out)
{
    double k26 = k26astro_lambert_tof_for_test(x, lam, n_rev);
    double abs_err = fabs(k26 - truth);
    double rel_err = (fabs(truth) > 1e-300) ? abs_err / fabs(truth) : abs_err;
    double tol = (fabs(x) <= 0.7) ? TOL_LAGRANGE_BAND : TOL_NEAR_BOUNDARY;
    if (rel_err > tol) {
        printf("MISS  λ=%+.3f x=%+.3f n_rev=%d : "
               "K26=%.17g truth=%.17g abs=%.3e rel=%.3e (tol=%.0e)\n",
               lam, x, n_rev, k26, truth, abs_err, rel_err, tol);
        (*misses_out)++;
        return 0;
    }
    return 1;
}

int main(void)
{
    FILE *f = fopen(CSV_PATH, "rb");
    if (!f) {
        fprintf(stderr,
            "test_lambert_tof_curve: cannot open %s (errno=%d %s)\n"
            "Regenerate via: python3 tests/data/lambert_tof_reference.py "
            "> tests/data/lambert_tof_reference.csv\n",
            CSV_PATH, errno, strerror(errno));
        return 1;
    }

    char line[512];
    int rows = 0, hits = 0, misses = 0;
    while (fgets(line, sizeof(line), f)) {
        /* Skip comment + header lines. */
        if (line[0] == '#' || line[0] == 'l') continue;
        double lam, x, truth;
        int n_rev;
        if (sscanf(line, "%lf,%lf,%d,%lf",
                   &lam, &x, &n_rev, &truth) != 4) {
            continue;
        }
        rows++;
        if (run_row(lam, x, n_rev, truth, &misses)) hits++;
    }
    fclose(f);

    printf("test_lambert_tof_curve: %d rows | %d hits | %d misses\n",
           rows, hits, misses);

    /* Hard gate: tof_curve_ must match the single-branch Lagrange
     * reference across the full grid. Any divergence here flags a
     * TOF-curve regression. */
    if (misses > 0) {
        printf("test_lambert_tof_curve: FAIL: tof_curve_ diverges "
               "from single-branch Lagrange reference at %d rows. "
               "Bisect the misses to localize the regression.\n",
               misses);
        return 1;
    }

    printf("test_lambert_tof_curve: OK (tof_curve_ agrees with "
           "single-branch Lagrange across the test grid)\n");
    return 0;
}
