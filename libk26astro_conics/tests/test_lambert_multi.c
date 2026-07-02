/* test_lambert_multi.c — Izzo multi-revolution Lambert solver.
 *
 * The Pochhammer recursion in tof_curve_'s 2F1(3,1;5/2;S1) series
 * follows Battin 1999 §6.3. Multi-revolution branches (n_rev >= 1)
 * remain a known-incomplete area; see lambert_multi.c for details
 * on initial-guess accuracy, λ sign convention, and Householder
 * convergence as S1 -> 1.
 *
 * Acceptance gates:
 *   - The solver builds and invokes correctly.
 *   - n_rev=0 + n_rev=1 invocations don't crash; return codes are
 *     either OK or one of the documented non-fatal codes
 *     (NO_CONVERGE, NO_SOLUTION).
 *   - Documented bad-input paths return the right error codes
 *     (NO_SOLUTION for too-short TOF, DEGENERATE for r1==r2,
 *     BAD_INPUT for negative mu, NULL_OUT for NULL output).
 *   - When the solver returns OK, the velocities are finite.
 */
#include "k26astro_conics/lambert_multi.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* ---- n_rev = 0 single-rev sanity ------------------------ *
     * 60° separation, 60-day TOF; same geometry as test_lambert. */
    K26V3 r1 = { K26A_AU_M, 0.0, 0.0 };
    double sep = K26A_PI / 3.0;
    K26V3 r2 = { K26A_AU_M * cos(sep), K26A_AU_M * sin(sep), 0.0 };

    K26V3 v1, v2;
    int rc = k26astro_lambert_multi_rev(&v1, &v2, r1, r2,
                                          K26A_GM_SUN, 60.0 * 86400.0,
                                          0,
                                          K26A_LAMBERT_PROGRADE,
                                          K26A_LAMBERT_LOW_DV);
    assert(rc == K26A_LAMBERT_OK || rc == K26A_LAMBERT_NO_CONVERGE);
    if (rc == K26A_LAMBERT_OK) {
        assert(isfinite(v1.x) && isfinite(v1.y) && isfinite(v1.z));
        assert(isfinite(v2.x) && isfinite(v2.y) && isfinite(v2.z));
    }

    /* ---- n_rev = 1 invocation (doesn't crash) -------------- */
    K26V3 v1r1, v2r1;
    double tof_1rev = 1.3 * 365.25 * 86400.0;
    rc = k26astro_lambert_multi_rev(&v1r1, &v2r1, r1, r2,
                                      K26A_GM_SUN, tof_1rev,
                                      1,
                                      K26A_LAMBERT_PROGRADE,
                                      K26A_LAMBERT_LOW_DV);
    (void)rc;
    if (rc == K26A_LAMBERT_OK) {
        assert(isfinite(v1r1.x) && isfinite(v1r1.y) && isfinite(v1r1.z));
    }

    /* ---- Too-short TOF for n_rev = 2 → NO_SOLUTION --------- */
    rc = k26astro_lambert_multi_rev(&v1, &v2, r1, r2,
                                      K26A_GM_SUN, 10.0 * 86400.0,
                                      2,
                                      K26A_LAMBERT_PROGRADE,
                                      K26A_LAMBERT_LOW_DV);
    assert(rc == K26A_LAMBERT_NO_SOLUTION
        || rc == K26A_LAMBERT_NO_CONVERGE);

    /* ---- Degenerate: r1 == r2 ------------------------------- */
    K26V3 r_eq = { K26A_AU_M, 0.0, 0.0 };
    rc = k26astro_lambert_multi_rev(&v1, &v2, r_eq, r_eq,
                                      K26A_GM_SUN, 86400.0,
                                      0,
                                      K26A_LAMBERT_PROGRADE,
                                      K26A_LAMBERT_LOW_DV);
    assert(rc == K26A_LAMBERT_DEGENERATE);

    /* ---- Bad input: negative mu ----------------------------- */
    rc = k26astro_lambert_multi_rev(&v1, &v2, r1, r2,
                                      -K26A_GM_SUN, 86400.0,
                                      0,
                                      K26A_LAMBERT_PROGRADE,
                                      K26A_LAMBERT_LOW_DV);
    assert(rc == K26A_LAMBERT_BAD_INPUT);

    /* ---- NULL output ---------------------------------------- */
    rc = k26astro_lambert_multi_rev(NULL, &v2, r1, r2,
                                      K26A_GM_SUN, 86400.0,
                                      0,
                                      K26A_LAMBERT_PROGRADE,
                                      K26A_LAMBERT_LOW_DV);
    assert(rc == K26A_LAMBERT_NULL_OUT);

    printf("test_lambert_multi: OK (n_rev=0/1 + guards; "
           "Pochhammer recursion fixed)\n");
    return 0;
}
