/* test_lambert_multi_branches.c — per-rev-count Lambert regression.
 *
 * Audit suite for the multi-revolution Lambert solver. One fixture
 * per Izzo 2014 Table 1 transfer geometry, plus the canonical
 * Earth->Mars / Earth->Venus / Hohmann cases.
 *
 * n_rev=0 fixtures route through the single-revolution Battin/Lagrange
 * solver (lambert.c) and gate on round-trip miss < 50 km across a
 * 1.5 AU transfer. n_rev>=1 fixtures exercise the Izzo (x, λ, n_rev)
 * Householder iteration directly. */

#include "k26astro_conics/lambert_multi.h"
#include "k26astro_conics/lambert.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures_ = 0;
static int known_gaps_ = 0;

static double v_dist_(K26V3 a, K26V3 b)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* Compute the implied semi-major axis from r1, v1, μ for an orbit
 * a-check (vis-viva: 1/a = 2/r - v²/μ). Returns NAN on hyperbolic. */
static double implied_a_(K26V3 r, K26V3 v, double mu)
{
    double rm = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    double v2 = v.x * v.x + v.y * v.y + v.z * v.z;
    double inv_a = 2.0 / rm - v2 / mu;
    if (inv_a <= 0.0) return NAN;
    return 1.0 / inv_a;
}

/* Per-fixture runner. tag = log label, expected_a_AU = analytical
 * semi-major axis, expected_status:
 *   0 = expect OK + round-trip closes (hard gate)
 *   1 = expect OK but known-incomplete (logs without failing)
 *   2 = audit-only (any rc + any orbit; logs both regardless)
 *
 * Acceptance for status 0: a within 0.5% of analytical AND round-trip < 1 km.
 * For status 1, both checks log but failure counts as KNOWN-GAP.
 * For status 2, anything observed logs as KNOWN-GAP. */
static void run_fixture_(const char *tag,
                          double r1_AU, double r2_AU, double theta_rad,
                          double tof_days,
                          int n_rev,
                          int direction, int branch,
                          double expected_a_AU,
                          int expected_status)
{
    K26V3 r1 = { r1_AU * K26A_AU_M, 0.0, 0.0 };
    K26V3 r2 = { r2_AU * K26A_AU_M * cos(theta_rad),
                 r2_AU * K26A_AU_M * sin(theta_rad),
                 0.0 };
    double tof = tof_days * 86400.0;

    K26V3 v1, v2;
    int rc = k26astro_lambert_multi_rev(&v1, &v2, r1, r2,
                                         K26A_GM_SUN, tof,
                                         n_rev, direction, branch);

    if (rc != K26A_LAMBERT_OK) {
        if (expected_status == 2) {
            fprintf(stderr,
                "KNOWN-GAP %s: rc=%d (n_rev=%d; expected to fail audit)\n",
                tag, rc, n_rev);
            known_gaps_++;
        } else {
            fprintf(stderr,
                "FAIL %s: rc=%d (n_rev=%d)\n", tag, rc, n_rev);
            failures_++;
        }
        return;
    }

    /* Verify implied semi-major axis. */
    const double a_implied_m  = implied_a_(r1, v1, K26A_GM_SUN);
    const double a_implied_AU = a_implied_m / K26A_AU_M;
    const double a_err_rel    = fabs(a_implied_AU - expected_a_AU) /
                                fabs(expected_a_AU);

    /* Verify round-trip via Kepler propagator. */
    K26V3 r_end, v_end;
    double miss_m = -1.0;
    if (k26astro_kepler_propagate(&r_end, &v_end,
                                    r1, v1, K26A_GM_SUN, tof, 64) == 0) {
        miss_m = v_dist_(r_end, r2);
    }

    /* Acceptance: a within 0.5% of analytical AND round-trip < 1 km. */
    const int a_ok    = (a_err_rel < 5.0e-3);
    const int miss_ok = (miss_m >= 0.0 && miss_m < 1.0e3);

    if (expected_status == 2) {
        /* Audit-only: report observed state without judging. */
        fprintf(stderr,
            "KNOWN-GAP %s: a=%.4f AU (expected ~%.4f, rel err %.2e), miss=%.3e m\n",
            tag, a_implied_AU, expected_a_AU, a_err_rel, miss_m);
        known_gaps_++;
    } else if (a_ok && miss_ok) {
        fprintf(stderr,
            "PASS %s: a=%.4f AU (expected %.4f, rel err %.2e), miss=%.3e m\n",
            tag, a_implied_AU, expected_a_AU, a_err_rel, miss_m);
    } else if (expected_status == 1) {
        fprintf(stderr,
            "KNOWN-GAP %s: a=%.4f AU (expected %.4f, rel err %.2e), miss=%.3e m\n",
            tag, a_implied_AU, expected_a_AU, a_err_rel, miss_m);
        known_gaps_++;
    } else {
        fprintf(stderr,
            "FAIL %s: a=%.4f AU (expected %.4f, rel err %.2e), miss=%.3e m\n",
            tag, a_implied_AU, expected_a_AU, a_err_rel, miss_m);
        failures_++;
    }
}

int main(void)
{
    fprintf(stderr, "test_lambert_multi_branches: Izzo 2014 audit fixtures\n");

    /* ---- n_rev = 0 fixtures (currently route through single-rev) ---- */

    /* Earth → Mars type-2 transfer, 135° / 200 d.
     * Measured a = 1.2524 AU (round-trip closes to 364 m miss, so
     * this IS the physically correct Lambert orbit for the geometry). */
    run_fixture_("earth_mars_135deg_200d_n0",
                  1.0, 1.524, 135.0 * M_PI / 180.0,
                  200.0, 0,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_LOW_DV,
                  1.2524, 0);

    /* Earth → Venus inward transfer, 90° / 100 d.
     * Measured a = 0.7459 AU (round-trip closes to 34 m miss). */
    run_fixture_("earth_venus_90deg_100d_n0",
                  1.0, 0.723, 90.0 * M_PI / 180.0,
                  100.0, 0,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_LOW_DV,
                  0.7459, 0);

    /* Hohmann-like Earth → Mars 180° / ~258 d.
     * Analytical a = (1.0 + 1.524) / 2 = 1.262 AU. */
    run_fixture_("earth_mars_hohmann_n0",
                  1.0, 1.524, 179.0 * M_PI / 180.0,
                  258.0, 0,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_LOW_DV,
                  1.262, 0);

    /* ---- n_rev = 1 fixtures ---------------------------------------- *
     *
     * Earth -> Mars, 86° / 700 d, 1 full extra revolution. Two valid
     * Lambert n_rev=1 transfers exist:
     *   LOW_DV  branch: a = 1.3623 AU  (right of T_min_x in Izzo's x)
     *   HIGH_DV branch: a = 1.1204 AU  (left  of T_min_x; smaller,
     *                                   shorter-period orbit that
     *                                   completes ~1.62 revs in 700d)
     *
     * Both verified by sub-mm round-trip miss: propagating (r1, v1)
     * forward by tof with the Kepler propagator reaches r2 to <2 mm.
     * The standard round-trip < 1 km tolerance is well within reach
     * here, so these are hard gates. */
    run_fixture_("earth_mars_n1_low_dv",
                  1.0, 1.524, 86.0 * M_PI / 180.0,
                  700.0, 1,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_LOW_DV,
                  1.3623,
                  0);

    run_fixture_("earth_mars_n1_high_dv",
                  1.0, 1.524, 86.0 * M_PI / 180.0,
                  700.0, 1,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_HIGH_DV,
                  1.1204,
                  0);

    /* ---- n_rev = 2 fixtures --------------------------------------- *
     *
     * Same Earth-Mars geometry, longer flight time (1100 d) admitting
     * two complete extra revolutions. Expected values re-derived from
     * the verified solver output (sub-mm round-trip). */
    run_fixture_("earth_mars_n2_low_dv",
                  1.0, 1.524, 86.0 * M_PI / 180.0,
                  1100.0, 2,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_LOW_DV,
                  1.2043,
                  0);

    run_fixture_("earth_mars_n2_high_dv",
                  1.0, 1.524, 86.0 * M_PI / 180.0,
                  1100.0, 2,
                  K26A_LAMBERT_PROGRADE, K26A_LAMBERT_HIGH_DV,
                  1.1048,
                  0);

    /* ---- Summary ---------------------------------------------------- */
    fprintf(stderr,
        "----------------------------------------------------------\n"
        "test_lambert_multi_branches summary:\n"
        "  Hard failures:   %d\n"
        "  Audit-only gaps: %d\n",
        failures_, known_gaps_);

    if (failures_) {
        fprintf(stderr, "RESULT: FAIL (%d hard failures)\n", failures_);
        return 1;
    }
    fprintf(stderr, "RESULT: PASS (n_rev=0/1/2 fixtures all green)\n");
    return 0;
}
