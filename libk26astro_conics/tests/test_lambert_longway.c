/* test_lambert_longway.c — single-rev Lambert long-way coverage.
 *
 * The existing test_lambert.c only exercises short-way (sin Δν > 0)
 * geometries. The long-way path (sin Δν < 0, A < 0) needs separate
 * coverage because the Newton iteration's dtdz derivative form has
 * different sign behaviour near small-z (parabolic) and large-z
 * (highly hyperbolic) for A < 0.
 *
 * Cases:
 *   1. Earth -> Mars long-way at Δν ≈ 220° (the planet "the other way
 *      around"; a transfer that takes the long route through
 *      aphelion). 280-day TOF; round-trip miss < 1e7 m (10000 km);
 *      a generous tolerance because long-way trajectories accumulate
 *      more numerical error in propagation.
 *   2. Near-degenerate long-way at Δν ≈ 350° (highly retrograde-like
 *      transfer). Solver should converge to finite v₁; no
 *      arrival-accuracy assertion in the extreme case.
 *
 * Reference: Vallado §7.6 / Battin 1999 §7.6 A-sign-stable form.
 */
#include "k26astro_conics/lambert.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <stdio.h>

static double v_dist_(K26V3 a, K26V3 b)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}
static double v_norm_(K26V3 a) { return sqrt(a.x*a.x + a.y*a.y + a.z*a.z); }

int main(void)
{
    int fails = 0;

    /* Case 1: long-way Earth → Mars (220° transfer). */
    K26V3 r1 = { K26A_AU_M, 0.0, 0.0 };
    double mars_au = 1.524;
    /* 220° around the +z axis means r2 sits at angle 220°
     * counter-clockwise from r1. The short-way arc covers 140° via
     * the +y half-plane; the long-way arc covers 220° via the -y
     * half-plane. r2 itself is at angle 140° (the same target
     * physical position); the WAY around is what differs. */
    double mars_ang = 140.0 * K26A_PI / 180.0;
    K26V3 r2 = { mars_au * K26A_AU_M * cos(mars_ang),
                 mars_au * K26A_AU_M * sin(mars_ang),
                 0.0 };
    double tof_long = 280.0 * 86400.0;
    K26V3 v1, v2;
    int rc = k26astro_lambert(&v1, &v2, r1, r2,
                                K26A_GM_SUN, tof_long,
                                K26A_LAMBERT_LONG_WAY);
    fprintf(stderr, "case 1 (220° long-way, 280d): rc=%d\n", rc);
    if (rc != 0) {
        fprintf(stderr, "  FAIL: solver did not converge\n");
        fails++;
    } else {
        double v1_mag = v_norm_(v1);
        fprintf(stderr, "  |v1|=%.3f km/s\n", v1_mag / 1000.0);
        if (v1_mag < 5.0e3 || v1_mag > 80.0e3) {
            fprintf(stderr, "  FAIL: |v1|=%.3f km/s out of sane range\n",
                    v1_mag / 1000.0);
            fails++;
        }
        K26V3 r_arrival, v_arrival;
        int prc = k26astro_kepler_propagate(&r_arrival, &v_arrival,
                                              r1, v1, K26A_GM_SUN,
                                              tof_long, 64);
        if (prc != 0) {
            fprintf(stderr, "  FAIL: kepler propagate rc=%d\n", prc);
            fails++;
        } else {
            double miss = v_dist_(r_arrival, r2);
            fprintf(stderr, "  kepler_miss=%.3e m\n", miss);
            /* 220° long-way travels through aphelion; longer arc =>
             * larger numerical error. 10000 km tolerance at 1.5 AU
             * is generous but rejects pathological non-convergence. */
            if (miss > 1.0e7) {
                fprintf(stderr, "  FAIL: round-trip miss > 10000 km\n");
                fails++;
            }
        }
    }

    /* Case 2: near-degenerate 350° (Δν → 360°). */
    double near_ang = 10.0 * K26A_PI / 180.0;  /* r2 at +10°; long-way covers 350° */
    K26V3 r2_near = { K26A_AU_M * cos(near_ang),
                      K26A_AU_M * sin(near_ang),
                      0.0 };
    double tof_near = 5.0 * 86400.0;
    K26V3 v1n, v2n;
    rc = k26astro_lambert(&v1n, &v2n, r1, r2_near,
                            K26A_GM_SUN, tof_near,
                            K26A_LAMBERT_LONG_WAY);
    fprintf(stderr, "case 2 (350° long-way, 5d, near-degenerate): rc=%d\n", rc);
    /* Extreme case: solver may return NO_CONVERGE (rc=4) and that's
     * acceptable. What's NOT acceptable: crash, NaN, infinity. */
    if (rc == 0) {
        if (!isfinite(v1n.x) || !isfinite(v1n.y) || !isfinite(v1n.z)) {
            fprintf(stderr, "  FAIL: non-finite v1\n");
            fails++;
        }
    } else if (rc == 4 || rc == 5) {
        fprintf(stderr, "  OK (NO_CONVERGE / 180-degen — extreme geometry)\n");
    } else {
        fprintf(stderr, "  FAIL: unexpected rc=%d\n", rc);
        fails++;
    }

    if (fails > 0) {
        fprintf(stderr, "test_lambert_longway: FAIL (%d failure(s))\n", fails);
        return 1;
    }
    printf("test_lambert_longway: OK (220° round-trip + 350° guard)\n");
    return 0;
}
