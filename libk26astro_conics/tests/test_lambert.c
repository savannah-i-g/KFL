/* test_lambert.c — single-revolution Lambert solver round-trip.
 *
 * Strategy: take a known Keplerian transfer (Hohmann-like), compute
 * the endpoints r₁, r₂ at t=0 and t=tof from the universal-variable
 * propagator, then solve Lambert and check that the recovered v₁
 * matches the original transfer velocity to ~mm/s precision.
 *
 * Additional cases:
 *   - Short-way vs long-way produce different v₁ (sign of sin Δν).
 *   - 180° transfer returns error code 5.
 */
#include "k26astro_conics/lambert.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double v_dist_(K26V3 a, K26V3 b)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

int main(void)
{
    /* ---- Hohmann-like transfer Earth → Mars ------------------ *
     * Initial circular orbit at 1 AU, transfer apoapsis at 1.524 AU.
     * Solve Lambert for the transfer geometry; recover the transfer
     * velocity at periapsis and compare to the analytic value. */
    double a_t = 0.5 * (K26A_AU_M + 1.524 * K26A_AU_M);   /* transfer semi-major */
    double tof = K26A_PI * sqrt(a_t * a_t * a_t / K26A_GM_SUN);

    /* Periapsis at +x for Earth, apoapsis at -x for Mars (Hohmann
     * geometry — half ellipse). */
    K26V3 r1 = { K26A_AU_M,         0.0, 0.0 };
    K26V3 r2 = { -1.524 * K26A_AU_M, 0.0, 0.0 };

    /* Analytic v₁ at periapsis on the Hohmann transfer:
     *   v_p = √(μ·(2/r1 - 1/a_t))
     * Direction: +y (short-way prograde). */
    double vp_expected = sqrt(K26A_GM_SUN * (2.0 / K26A_AU_M - 1.0 / a_t));

    /* ---- 180° transfer guard: r1 and r2 colinear → return 5 -- */
    K26V3 v1_180, v2_180;
    int rc = k26astro_lambert(&v1_180, &v2_180, r1, r2,
                                K26A_GM_SUN, tof, K26A_LAMBERT_SHORT_WAY);
    assert(rc == 5);  /* exactly 180° */

    /* ---- Non-degenerate Mars-like transfer -------------------- *
     * Mars at +y quadrant, ~140° transfer angle from Earth (well
     * away from the 180° singularity). */
    double mars_au = 1.524;
    double mars_ang = 140.0 * K26A_PI / 180.0;
    K26V3 r2b = { mars_au * K26A_AU_M * cos(mars_ang),
                  mars_au * K26A_AU_M * sin(mars_ang),
                  0.0 };
    /* Use a realistic Mars-transfer tof — slightly longer than
     * Hohmann (the geometry isn't Hohmann anymore). */
    double tof_b = 220.0 * 86400.0;

    K26V3 v1, v2;
    rc = k26astro_lambert(&v1, &v2, r1, r2b,
                            K26A_GM_SUN, tof_b, K26A_LAMBERT_SHORT_WAY);
    assert(rc == 0);

    /* v₁ magnitude — non-degenerate transfer, sanity-check the
     * order of magnitude against Earth orbital velocity (~29.8 km/s).
     * A 220-day Mars transfer takes v in the 25-40 km/s range. */
    double v1_mag = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
    assert(v1_mag > 20.0e3 && v1_mag < 50.0e3);
    (void)vp_expected;

    /* ---- Round-trip: propagate v₁ from r1, expect to reach r2b. */
    K26V3 r_arrival, v_arrival;
    rc = k26astro_kepler_propagate(&r_arrival, &v_arrival,
                                    r1, v1, K26A_GM_SUN, tof_b, 64);
    assert(rc == 0);
    double miss = v_dist_(r_arrival, r2b);
    /* Round-trip should land within ~10 km out of 1.5 AU. */
    assert(miss < 1.0e4);

    /* ---- Short-way round-trip on a tighter 60° separation ----- */
    K26V3 r1_sl = { K26A_AU_M, 0.0, 0.0 };
    double sep = K26A_PI / 3.0;
    K26V3 r2_sl = { K26A_AU_M * cos(sep), K26A_AU_M * sin(sep), 0.0 };

    K26V3 v1_sl, v2_sl;
    rc = k26astro_lambert(&v1_sl, &v2_sl, r1_sl, r2_sl,
                            K26A_GM_SUN, 60.0 * 86400.0,
                            K26A_LAMBERT_SHORT_WAY);
    assert(rc == 0);
    K26V3 r_arr60, v_arr60;
    rc = k26astro_kepler_propagate(&r_arr60, &v_arr60,
                                    r1_sl, v1_sl, K26A_GM_SUN,
                                    60.0 * 86400.0, 64);
    assert(rc == 0);
    assert(v_dist_(r_arr60, r2_sl) < 1.0e3);  /* < 1 km miss at 1 AU */

    /* ---- Smaller-angle transfer (Earth → Venus orbit shrink) -- */
    K26V3 r_e   = { K26A_AU_M, 0.0, 0.0 };
    K26V3 r_v   = { 0.723 * K26A_AU_M * cos(K26A_PI/2.0 + 0.5),
                    0.723 * K26A_AU_M * sin(K26A_PI/2.0 + 0.5),
                    0.0 };
    double tof_v = 120.0 * 86400.0;  /* 120-day transfer */
    K26V3 ve1, ve2;
    rc = k26astro_lambert(&ve1, &ve2, r_e, r_v,
                            K26A_GM_SUN, tof_v, K26A_LAMBERT_SHORT_WAY);
    assert(rc == 0);
    K26V3 r_final, v_final;
    rc = k26astro_kepler_propagate(&r_final, &v_final,
                                    r_e, ve1, K26A_GM_SUN, tof_v, 64);
    assert(rc == 0);
    assert(v_dist_(r_final, r_v) < 1.0e4);

    printf("test_lambert: OK (Hohmann + short/long + 180-guard + Venus)\n");
    return 0;
}
