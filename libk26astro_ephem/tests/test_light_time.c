/* test_light_time.c — light-time-correction fixed-point convergence.
 *
 * Setup:
 *   - Synthetic SPK with the same one-day circular orbit used by
 *     test_spk.c (target body 399, ~1 AU radius)
 *   - Observer fixed at the origin (the heliocentre, with the SPK's
 *     SSB reference)
 *
 * For a target at ~1 AU, light-time is ~500 s. Two-pass iteration
 * should give nanosecond-level convergence. */
#include "k26astro_ephem/ephem.h"
#include "k26astro_ephem/spk.h"

#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ORDER 14
#define K     ORDER
#define RSIZE (2 + 3 * K)

static double R_KM;
static double OMEGA;

typedef struct { double mid, rad; int axis; } OrbitCtx;
static double orbit_axis_(double s, void *ud_)
{
    OrbitCtx *u = (OrbitCtx *)ud_;
    double et = u->mid + s * u->rad;
    double theta = OMEGA * et;
    if (u->axis == 0) return R_KM * cos(theta);
    if (u->axis == 1) return R_KM * sin(theta);
    return 0.0;
}

static void cheby_fit_(double (*f)(double, void *), void *ud,
                       double *out)
{
    int N = ORDER - 1;
    double samples[ORDER];
    for (int j = 0; j <= N; j++) {
        double s = cos(M_PI * (double)j / (double)N);
        samples[j] = f(s, ud);
    }
    for (int k = 0; k < ORDER; k++) {
        double sum = 0.0;
        for (int j = 0; j <= N; j++) {
            double w = (j == 0 || j == N) ? 0.5 : 1.0;
            sum += w * samples[j] * cos(M_PI * (double)j * (double)k / (double)N);
        }
        double scale = (k == 0 || k == N) ? 1.0 / (double)N : 2.0 / (double)N;
        out[k] = sum * scale;
    }
}

int main(void)
{
    R_KM  = K26A_AU_M / 1000.0;
    OMEGA = 2.0 * M_PI / K26A_YEAR_JULIAN_S;

    const double rec_dur = 43200.0;       /* half-day */
    const int    n_recs  = 4;             /* two-day total */
    const double init_et = 0.0;
    double *records = (double *)calloc((size_t)n_recs * (size_t)RSIZE,
                                        sizeof(double));
    assert(records != NULL);
    for (int r = 0; r < n_recs; r++) {
        double mid = init_et + ((double)r + 0.5) * rec_dur;
        double rad = rec_dur * 0.5;
        double *rec = records + (size_t)r * (size_t)RSIZE;
        rec[0] = mid; rec[1] = rad;
        for (int ax = 0; ax < 3; ax++) {
            OrbitCtx ctx = { mid, rad, ax };
            cheby_fit_(orbit_axis_, &ctx, &rec[2 + ax * K]);
        }
    }

    K26AstroSpkWriteSegment seg = {
        .target_body = 399, .center_body = 0,
        .start_et = init_et, .end_et = init_et + n_recs * rec_dur,
        .interval_seconds = rec_dur,
        .records = records, .n_records = n_recs, .coeffs_per_axis = K
    };
    const char *path = "tests/fixtures/lt_orbit.spk";
    int rc = k26astro_spk_write_synthetic(path, &seg, 1);
    assert(rc == 0);
    free(records);

    K26AstroEphem *e = k26astro_ephem_load(path);
    assert(e != NULL);

    /* Observer at heliocentre origin, observation epoch = J2000 + 1h
     * to keep ET clearly inside [start_et, end_et]. */
    K26AstroPos obs = k26astro_pos_zero();
    K26AstroEpoch t_obs = k26astro_epoch_j2000_tt();
    k26astro_epoch_add_seconds(&t_obs, 3600.0);
    k26astro_epoch_convert(&t_obs, K26A_TS_TDB);

    /* Geometric baseline — body position at t_obs, no light-time
     * correction. Using ephem_body_pos directly rather than
     * ephem_observe(.., 0) makes the comparison's two endpoints
     * explicit: a pure ephemeris query vs. a light-time-corrected
     * observation. */
    K26AstroPosTagged geom = k26astro_ephem_body_pos(e, 399, &t_obs);
    assert(geom.frame_id == K26A_FRAME_ICRF);

    /* Apparent (max_iter = 4) — should converge cleanly. */
    K26AstroPosTagged appar = k26astro_ephem_observe(e, 399, &obs, &t_obs, 4);
    assert(appar.frame_id == K26A_FRAME_ICRF);

    /* The apparent position should differ from the geometric one by
     * the body's motion during the light-time τ ≈ R/c ≈ 500 s.
     * For a 1-year-period circular orbit at 1 AU, the body moves
     * R * ω * τ ≈ R * (2π/yr) * 500s ≈ 1.5e8 m ≈ 50 km/s × 500 s ≈
     * 25 000 km. */
    K26V3 d = k26astro_pos_sub(&geom.p, &appar.p);
    double mag = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    /* Tolerance loose — exact analytic motion is omega * R * τ,
     * τ from R/c. */
    double tau = K26A_AU_M / K26A_C;
    double expected_motion = K26A_AU_M * OMEGA * tau;
    assert(fabs(mag - expected_motion) < 1.0e4);  /* within 10 km */

    /* Re-running with more iterations should not move the result
     * by more than a metre (the v0.1 convergence is sub-nanosecond
     * for inner solar system). */
    K26AstroPosTagged appar2 = k26astro_ephem_observe(e, 399, &obs, &t_obs, 8);
    K26V3 dd = k26astro_pos_sub(&appar.p, &appar2.p);
    double resid = sqrt(dd.x * dd.x + dd.y * dd.y + dd.z * dd.z);
    assert(resid < 1.0);   /* < 1 metre */

    k26astro_ephem_close(e);

    printf("test_light_time: OK (light-time motion %.3f km, "
           "iteration residual %.6g m)\n",
           mag / 1000.0, resid);
    return 0;
}
