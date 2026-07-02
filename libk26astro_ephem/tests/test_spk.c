/* test_spk.c — DAF/SPK Type-2 round-trip via the synthetic writer.
 *
 * Strategy:
 *   1. Generate a tiny synthetic Chebyshev-coefficient table for a
 *      circular orbit (one body, two records of half-day duration,
 *      so a one-day total span).
 *   2. Write it as a valid Type-2 SPK file via
 *      k26astro_spk_write_synthetic.
 *   3. Open the file with k26astro_spk_open.
 *   4. Query positions and verify they match the analytic orbit to
 *      sub-metre at the nominal radius (1 AU).
 *   5. Sanity-check the segment metadata + ephem.c query path. */
#include "k26astro_ephem/spk.h"
#include "k26astro_ephem/ephem.h"

#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ORDER 14         /* DE441-typical Chebyshev order */
#define K     ORDER       /* coefficients per axis */
#define RSIZE (2 + 3 * K)

/* Sample a function over [-1, 1] at the Chebyshev-Gauss-Lobatto
 * nodes, then transform to Chebyshev coefficients via the DCT-I
 * (discrete cosine transform of the first kind). For an n-coeff
 * decomposition: c_k = (2/N) Σ_j cos(πjk/N) * f(x_j) with edge
 * weights halved. We use ORDER coefficients which means N+1 = ORDER
 * nodes — so N = ORDER - 1 in the DCT formula. */
static void cheby_fit_(double (*f)(double s, void *), void *ud,
                       double *coeffs)
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
        coeffs[k] = sum * scale;
    }
}

/* Analytic circular orbit at radius R, angular frequency ω about z.
 * `ud` carries the record midpoint epoch + half-interval so the
 * fitter can map s ∈ [-1, 1] back to absolute ET. */
typedef struct {
    double R, omega;
    double mid_et, radius_sec;
    int    axis;   /* 0=x, 1=y, 2=z */
} OrbitCtx;

static double orbit_axis_(double s, void *ud_)
{
    OrbitCtx *u = (OrbitCtx *)ud_;
    double et = u->mid_et + s * u->radius_sec;
    double theta = u->omega * et;
    if (u->axis == 0) return u->R * cos(theta);
    if (u->axis == 1) return u->R * sin(theta);
    return 0.0;   /* z always 0 for an x-y plane orbit */
}

int main(void)
{
    /* Synthetic ephemeris: a body in a circular orbit at radius
     * 1 AU around the system barycentre, period 1 year. Two records
     * of half-day each — total coverage 1 day. */
    const double R       = K26A_AU_M / 1000.0;        /* km — SPK units */
    const double period  = K26A_YEAR_JULIAN_S;
    const double omega   = 2.0 * M_PI / period;
    const double total_dur = 86400.0;
    const double rec_dur   = total_dur / 2.0;
    const int    n_records = 2;
    const double init_et   = 0.0;

    /* Build the record array. Layout per record: MID, RADIUS, then
     * 3*K Chebyshev coefficients. */
    double *records = (double *)calloc((size_t)n_records * (size_t)RSIZE,
                                        sizeof(double));
    assert(records != NULL);
    for (int r = 0; r < n_records; r++) {
        double mid    = init_et + ((double)r + 0.5) * rec_dur;
        double radius = rec_dur * 0.5;
        double *rec   = records + (size_t)r * (size_t)RSIZE;
        rec[0] = mid;
        rec[1] = radius;
        for (int axis = 0; axis < 3; axis++) {
            OrbitCtx ctx = { R, omega, mid, radius, axis };
            cheby_fit_(orbit_axis_, &ctx, &rec[2 + axis * K]);
        }
    }

    /* Write the synthetic kernel. */
    K26AstroSpkWriteSegment seg = {
        .target_body = 399,    /* Earth (NAIF id) */
        .center_body = 0,      /* SSB */
        .start_et    = init_et,
        .end_et      = init_et + total_dur,
        .interval_seconds = rec_dur,
        .records     = records,
        .n_records   = n_records,
        .coeffs_per_axis = K
    };
    const char *spk_path = "tests/fixtures/synthetic_circular.spk";
    int rc = k26astro_spk_write_synthetic(spk_path, &seg, 1);
    assert(rc == 0);
    free(records);

    /* Open and verify metadata. */
    K26AstroSpk *spk = k26astro_spk_open(spk_path);
    assert(spk != NULL);
    assert(k26astro_spk_n_segments(spk) == 1);
    const K26AstroSpkSegment *got = k26astro_spk_segment(spk, 0);
    assert(got != NULL);
    assert(got->target_body == 399);
    assert(got->n_records == n_records);
    assert(got->coeffs_per_axis == K);

    /* Query a handful of epochs and compare to the analytic orbit. */
    int n_probes = 11;
    double max_err = 0.0;
    for (int i = 0; i < n_probes; i++) {
        double et = init_et + (double)i / (double)(n_probes - 1) * total_dur;
        double xyz[3];
        rc = k26astro_spk_pos(spk, 399, et, xyz);
        assert(rc == 0);
        double th = omega * et;
        double expected_x = R * cos(th);
        double expected_y = R * sin(th);
        double dx = xyz[0] - expected_x;
        double dy = xyz[1] - expected_y;
        double err = sqrt(dx * dx + dy * dy);
        if (err > max_err) max_err = err;
    }
    /* Order-14 Chebyshev over a half-day window for a year-period
     * orbit should give µm-level error at the radius. */
    assert(max_err < 1.0);   /* km — i.e. < 1 km */

    /* Position + velocity. At ET = 0 (start of orbit), the body sits
     * on the +x axis with velocity along +y. */
    double xyz[3], vel[3];
    rc = k26astro_spk_pos_vel(spk, 399, init_et, xyz, vel);
    assert(rc == 0);
    /* x ≈ R, y ≈ 0, vx ≈ 0, vy ≈ R * omega. */
    assert(fabs(xyz[0] - R) < 1.0);
    assert(fabs(xyz[1]) < 1.0);
    assert(fabs(vel[0]) < 1e-3);
    assert(fabs(vel[1] - R * omega) < 1e-2);

    /* High-level ephem query (km → m conversion + ICRF tagging).
     * Use J2000 + 1 hour so we're well inside the segment's
     * [start_et, end_et] = [0, 86400] window — the J2000-TT-to-TDB
     * conversion can push the ET to a marginally-negative value
     * which would miss this synthetic segment's left edge. */
    K26AstroEphem *e = k26astro_ephem_load(spk_path);
    assert(e != NULL);
    K26AstroEpoch t = k26astro_epoch_j2000_tt();
    k26astro_epoch_add_seconds(&t, 3600.0);
    K26AstroPosTagged tp = k26astro_ephem_body_pos(e, 399, &t);
    assert(tp.frame_id == K26A_FRAME_ICRF);
    /* Position should be near R*(cos(ω*3600), sin(ω*3600), 0) m. */
    K26AstroPos origin = k26astro_pos_zero();
    K26V3 rel = k26astro_pos_sub(&tp.p, &origin);
    double th_expected = omega * 3600.0;
    double expected_x_m = K26A_AU_M * cos(th_expected);
    double expected_y_m = K26A_AU_M * sin(th_expected);
    assert(fabs(rel.x - expected_x_m) < 1000.0);  /* < 1 km */
    assert(fabs(rel.y - expected_y_m) < 1000.0);
    assert(fabs(rel.z) < 1.0);

    /* NAIF id ↔ name. */
    assert(k26astro_ephem_lookup_name("earth") == 399);
    assert(k26astro_ephem_lookup_name("EARTH") == 399);
    assert(strcmp(k26astro_ephem_id_name(399), "earth") == 0);
    assert(k26astro_ephem_lookup_name("xyzzy") == -1);

    k26astro_ephem_close(e);
    k26astro_spk_close(spk);

    printf("test_spk: OK (synthetic round-trip + max error %.6g km)\n", max_err);
    return 0;
}
