/* test_chain_walker.c - SPK observer-chain walker.
 *
 * DE-series SPK kernels store each body relative to its centre, not
 * relative to the Solar System Barycentre. DE441's Earth segment
 * (NAIF 399) carries center_body=3 (Earth-Moon barycentre), so a
 * naive reader that ignores center_body returns Earth-relative-to-
 * EMB and treats it as Earth-relative-to-SSB; resolving the chain
 * back to the requested observer is required to avoid large
 * residuals.
 *
 * This test builds a synthetic 3-hop chain (target=99 → 2 → 1 → 0
 * SSB) with constant-offset segments and asserts the resolved
 * position equals the vector sum. Constants make the test
 * tolerance-free (no Chebyshev approximation error to budget).
 */
#include "k26astro_ephem/spk.h"
#include "k26astro_ephem/ephem.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Coeffs-per-axis = 2 (c[0] = offset, c[1] = 0 → constant). */
#define K     2
#define RSIZE (2 + 3 * K)

/* Build one constant-offset record at midpoint `mid` with half-
 * interval `radius`. SPK record layout: [MID, RADIUS, cx[K],
 * cy[K], cz[K]]. */
static void
build_constant_record_(double *out, double mid, double radius,
                       double off_x_km, double off_y_km, double off_z_km)
{
    out[0]                 = mid;
    out[1]                 = radius;
    out[2 + 0 * K + 0]     = off_x_km;  out[2 + 0 * K + 1] = 0.0;
    out[2 + 1 * K + 0]     = off_y_km;  out[2 + 1 * K + 1] = 0.0;
    out[2 + 2 * K + 0]     = off_z_km;  out[2 + 2 * K + 1] = 0.0;
}

int main(void)
{
    /* Three-hop chain. Total span = 1 day, single record per
     * segment, midpoint at 43200 s = noon-after-J2000. */
    const double mid    = 43200.0;
    const double radius = 43200.0;
    const double start_et = 0.0;
    const double end_et   = 86400.0;

    /* Hop offsets (km). Chosen so the resolved SSB-relative for
     * target=99 is non-trivial in all three axes and the
     * individual hops don't cancel. */
    const double off99_km[3]  = { 1.0e6,  0.0,    0.0   };
    const double off2_km[3]   = { 0.0,    2.0e6,  0.0   };
    const double off1_km[3]   = { 0.0,    0.0,    3.0e6 };

    /* Build per-segment record arrays. */
    double rec99[RSIZE], rec2[RSIZE], rec1[RSIZE];
    build_constant_record_(rec99, mid, radius,
                           off99_km[0], off99_km[1], off99_km[2]);
    build_constant_record_(rec2,  mid, radius,
                           off2_km[0],  off2_km[1],  off2_km[2]);
    build_constant_record_(rec1,  mid, radius,
                           off1_km[0],  off1_km[1],  off1_km[2]);

    K26AstroSpkWriteSegment segs[3] = {
        { .target_body = 99, .center_body = 2,
          .start_et = start_et, .end_et = end_et,
          .interval_seconds = 86400.0,
          .records = rec99, .n_records = 1, .coeffs_per_axis = K },
        { .target_body = 2,  .center_body = 1,
          .start_et = start_et, .end_et = end_et,
          .interval_seconds = 86400.0,
          .records = rec2,  .n_records = 1, .coeffs_per_axis = K },
        { .target_body = 1,  .center_body = 0,
          .start_et = start_et, .end_et = end_et,
          .interval_seconds = 86400.0,
          .records = rec1,  .n_records = 1, .coeffs_per_axis = K },
    };

    const char *spk_path = "tests/fixtures/synthetic_chain.spk";
    int rc = k26astro_spk_write_synthetic(spk_path, segs, 3);
    assert(rc == 0);

    K26AstroEphem *e = k26astro_ephem_load(spk_path);
    assert(e != NULL);

    /* Query target=99 at the segment midpoint. Resolved SSB-
     * relative position should be the vector sum of all three
     * hop offsets, in metres. */
    K26AstroEpoch t = k26astro_epoch_j2000_tt();
    /* The segments cover ET in [0, 86400] TDB. J2000_TT-to-TDB
     * lands close to 0; nudge into the segment by adding 1 hr. */
    k26astro_epoch_add_seconds(&t, 3600.0);
    K26AstroPosTagged tp = k26astro_ephem_body_pos(e, 99, &t);
    assert(tp.frame_id == K26A_FRAME_ICRF);

    K26AstroPos origin = k26astro_pos_zero();
    K26V3 rel = k26astro_pos_sub(&tp.p, &origin);

    const double expected_x_m = (off99_km[0] + off2_km[0] + off1_km[0]) * 1000.0;
    const double expected_y_m = (off99_km[1] + off2_km[1] + off1_km[1]) * 1000.0;
    const double expected_z_m = (off99_km[2] + off2_km[2] + off1_km[2]) * 1000.0;

    /* Constants are exactly representable through Chebyshev with
     * c[0]=offset, c[1]=0; the chain sum is three exact additions.
     * The only error path is the s_norm → cheby_eval, which for a
     * constant series returns c[0] regardless of s. So the
     * tolerance can be tight. */
    assert(fabs(rel.x - expected_x_m) < 1.0e-3);  /* < 1 mm */
    assert(fabs(rel.y - expected_y_m) < 1.0e-3);
    assert(fabs(rel.z - expected_z_m) < 1.0e-3);

    /* Verify single-hop case is unchanged: query target=1 directly.
     * Its centre is SSB; only one hop. Expected = off1_km in m. */
    K26AstroPosTagged tp1 = k26astro_ephem_body_pos(e, 1, &t);
    assert(tp1.frame_id == K26A_FRAME_ICRF);
    K26V3 rel1 = k26astro_pos_sub(&tp1.p, &origin);
    assert(fabs(rel1.x - off1_km[0] * 1000.0) < 1.0e-3);
    assert(fabs(rel1.y - off1_km[1] * 1000.0) < 1.0e-3);
    assert(fabs(rel1.z - off1_km[2] * 1000.0) < 1.0e-3);

    /* Mid-chain case: query target=2. Chain is 2 → 1 → 0. Sum is
     * off2 + off1. */
    K26AstroPosTagged tp2 = k26astro_ephem_body_pos(e, 2, &t);
    assert(tp2.frame_id == K26A_FRAME_ICRF);
    K26V3 rel2 = k26astro_pos_sub(&tp2.p, &origin);
    const double exp2_x = (off2_km[0] + off1_km[0]) * 1000.0;
    const double exp2_y = (off2_km[1] + off1_km[1]) * 1000.0;
    const double exp2_z = (off2_km[2] + off1_km[2]) * 1000.0;
    assert(fabs(rel2.x - exp2_x) < 1.0e-3);
    assert(fabs(rel2.y - exp2_y) < 1.0e-3);
    assert(fabs(rel2.z - exp2_z) < 1.0e-3);

    /* Missing body: query an id with no segment. Should return
     * INVALID frame cleanly (no crash). */
    K26AstroPosTagged tm = k26astro_ephem_body_pos(e, 12345, &t);
    assert(tm.frame_id == K26A_FRAME_INVALID);

    k26astro_ephem_close(e);

    printf("test_chain_walker: OK (3-hop SSB resolution exact "
           "to < 1 mm)\n");
    return 0;
}
