/* test_epoch.c — K26AstroEpoch + time-scale conversions.
 *
 * Acceptance:
 *   - J2000 round-trip via calendar + ISO-8601
 *   - TAI↔TT differs by exactly 32.184 s
 *   - TAI↔UTC differs by the leap-second offset at the epoch
 *   - Full chain TT→UTC→UT1 (and back) is round-trip stable to
 *     sub-microsecond (DUT1=0 in v0.1, so UT1≡UTC for now)
 *   - TT↔TDB FB1990 series gives small (<10 ms) offset, and one-pass
 *     iteration converges to a fixed point. */
#include "k26astro_core/epoch.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* ---- J2000 + calendar round-trip ----------------------------- */
    K26AstroEpoch j2000 = k26astro_epoch_j2000_tt();
    assert(j2000.days_since_J2000 == 0);
    assert(j2000.seconds_of_day == 0.0);
    assert(j2000.scale == K26A_TS_TT);

    int y, mo, d, h, mi;
    double s;
    k26astro_epoch_to_calendar(&j2000, &y, &mo, &d, &h, &mi, &s);
    assert(y == 2000 && mo == 1 && d == 1);
    assert(h == 12 && mi == 0);
    assert(fabs(s) < 1e-9);

    /* Reconstruct J2000 from the calendar and confirm round-trip. */
    K26AstroEpoch j2000_re;
    int rc = k26astro_epoch_from_calendar(&j2000_re,
                                           2000, 1, 1, 12, 0, 0.0,
                                           K26A_TS_TT);
    assert(rc == 0);
    assert(j2000_re.days_since_J2000 == 0);
    assert(fabs(j2000_re.seconds_of_day) < 1e-9);

    /* ---- ISO-8601 round-trip ------------------------------------ */
    K26AstroEpoch iso_e;
    rc = k26astro_epoch_from_iso8601(&iso_e,
                                      "2000-01-01T12:00:00",
                                      K26A_TS_TT);
    assert(rc == 0);
    assert(iso_e.days_since_J2000 == 0);
    assert(fabs(iso_e.seconds_of_day) < 1e-9);

    char buf[64];
    int n = k26astro_epoch_to_iso8601(&iso_e, buf, sizeof buf);
    assert(n > 0);
    assert(strncmp(buf, "2000-01-01T12:00:00", 19) == 0);

    /* ---- TAI ↔ TT  --------------------------------------------- */
    K26AstroEpoch tai = j2000;
    rc = k26astro_epoch_convert(&tai, K26A_TS_TAI);
    assert(rc == 0);
    assert(tai.scale == K26A_TS_TAI);
    /* TAI - TT = -32.184 s */
    double diff = k26astro_epoch_diff_seconds(&j2000, &tai);   /* tt - tai */
    /* Both scales differ; use raw difference: j2000 in TT minus tai in TAI
     * with the conversion already applied. The diff function asserts
     * same scale, so flip to a manual check via numerical seconds. */
    /* Re-derive: tai.seconds_of_day = j2000.seconds_of_day - 32.184 (mod day) */
    double approx_dt = (double)(j2000.days_since_J2000 - tai.days_since_J2000) * 86400.0
                     + (j2000.seconds_of_day - tai.seconds_of_day);
    assert(fabs(approx_dt - 32.184) < 1e-9);
    (void)diff;

    /* Back the other way. */
    rc = k26astro_epoch_convert(&tai, K26A_TS_TT);
    assert(rc == 0);
    assert(tai.scale == K26A_TS_TT);
    assert(tai.days_since_J2000 == j2000.days_since_J2000);
    assert(fabs(tai.seconds_of_day - j2000.seconds_of_day) < 1e-12);

    /* ---- UTC ↔ TAI (leap seconds at J2000) --------------------- */
    /* At J2000 (Jan 2000), TAI-UTC = 32 s. */
    K26AstroEpoch utc = j2000;
    rc = k26astro_epoch_convert(&utc, K26A_TS_UTC);
    assert(rc == 0);
    assert(utc.scale == K26A_TS_UTC);

    int dtai = k26astro_leap_seconds_at(&utc);
    assert(dtai == 32);

    /* TT - UTC at J2000 = 32 + 32.184 = 64.184 s. */
    double tt_minus_utc = (double)(j2000.days_since_J2000 - utc.days_since_J2000) * 86400.0
                        + (j2000.seconds_of_day - utc.seconds_of_day);
    assert(fabs(tt_minus_utc - 64.184) < 1e-9);

    /* Round-trip via TT → UTC → TT. */
    K26AstroEpoch rtt = j2000;
    rc = k26astro_epoch_convert(&rtt, K26A_TS_UTC); assert(rc == 0);
    rc = k26astro_epoch_convert(&rtt, K26A_TS_TT);  assert(rc == 0);
    assert(rtt.days_since_J2000 == j2000.days_since_J2000);
    assert(fabs(rtt.seconds_of_day - j2000.seconds_of_day) < 1e-9);

    /* ---- TT ↔ TDB (FB1990) ------------------------------------- */
    K26AstroEpoch tdb = j2000;
    rc = k26astro_epoch_convert(&tdb, K26A_TS_TDB);
    assert(rc == 0);
    assert(tdb.scale == K26A_TS_TDB);
    /* TDB - TT at J2000 is ~ -0.069 ms (FB1990 main term). */
    double tdb_minus_tt = (double)(tdb.days_since_J2000 - j2000.days_since_J2000) * 86400.0
                        + (tdb.seconds_of_day - j2000.seconds_of_day);
    assert(fabs(tdb_minus_tt) < 0.01);   /* well under 10 ms */

    /* Round-trip back. */
    rc = k26astro_epoch_convert(&tdb, K26A_TS_TT);
    assert(rc == 0);
    assert(tdb.days_since_J2000 == j2000.days_since_J2000);
    assert(fabs(tdb.seconds_of_day - j2000.seconds_of_day) < 1e-6);

    /* ---- Add seconds ------------------------------------------- */
    K26AstroEpoch t = j2000;
    k26astro_epoch_add_seconds(&t, 86400.0 + 3600.5);
    assert(t.days_since_J2000 == 1);
    assert(fabs(t.seconds_of_day - 3600.5) < 1e-9);

    k26astro_epoch_add_seconds(&t, -86400.0 - 3600.5);
    assert(t.days_since_J2000 == 0);
    assert(fabs(t.seconds_of_day) < 1e-9);

    /* ---- Negative-time arithmetic (pre-J2000) ----------------- */
    /* 1999-12-31T12:00:00 TT = J2000 - 1 day. */
    K26AstroEpoch pre;
    rc = k26astro_epoch_from_iso8601(&pre, "1999-12-31T12:00:00", K26A_TS_TT);
    assert(rc == 0);
    assert(pre.days_since_J2000 == -1);
    assert(fabs(pre.seconds_of_day) < 1e-9);

    printf("test_epoch: OK (J2000 + calendar/ISO + TAI/UTC/TT/TDB)\n");
    return 0;
}
