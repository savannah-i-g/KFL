/* epoch.c — K26AstroEpoch + time-scale arithmetic.
 *
 * Day-boundary convention: noon TT (Julian-Date aligned). See
 * include/k26astro_core/epoch.h for the rationale.
 *
 * Leap-second table: hardcoded for the full IERS history (1972 →
 * 2026). The post-2035 phaseout (28th CGPM Resolution) doesn't
 * remove the historical record; it just means no new entries.
 *
 * TT-TDB: FB1990 single-term approximation, accurate to ~10 μs over
 * 1900-2100. A TE405-class sub-nanosecond backend is not currently
 * provided. */
#include "k26astro_core/epoch.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* ---- Julian Date arithmetic ------------------------------------ */

/* J2000 = JD 2451545.0 = 2000-01-01T12:00:00 TT. The day before
 * J2000 noon (i.e. 2000-01-01T00:00:00 TT) has JD 2451544.5. */
#define K26A_JD_J2000 2451545.0

/* Gregorian → JD via Fliegel-Van Flandern (1968). Valid for any
 * proleptic Gregorian date that fits in int64 (which is everywhere
 * we care about — pre-1582 history is treated as if Gregorian, the
 * "proleptic" calendar). */
static double calendar_to_jd_(int year, int month, int day,
                              int hour, int min, double sec)
{
    int64_t y = (int64_t)year;
    int64_t m = (int64_t)month;
    int64_t d = (int64_t)day;
    /* FvF: jdn at noon. */
    int64_t a = (14 - m) / 12;
    int64_t y_ = y + 4800 - a;
    int64_t m_ = m + 12 * a - 3;
    int64_t jdn = d + (153 * m_ + 2) / 5 + 365 * y_
                + y_ / 4 - y_ / 100 + y_ / 400 - 32045;
    double  jd_at_noon = (double)jdn;
    double  frac_of_day_from_noon =
        ((double)hour - 12.0) / 24.0
      +  (double)min  / 1440.0
      +          sec   / 86400.0;
    return jd_at_noon + frac_of_day_from_noon;
}

/* JD → Gregorian (also Fliegel-Van Flandern, inverse). */
static void jd_to_calendar_(double jd,
                            int *out_y, int *out_m, int *out_d,
                            int *out_hr, int *out_mn, double *out_sec)
{
    /* Day fraction relative to JD noon. */
    int64_t jdn = (int64_t)floor(jd + 0.5);
    double  day_frac = (jd + 0.5) - (double)jdn;   /* in [0, 1) */

    int64_t a = jdn + 32044;
    int64_t b = (4 * a + 3) / 146097;
    int64_t c = a - (146097 * b) / 4;
    int64_t d = (4 * c + 3) / 1461;
    int64_t e = c - (1461 * d) / 4;
    int64_t m = (5 * e + 2) / 153;

    int day   = (int)(e - (153 * m + 2) / 5 + 1);
    int month = (int)(m + 3 - 12 * (m / 10));
    int year  = (int)(100 * b + d - 4800 + m / 10);

    if (out_y) *out_y = year;
    if (out_m) *out_m = month;
    if (out_d) *out_d = day;

    double secs = day_frac * 86400.0;
    int hh = (int)floor(secs / 3600.0);
    secs -= (double)hh * 3600.0;
    int mm = (int)floor(secs / 60.0);
    secs -= (double)mm * 60.0;
    if (out_hr)  *out_hr  = hh;
    if (out_mn)  *out_mn  = mm;
    if (out_sec) *out_sec = secs;
}

/* ---- Epoch <-> JD --------------------------------------------- */

K26AstroEpoch k26astro_epoch_j2000_tt(void)
{
    K26AstroEpoch e;
    e.days_since_J2000 = 0;
    e.seconds_of_day   = 0.0;
    e.scale            = K26A_TS_TT;
    memset(e._pad, 0, sizeof(e._pad));
    return e;
}

double k26astro_epoch_to_jd(const K26AstroEpoch *e)
{
    if (!e) return 0.0;
    return K26A_JD_J2000
         + (double)e->days_since_J2000
         + e->seconds_of_day / 86400.0;
}

K26AstroEpoch k26astro_epoch_from_jd(double jd, K26AstroTimeScale scale)
{
    K26AstroEpoch e;
    double offset_days = jd - K26A_JD_J2000;
    int64_t whole = (int64_t)floor(offset_days);
    double  frac  = offset_days - (double)whole;   /* [0, 1) */
    e.days_since_J2000 = whole;
    e.seconds_of_day   = frac * 86400.0;
    /* Defensive normalisation against rounding in the
     * floor/subtract pair landing exactly at 86400. */
    while (e.seconds_of_day >= 86400.0) {
        e.seconds_of_day -= 86400.0;
        e.days_since_J2000++;
    }
    while (e.seconds_of_day < 0.0) {
        e.seconds_of_day += 86400.0;
        e.days_since_J2000--;
    }
    e.scale = (uint8_t)scale;
    memset(e._pad, 0, sizeof(e._pad));
    return e;
}

/* ---- Calendar ------------------------------------------------- */

int k26astro_epoch_from_calendar(K26AstroEpoch *out,
                                  int year, int month, int day,
                                  int hour, int min, double sec,
                                  K26AstroTimeScale scale)
{
    if (!out) return 1;
    if (month < 1 || month > 12) return 2;
    if (day   < 1 || day   > 31) return 2;
    if (hour  < 0 || hour  > 23) return 2;
    if (min   < 0 || min   > 59) return 2;
    if (sec   < 0.0 || sec >= 60.0) return 2;
    double jd = calendar_to_jd_(year, month, day, hour, min, sec);
    *out = k26astro_epoch_from_jd(jd, scale);
    return 0;
}

void k26astro_epoch_to_calendar(const K26AstroEpoch *e,
                                 int *out_year, int *out_month, int *out_day,
                                 int *out_hour, int *out_min, double *out_sec)
{
    if (!e) return;
    jd_to_calendar_(k26astro_epoch_to_jd(e),
                    out_year, out_month, out_day,
                    out_hour, out_min, out_sec);
}

/* ---- ISO-8601 ------------------------------------------------- */

int k26astro_epoch_from_iso8601(K26AstroEpoch *out,
                                 const char *iso,
                                 K26AstroTimeScale scale)
{
    if (!out || !iso) return 1;
    int  y, mo, d, h, mi;
    double s = 0.0;
    /* Accept `YYYY-MM-DDTHH:MM:SS` and `YYYY-MM-DDTHH:MM:SS.ffffff`.
     * No time-zone suffix — scale supplied via the `scale` arg. */
    int n = sscanf(iso, "%4d-%2d-%2dT%2d:%2d:%lf", &y, &mo, &d, &h, &mi, &s);
    if (n != 6) {
        /* Allow date-only (assumed midnight). */
        n = sscanf(iso, "%4d-%2d-%2d", &y, &mo, &d);
        if (n != 3) return 2;
        h = 0; mi = 0; s = 0.0;
    }
    return k26astro_epoch_from_calendar(out, y, mo, d, h, mi, s, scale);
}

int k26astro_epoch_to_iso8601(const K26AstroEpoch *e,
                               char *buf, size_t bufsz)
{
    if (!e || !buf || bufsz < 27) return -1;
    int y, mo, d, h, mi;
    double s;
    k26astro_epoch_to_calendar(e, &y, &mo, &d, &h, &mi, &s);
    int n = snprintf(buf, bufsz, "%04d-%02d-%02dT%02d:%02d:%09.6f",
                     y, mo, d, h, mi, s);
    return n;
}

/* ---- Difference / addition ----------------------------------- */

double k26astro_epoch_diff_seconds(const K26AstroEpoch *a,
                                    const K26AstroEpoch *b)
{
    if (!a || !b) return 0.0;
    double dd = (double)(a->days_since_J2000 - b->days_since_J2000);
    return dd * 86400.0 + (a->seconds_of_day - b->seconds_of_day);
}

void k26astro_epoch_add_seconds(K26AstroEpoch *e, double delta_s)
{
    if (!e) return;
    e->seconds_of_day += delta_s;
    /* Normalise into [0, 86400). For large |delta_s|, fold via
     * truncation; the int64 day count gives plenty of headroom. */
    if (e->seconds_of_day >= 86400.0 || e->seconds_of_day < 0.0) {
        double whole_days = floor(e->seconds_of_day / 86400.0);
        e->seconds_of_day -= whole_days * 86400.0;
        e->days_since_J2000 += (int64_t)whole_days;
    }
    /* Final safety pass for boundary cases (seconds_of_day landing
     * at exactly 86400 due to rounding). */
    while (e->seconds_of_day >= 86400.0) {
        e->seconds_of_day -= 86400.0;
        e->days_since_J2000++;
    }
    while (e->seconds_of_day < 0.0) {
        e->seconds_of_day += 86400.0;
        e->days_since_J2000--;
    }
}

/* ---- Leap-second table --------------------------------------- *
 *
 * Each entry: { JD of the day on which TAI-UTC took effect, value }.
 * The table is keyed by the JD of the START of the UTC day on which
 * the new offset applies (so 1972-01-01T00:00:00 UTC at JD 2441317.5
 * inaugurates TAI-UTC=10s).
 *
 * History: pre-1972, UTC drifted continuously; the table starts at
 * the 1972 reform when leap seconds became the discontinuity model.
 * Post-2017 no leap seconds have been added; the 28th CGPM (2026)
 * resolved to phase out leap seconds entirely by 2035, but the
 * historical entries remain authoritative. */
typedef struct { double jd_start; int dtai; } LeapEntry;

static const LeapEntry LEAP_TABLE[] = {
    { 2441317.5, 10 },   /* 1972-01-01 */
    { 2441499.5, 11 },   /* 1972-07-01 */
    { 2441683.5, 12 },   /* 1973-01-01 */
    { 2442048.5, 13 },   /* 1974-01-01 */
    { 2442413.5, 14 },   /* 1975-01-01 */
    { 2442778.5, 15 },   /* 1976-01-01 */
    { 2443144.5, 16 },   /* 1977-01-01 */
    { 2443509.5, 17 },   /* 1978-01-01 */
    { 2443874.5, 18 },   /* 1979-01-01 */
    { 2444239.5, 19 },   /* 1980-01-01 */
    { 2444786.5, 20 },   /* 1981-07-01 */
    { 2445151.5, 21 },   /* 1982-07-01 */
    { 2445516.5, 22 },   /* 1983-07-01 */
    { 2446247.5, 23 },   /* 1985-07-01 */
    { 2447161.5, 24 },   /* 1988-01-01 */
    { 2447892.5, 25 },   /* 1990-01-01 */
    { 2448257.5, 26 },   /* 1991-01-01 */
    { 2448804.5, 27 },   /* 1992-07-01 */
    { 2449169.5, 28 },   /* 1993-07-01 */
    { 2449534.5, 29 },   /* 1994-07-01 */
    { 2450083.5, 30 },   /* 1996-01-01 */
    { 2450630.5, 31 },   /* 1997-07-01 */
    { 2451179.5, 32 },   /* 1999-01-01 */
    { 2453736.5, 33 },   /* 2006-01-01 */
    { 2454832.5, 34 },   /* 2009-01-01 */
    { 2456109.5, 35 },   /* 2012-07-01 */
    { 2457204.5, 36 },   /* 2015-07-01 */
    { 2457754.5, 37 }    /* 2017-01-01 — most recent */
};
static const int LEAP_N = (int)(sizeof(LEAP_TABLE) / sizeof(LEAP_TABLE[0]));

int k26astro_leap_seconds_at(const K26AstroEpoch *utc_epoch)
{
    if (!utc_epoch) return 0;
    double jd = k26astro_epoch_to_jd(utc_epoch);
    /* Pre-1972: not a regime this library claims to support; return
     * the 1972 baseline as the closest sensible answer. */
    if (jd < LEAP_TABLE[0].jd_start) return LEAP_TABLE[0].dtai;
    /* Linear scan suffices — table is ~30 entries. */
    int dtai = LEAP_TABLE[0].dtai;
    for (int i = 0; i < LEAP_N; i++) {
        if (LEAP_TABLE[i].jd_start <= jd) dtai = LEAP_TABLE[i].dtai;
        else break;
    }
    return dtai;
}

/* ---- DUT1 stub ------------------------------------------------ */

double k26astro_dut1_at(const K26AstroEpoch *utc_epoch)
{
    (void)utc_epoch;
    /* v0.1: zero. The runtime layer can wrap a real IERS series. */
    return 0.0;
}

/* ---- TT - TDB (FB1990 series) -------------------------------- *
 *
 * Fairhead & Bretagnon 1990 single-term approximation, sufficient
 * to ~10 μs accuracy across 1900-2100:
 *
 *   TDB - TT ≈ 0.001658 * sin(g) + 0.000014 * sin(2g)
 *
 * where g = (357.53 + 0.9856003 * d) * pi/180, and d is days from
 * J2000 in TT. The runtime layer can replace this with a TE405
 * Chebyshev series for sub-nanosecond accuracy when a regime demands
 * it. */
static double tdb_minus_tt_seconds_(const K26AstroEpoch *tt)
{
    double d = (double)tt->days_since_J2000 + tt->seconds_of_day / 86400.0;
    double g_rad = (357.53 + 0.9856003 * d) * K26A_RAD_PER_DEG;
    return 0.001658 * sin(g_rad) + 0.000014 * sin(2.0 * g_rad);
}

/* ---- Time-scale conversion ----------------------------------- */

/* Convert e in place: e->scale → target. Returns 0 on success.
 *
 * Pairwise conversions go through TAI as the canonical hub for
 * UTC/UT1/TT (and via TT for TDB). The graph:
 *
 *   UT1 ↔ UTC ↔ TAI ↔ TT ↔ TDB
 *
 * Each adjacent edge is a small adjustment; we walk the graph rather
 * than burning code on every direct pair. */
static int convert_via_tai_(K26AstroEpoch *e, K26AstroTimeScale target);

int k26astro_epoch_convert(K26AstroEpoch *e, K26AstroTimeScale target)
{
    if (!e) return 1;
    if (e->scale == (uint8_t)target) return 0;
    return convert_via_tai_(e, target);
}

/* Single-step adjacent-scale conversion. After this routine,
 * e->scale equals target (if the requested step is adjacent), and
 * the seconds-of-day / day fields carry the offset. */
static int step_(K26AstroEpoch *e, K26AstroTimeScale target)
{
    K26AstroTimeScale cur = (K26AstroTimeScale)e->scale;

    if (cur == K26A_TS_UT1 && target == K26A_TS_UTC) {
        /* UTC = UT1 - DUT1.  We re-evaluate DUT1 at the UT1 epoch
         * because the IERS series is keyed by UTC, but the
         * difference at boundary effects is negligible (DUT1 is
         * smooth, < 1 s). */
        k26astro_epoch_add_seconds(e, -k26astro_dut1_at(e));
        e->scale = K26A_TS_UTC;
        return 0;
    }
    if (cur == K26A_TS_UTC && target == K26A_TS_UT1) {
        k26astro_epoch_add_seconds(e, k26astro_dut1_at(e));
        e->scale = K26A_TS_UT1;
        return 0;
    }
    if (cur == K26A_TS_UTC && target == K26A_TS_TAI) {
        k26astro_epoch_add_seconds(e, (double)k26astro_leap_seconds_at(e));
        e->scale = K26A_TS_TAI;
        return 0;
    }
    if (cur == K26A_TS_TAI && target == K26A_TS_UTC) {
        /* The leap-second count needs to be evaluated at the UTC
         * epoch — but we only have TAI. The two differ by the very
         * count we want. Solve by iterating once: estimate UTC from
         * the TAI epoch using the table keyed by TAI, then refine
         * if the boundary moved. In practice the table boundaries
         * are at integer-second offsets so one pass converges. */
        K26AstroEpoch probe = *e;
        probe.scale = K26A_TS_UTC;
        int dtai = k26astro_leap_seconds_at(&probe);
        k26astro_epoch_add_seconds(e, -(double)dtai);
        e->scale = K26A_TS_UTC;
        return 0;
    }
    if (cur == K26A_TS_TAI && target == K26A_TS_TT) {
        k26astro_epoch_add_seconds(e, K26A_TT_TAI_OFFSET_S);
        e->scale = K26A_TS_TT;
        return 0;
    }
    if (cur == K26A_TS_TT && target == K26A_TS_TAI) {
        k26astro_epoch_add_seconds(e, -K26A_TT_TAI_OFFSET_S);
        e->scale = K26A_TS_TAI;
        return 0;
    }
    if (cur == K26A_TS_TT && target == K26A_TS_TDB) {
        k26astro_epoch_add_seconds(e, tdb_minus_tt_seconds_(e));
        e->scale = K26A_TS_TDB;
        return 0;
    }
    if (cur == K26A_TS_TDB && target == K26A_TS_TT) {
        /* Iterate once — the TDB-TT term is tiny so the fixed point
         * converges in a single pass. */
        K26AstroEpoch probe = *e;
        probe.scale = K26A_TS_TT;
        double delta = tdb_minus_tt_seconds_(&probe);
        k26astro_epoch_add_seconds(e, -delta);
        e->scale = K26A_TS_TT;
        return 0;
    }
    return 2;   /* unsupported adjacent step */
}

static int convert_via_tai_(K26AstroEpoch *e, K26AstroTimeScale target)
{
    /* Walk the chain UT1 ↔ UTC ↔ TAI ↔ TT ↔ TDB in the direction
     * of `target`. Each pass advances the scale by one. */
    static const K26AstroTimeScale chain[] = {
        K26A_TS_UT1, K26A_TS_UTC, K26A_TS_TAI, K26A_TS_TT, K26A_TS_TDB
    };
    static const int N = (int)(sizeof(chain) / sizeof(chain[0]));
    int cur_idx = -1, tgt_idx = -1;
    for (int i = 0; i < N; i++) {
        if (chain[i] == (K26AstroTimeScale)e->scale) cur_idx = i;
        if (chain[i] == target) tgt_idx = i;
    }
    if (cur_idx < 0 || tgt_idx < 0) return 3;
    while (cur_idx != tgt_idx) {
        int step = (tgt_idx > cur_idx) ? +1 : -1;
        if (step_(e, chain[cur_idx + step]) != 0) return 4;
        cur_idx += step;
    }
    return 0;
}
