/* k26astro_core/epoch.h — precision-preserving epoch type +
 * time-scale arithmetic.
 *
 * Why this is not just `double seconds_since_epoch`
 * -------------------------------------------------
 * A century of integrator substeps at 1-minute resolution is ~5.3
 * million seconds — already at the binary64 precision wall for a
 * 1-microsecond resolution start. Long-running astronomical integration
 * lives in five separate time scales (TAI, UTC, UT1, TT, TDB) and
 * mixes them constantly (Earth-based ephemerides are TT-based; JPL
 * DE441 ephemerides are TDB-based; UTC is the civil time the user
 * sees; UT1 drives sidereal time, hence ECI↔ECEF). A single
 * timestamp must be cheap to convert between scales without
 * accumulating drift over centuries.
 *
 * Internal layout
 * ---------------
 *   int64 days_since_J2000   — signed; J2000 itself is days=0
 *   double seconds_of_day    — in [0, 86400); 0 at J2000-aligned noon
 *   uint8 scale              — K26A_TS_* enum
 *
 * Day-boundary convention: the integer day index ticks over at noon
 * TT (matching the Julian Date convention). J2000 itself
 * (2000-01-01T12:00:00 TT) is days=0, seconds=0. The next noon is
 * days=1, seconds=0; the midnight at 2000-01-02 (i.e. 12 hours later)
 * is days=1, seconds=43200.
 *
 * This convention preserves bit-equivalence with JD-based ephemerides
 * (JD = 2451545.0 + days + seconds/86400) and avoids the half-day-
 * offset bug that any "midnight as day boundary" choice introduces
 * when the time scale is TT or TDB.
 */
#ifndef K26ASTRO_CORE_EPOCH_H
#define K26ASTRO_CORE_EPOCH_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Time scale identifiers ------------------------------------- */
typedef enum {
    K26A_TS_TAI = 0,    /* International Atomic Time — monotonic */
    K26A_TS_UTC = 1,    /* Coordinated Universal Time — has leap seconds */
    K26A_TS_UT1 = 2,    /* Earth-rotation time — tabulated ΔUT1 from IERS */
    K26A_TS_TT  = 3,    /* Terrestrial Time — TAI + 32.184 s exact */
    K26A_TS_TDB = 4     /* Barycentric Dynamical Time — used by JPL DE   */
} K26AstroTimeScale;

/* ---- Epoch ------------------------------------------------------- */
typedef struct {
    int64_t days_since_J2000;   /* signed; J2000 itself is days=0  */
    double  seconds_of_day;     /* in [0, 86400) */
    uint8_t scale;              /* K26AstroTimeScale */
    uint8_t _pad[7];            /* explicit; struct size = 24 bytes */
} K26AstroEpoch;

/* J2000.0 as a constant: days=0, seconds=0, TT. */
K26AstroEpoch k26astro_epoch_j2000_tt(void);

/* ---- Calendar (Gregorian) ↔ epoch round-trips ------------------ */
/* Build an epoch from a Gregorian (proleptic for pre-1582) calendar
 * date in the given time scale. Returns 0 on success, non-zero on
 * out-of-range inputs. month is 1-12; day 1-31; hour 0-23; min 0-59;
 * sec [0, 60) — exact 60.0 is rejected (allows leap seconds via UTC
 * scale but the caller must encode separately). */
int k26astro_epoch_from_calendar(K26AstroEpoch *out,
                                  int year, int month, int day,
                                  int hour, int min, double sec,
                                  K26AstroTimeScale scale);

/* Inverse — split an epoch into calendar fields. The scale field is
 * NOT converted; the input is interpreted in its declared scale. */
void k26astro_epoch_to_calendar(const K26AstroEpoch *e,
                                 int *out_year, int *out_month, int *out_day,
                                 int *out_hour, int *out_min, double *out_sec);

/* ---- ISO-8601 parsing + formatting ----------------------------- */
/* Parse `YYYY-MM-DDTHH:MM:SS[.fff]` (no time-zone suffix — the scale
 * argument supplies that). Returns 0 on success, non-zero on parse
 * failure. */
int k26astro_epoch_from_iso8601(K26AstroEpoch *out,
                                 const char *iso,
                                 K26AstroTimeScale scale);

/* Format as `YYYY-MM-DDTHH:MM:SS.SSSSSS` (microsecond precision).
 * `bufsz` must be at least 27 bytes. Writes a NUL-terminated string.
 * Returns the number of characters written (excluding NUL) on
 * success, or -1 on undersized buffer. */
int k26astro_epoch_to_iso8601(const K26AstroEpoch *e,
                               char *buf, size_t bufsz);

/* ---- Julian Date conversions ----------------------------------- */
/* JD = 2451545.0 + days_since_J2000 + seconds_of_day / 86400.
 * Returned as a double — loses ~10μs precision past 2300 CE because
 * the magnitude grows past 2,600,000 — use the (int64, double)
 * representation for storage. */
double k26astro_epoch_to_jd(const K26AstroEpoch *e);
K26AstroEpoch k26astro_epoch_from_jd(double jd, K26AstroTimeScale scale);

/* ---- Difference / addition ------------------------------------- */
/* Signed `a - b` in seconds. The two epochs MUST be in the same
 * scale; the function asserts this in debug builds. For mixed scales,
 * convert one to the other's scale first. */
double k26astro_epoch_diff_seconds(const K26AstroEpoch *a,
                                    const K26AstroEpoch *b);

/* Add a signed seconds delta to an epoch in place. Carries the
 * seconds-of-day overflow into the integer day count. */
void k26astro_epoch_add_seconds(K26AstroEpoch *e, double delta_s);

/* ---- Time-scale conversion ------------------------------------- */
/* Convert `e` in place from its current scale to `target`. Returns 0
 * on success; non-zero if a needed conversion is unsupported (e.g.
 * UT1↔UTC when no IERS DUT1 table has been provided — falls back to
 * the v0.1 default of DUT1 = 0). */
int k26astro_epoch_convert(K26AstroEpoch *e, K26AstroTimeScale target);

/* Inquire the integer offset TAI - UTC (the leap-second count) at
 * the moment of `utc_epoch`. Used by the UTC↔TAI conversion and by
 * UTC-aware downstream code. */
int k26astro_leap_seconds_at(const K26AstroEpoch *utc_epoch);

/* ---- DUT1 ------------------------------------------------------- */
/* UT1 - UTC offset in seconds at the given epoch. v0.1 returns 0.0
 * unless the caller has registered a DUT1 series (the runtime
 * loads IERS Bulletin A / C04 data). Future API:
 *   int k26astro_epoch_register_dut1_series(const K26AstroDut1Pt *pts,
 *                                            size_t n);
 * Not in v0.1; the stub returns 0.0 and the conversion is no-op. */
double k26astro_dut1_at(const K26AstroEpoch *utc_epoch);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CORE_EPOCH_H */
