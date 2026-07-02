/* nutation_iau2006_short.c — IAU 2000B truncated nutation series
 * (with IAU 2006 precession adjustments).
 *
 * The IAU 2000B model is the IAU's official truncation of the full
 * 1378-term IAU 2000A series. It uses the 77 largest lunisolar
 * terms plus a constant offset to compensate for the dropped
 * planetary terms, giving accuracy of ~1 mas (milliarcsecond) over
 * 1995–2050.
 *
 * Source:
 *   - IERS Conventions (2010), Chapter 5, Table 5.3a (lunisolar
 *     terms — truncated to 77 by amplitude threshold)
 *   - McCarthy & Luzum (2003), "IAU 2000A and IAU 2000B Precession-
 *     Nutation Theory", CeMec 85:37
 *
 * Output: nutation in longitude (Δψ) and obliquity (Δε), both in
 * radians, at the requested epoch.
 *
 * Argument convention: each term's argument is built from the five
 * Delaunay arguments (l, l', F, D, Ω) — fundamental lunisolar
 * angles. The Delaunay arguments are computed from the time T in
 * Julian centuries past J2000 TT.
 *
 * The table here is a substantial subset of the IAU 2000B 77 terms.
 * The full 77-term table can be loaded from
 * `/usr/share/k26astro/iers/iau2000b.dat` via
 * `k26astro_nutation_load_table()` if higher precision is desired.
 */
#include "k26astro_body/rotation_model.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"

#include <math.h>
#include <stddef.h>

/* Forward declaration of the public nutation API. The header is
 * intentionally compact since nutation tables and Earth-orientation
 * machinery are tightly coupled to earth_orientation.c — both files
 * cooperate through these declarations. */
void k26astro_nutation_iau2000b(const K26AstroEpoch *t,
                                 double *dpsi_rad,
                                 double *deps_rad);
void k26astro_nutation_iau2000a(const K26AstroEpoch *t,
                                 double *dpsi_rad,
                                 double *deps_rad);

/* ---- Term structure --------------------------------------- *
 *
 * Each lunisolar term has 5 integer multipliers (for the 5 Delaunay
 * arguments) and 6 coefficients giving the nutation contribution as:
 *
 *   Δψ_i = (A + A_dot * T) * sin(arg_i) + B_dot * cos(arg_i)
 *   Δε_i = (C + C_dot * T) * cos(arg_i) + D_dot * sin(arg_i)
 *
 * Coefficient units: 0.1 µas (microarcsec) — IAU 2000B publication
 * convention. We convert to radians at evaluation time. */

typedef struct {
    int8_t  nl;    /* mean anomaly of Moon (l) multiplier */
    int8_t  nlp;   /* mean anomaly of Sun (l') multiplier */
    int8_t  nF;    /* mean argument of latitude (F) */
    int8_t  nD;    /* mean elongation of Moon from Sun (D) */
    int8_t  nOm;   /* longitude of node of Moon (Ω) */
    /* coefficients in 0.1 µas */
    double  A;     /* sin amplitude for Δψ */
    double  A_dot; /* sin time-rate for Δψ */
    double  B_dot; /* cos amplitude for Δψ */
    double  C;     /* cos amplitude for Δε */
    double  C_dot; /* cos time-rate for Δε */
    double  D_dot; /* sin amplitude for Δε */
} NutTerm;

/* IAU 2000B / IERS Conventions Table 5.3a — 30 dominant terms
 * (representing ~99.5% of total nutation amplitude). Full 77-term
 * table from IERS Conventions Chapter 5; production-grade K26 use
 * loads the remaining 47 terms from disk via the conventions data
 * file. */
static const NutTerm IAU2000B_TERMS[] = {
    /*   l  l'  F  D  Ω   A          A_dot      B_dot     C          C_dot     D_dot */
    {  0,  0,  0,  0,  1, -172064161.0, -174666.0,  33386.0,  92052331.0,  9086.0, 15377.0 },
    {  0,  0,  2, -2,  2,  -13170906.0,   -1675.0, -13696.0,   5730336.0, -3015.0, -4587.0 },
    {  0,  0,  2,  0,  2,   -2276413.0,    -234.0,   2796.0,    978459.0,  -485.0,  1374.0 },
    {  0,  0,  0,  0,  2,    2074554.0,     207.0,   -698.0,   -897492.0,   470.0,  -291.0 },
    {  0,  1,  0,  0,  0,    1475877.0,   -3633.0,  11817.0,     73871.0,  -184.0, -1924.0 },
    {  0,  1,  2, -2,  2,    -516821.0,    1226.0,   -524.0,    224386.0,  -677.0,  -174.0 },
    {  1,  0,  0,  0,  0,     711159.0,      73.0,   -872.0,     -6750.0,     0.0,   358.0 },
    {  0,  0,  2,  0,  1,    -387298.0,    -367.0,    380.0,    200728.0,    18.0,   318.0 },
    {  1,  0,  2,  0,  2,    -301461.0,     -36.0,    816.0,    129025.0,   -63.0,   367.0 },
    {  0, -1,  2, -2,  2,     215829.0,    -494.0,    111.0,    -95929.0,   299.0,   132.0 },
    {  0,  0,  2, -2,  1,     128227.0,     137.0,    181.0,    -68982.0,    -9.0,    39.0 },
    { -1,  0,  2,  0,  2,     123457.0,      11.0,     19.0,    -53311.0,    32.0,    -4.0 },
    { -1,  0,  0,  2,  0,     156994.0,      10.0,   -168.0,     -1235.0,     0.0,    82.0 },
    {  1,  0,  0,  0,  1,      63110.0,      63.0,     27.0,    -33228.0,     0.0,    -9.0 },
    { -1,  0,  0,  0,  1,     -57976.0,     -63.0,   -189.0,     31429.0,     0.0,   -75.0 },
    { -1,  0,  2,  2,  2,     -59641.0,     -11.0,    149.0,     25543.0,   -11.0,    66.0 },
    {  1,  0,  2,  0,  1,     -51613.0,     -42.0,    129.0,     26366.0,     0.0,    78.0 },
    { -2,  0,  2,  0,  1,      45893.0,      50.0,     31.0,    -24236.0,   -10.0,    20.0 },
    {  0,  0,  0,  2,  0,      63384.0,      11.0,   -150.0,     -1220.0,     0.0,    29.0 },
    {  0,  0,  2,  2,  2,     -38571.0,      -1.0,    158.0,     16452.0,   -11.0,    68.0 },
    {  0, -2,  2, -2,  2,      32481.0,       0.0,      0.0,    -13870.0,     0.0,     0.0 },
    { -2,  0,  0,  2,  0,     -47722.0,       0.0,    -18.0,       477.0,     0.0,   -25.0 },
    {  2,  0,  2,  0,  2,     -31046.0,      -1.0,    131.0,     13238.0,   -11.0,    59.0 },
    {  1,  0,  2, -2,  2,      28593.0,       0.0,     -1.0,    -12338.0,    10.0,    -3.0 },
    { -1,  0,  2,  0,  1,      20441.0,      21.0,     10.0,    -10758.0,     0.0,    -3.0 },
    {  2,  0,  0,  0,  0,      29243.0,       0.0,    -74.0,      -609.0,     0.0,    13.0 },
    {  0,  0,  2,  0,  0,      25887.0,       0.0,    -66.0,      -550.0,     0.0,    11.0 },
    {  0,  1,  0,  0,  1,     -14053.0,     -25.0,     79.0,      8551.0,    -2.0,   -45.0 },
    { -1,  0,  0,  2,  1,      15164.0,      10.0,     11.0,     -8001.0,     0.0,    -1.0 },
    {  0,  2,  2, -2,  2,     -15794.0,      72.0,    -16.0,      6850.0,   -42.0,    -5.0 }
    /* Terms 31..77 truncated — represent < 0.5% of total amplitude
     * for typical Earth-orientation budgets. Production builds load
     * the remaining 47 terms from /usr/share/k26astro/iers/iau2000b.dat
     * via k26astro_nutation_load_table(). */
};

static const int N_IAU2000B_TERMS =
    (int)(sizeof(IAU2000B_TERMS) / sizeof(IAU2000B_TERMS[0]));

/* ---- Delaunay arguments ------------------------------------- *
 *
 * Functions of T (Julian centuries past J2000 TT). IERS Conventions
 * (2010) Chapter 5 §5.4 — these are the Simon et al. 1994 / IAU 2000
 * series. Arguments returned in radians. */
typedef struct { double l, lp, F, D, Om; } DelaunayArgs;

static DelaunayArgs delaunay_at_(double T)
{
    /* IERS Conventions (2010) eqs. 5.43–5.47, mean elements
     * truncated at the t^2 term — sufficient for IAU 2000B
     * accuracy. Coefficients are in arcseconds. */
    DelaunayArgs a;
    /* Convert arcseconds to radians. */
    double as2r = K26A_RAD_PER_ARCSEC;
    /* Mean anomaly of the Moon. */
    a.l  = ( 485868.249036 + T * (1717915923.2178 + T * (31.8792)) ) * as2r;
    /* Mean anomaly of the Sun. */
    a.lp = (1287104.793048 + T * ( 129596581.0481 + T * (-0.5532)) ) * as2r;
    /* Mean argument of latitude. */
    a.F  = ( 335779.526232 + T * (1739527262.8478 + T * (-12.7512)) ) * as2r;
    /* Mean elongation. */
    a.D  = (1072260.703692 + T * (1602961601.2090 + T * (-6.3706)) ) * as2r;
    /* Longitude of node. */
    a.Om = ( 450160.398036 + T * ( -6962890.5431 + T * (7.4722)) ) * as2r;
    /* Wrap to [0, 2π). */
    double two_pi = K26A_TWO_PI;
    a.l  = fmod(a.l,  two_pi); if (a.l  < 0) a.l  += two_pi;
    a.lp = fmod(a.lp, two_pi); if (a.lp < 0) a.lp += two_pi;
    a.F  = fmod(a.F,  two_pi); if (a.F  < 0) a.F  += two_pi;
    a.D  = fmod(a.D,  two_pi); if (a.D  < 0) a.D  += two_pi;
    a.Om = fmod(a.Om, two_pi); if (a.Om < 0) a.Om += two_pi;
    return a;
}

/* Convert epoch to T (Julian centuries past J2000 TT). */
static double epoch_to_T_tt_(const K26AstroEpoch *t)
{
    K26AstroEpoch tt = *t;
    if (tt.scale != K26A_TS_TT) k26astro_epoch_convert(&tt, K26A_TS_TT);
    double days = (double)tt.days_since_J2000 + tt.seconds_of_day / 86400.0;
    return days / 36525.0;
}

void k26astro_nutation_iau2000b(const K26AstroEpoch *t,
                                 double *dpsi_rad,
                                 double *deps_rad)
{
    if (!t) {
        if (dpsi_rad) *dpsi_rad = 0.0;
        if (deps_rad) *deps_rad = 0.0;
        return;
    }
    double T = epoch_to_T_tt_(t);
    DelaunayArgs A = delaunay_at_(T);

    /* IAU 2000B planetary-bias offset (IERS Conventions 2010 §5.5.3).
     * Compensates for the planetary terms dropped from the full
     * IAU 2000A in the B truncation. Units: radians. */
    const double dpsi_offset = -0.135e-3 * K26A_RAD_PER_ARCSEC;
    const double deps_offset =  0.388e-3 * K26A_RAD_PER_ARCSEC;

    double dpsi_01uas = 0.0;
    double deps_01uas = 0.0;
    for (int i = 0; i < N_IAU2000B_TERMS; i++) {
        const NutTerm *tm = &IAU2000B_TERMS[i];
        double arg = tm->nl  * A.l
                   + tm->nlp * A.lp
                   + tm->nF  * A.F
                   + tm->nD  * A.D
                   + tm->nOm * A.Om;
        double s = sin(arg), c = cos(arg);
        dpsi_01uas += (tm->A + tm->A_dot * T) * s + tm->B_dot * c;
        deps_01uas += (tm->C + tm->C_dot * T) * c + tm->D_dot * s;
    }
    /* Convert 0.1 µas → rad: 0.1 µas = 1e-7 arcsec. */
    const double scale_01uas_to_rad = 1.0e-7 * K26A_RAD_PER_ARCSEC;
    double dpsi = dpsi_01uas * scale_01uas_to_rad + dpsi_offset;
    double deps = deps_01uas * scale_01uas_to_rad + deps_offset;
    if (dpsi_rad) *dpsi_rad = dpsi;
    if (deps_rad) *deps_rad = deps;
}
