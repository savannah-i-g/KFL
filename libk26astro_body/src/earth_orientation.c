/* earth_orientation.c — GCRS ↔ ITRS (ECEF) full IAU 2006/2000A pipeline.
 *
 * The Earth's orientation in inertial space evolves under four
 * superposed effects:
 *
 *   1. Precession    — long-period drift of the spin axis (26 ka
 *                      cycle). IAU 2006 (Capitaine et al. 2003).
 *   2. Nutation      — short-period wobble of the spin axis driven
 *                      by lunisolar + planetary torques. Delegated to
 *                      nutation_iau2000a.c (full series) or
 *                      nutation_iau2006_short.c (~1 mas truncation).
 *   3. Rotation      — Earth Rotation Angle (ERA) about the celestial
 *                      intermediate pole. Linear in UT1.
 *   4. Polar motion  — small (~10 m at Earth surface) wobble of the
 *                      pole relative to the body-fixed frame. IERS
 *                      Bulletin A provides daily xp, yp.
 *
 * The full transform from GCRS (= ICRF re-centred at the geocentre)
 * to ITRS (= ECEF, the Earth-fixed reference frame) at TT epoch
 * t_TT and UT1 epoch t_UT1 is:
 *
 *   ITRS = W(t_TT) · R3(ERA(t_UT1)) · Q(t_TT) · GCRS
 *
 * where:
 *   Q(t_TT)  = bias-precession-nutation matrix (CIO-based)
 *   ERA      = Earth Rotation Angle (IAU 2000 / IERS Conv. eq. 5.15)
 *   W(t_TT)  = polar motion matrix (small corrections)
 *
 * Reference:
 *   IERS Conventions (2010), Chapter 5 (Petit & Luzum, eds.).
 *   IAU 2006 Resolution B1 (precession), B2 (nutation).
 *
 * Implementation note:
 *   The CIO-based formulation (using s' and ERA) is the IAU 2006
 *   preferred convention. The earlier IAU 1976/1980 equinox-based
 *   formulation (using GMST and equation of the equinoxes) is
 *   mathematically equivalent at sub-microarcsecond precision but
 *   conceptually denser. K26 uses CIO-based throughout.
 */
#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#include <math.h>
#include <stddef.h>

/* Forward declarations from the nutation modules. */
void k26astro_nutation_iau2000a(const K26AstroEpoch *t,
                                 double *dpsi_rad, double *deps_rad);
void k26astro_nutation_iau2000b(const K26AstroEpoch *t,
                                 double *dpsi_rad, double *deps_rad);

/* ---- Public API ------------------------------------------- */

/* Compute the 3×3 rotation matrix GCRS → ITRS (column-major) at
 * the given epoch. `t` is treated as TT for the precession-nutation
 * pieces and as UT1 for the rotation angle (typical caller has both;
 * we lump them into a single TT-aligned input and add UT1-TT via
 * the DUT1 lookup at epoch).
 *
 * `mode` selects nutation precision:
 *   0 = IAU 2006/2000B (truncated, ~1 mas)
 *   1 = IAU 2006/2000A (full series — uses the active table in
 *       nutation_iau2000a.c which may have been replaced by
 *       k26astro_nutation_load_table) */
void k26astro_earth_orientation_matrix(const K26AstroEpoch *t,
                                       int mode,
                                       double out_R[9]);

/* Apply that matrix to a vector. */
void k26astro_earth_orientation_apply(const K26AstroEpoch *t,
                                      int mode,
                                      const double in_xyz[3],
                                      double out_xyz[3]);

/* Output the cardinal components of the pipeline (precession-nutation
 * matrix Q, ERA, polar motion matrix W) separately — useful for
 * tests + diagnostic comparisons against SOFA reference values. */
void k26astro_earth_orientation_components(const K26AstroEpoch *t,
                                           int mode,
                                           double out_Q[9],
                                           double *out_era,
                                           double out_W[9]);

/* Polar motion coordinate setter. Default is (0, 0) — no wobble.
 * IERS Bulletin A subscribers can update these per call:
 *   k26astro_polar_motion_set(xp_arcsec, yp_arcsec);
 * Coordinates are in arcseconds. */
void k26astro_polar_motion_set(double xp_arcsec, double yp_arcsec);
void k26astro_polar_motion_get(double *out_xp_arcsec, double *out_yp_arcsec);

/* ---- Implementation ---------------------------------------- */

static double g_xp_arcsec = 0.0;
static double g_yp_arcsec = 0.0;

void k26astro_polar_motion_set(double xp_arcsec, double yp_arcsec)
{
    g_xp_arcsec = xp_arcsec;
    g_yp_arcsec = yp_arcsec;
}

void k26astro_polar_motion_get(double *out_xp, double *out_yp)
{
    if (out_xp) *out_xp = g_xp_arcsec;
    if (out_yp) *out_yp = g_yp_arcsec;
}

/* Convert epoch to T (Julian centuries past J2000 TT). */
static double epoch_to_T_tt_(const K26AstroEpoch *t)
{
    K26AstroEpoch tt = *t;
    if (tt.scale != K26A_TS_TT) k26astro_epoch_convert(&tt, K26A_TS_TT);
    double days = (double)tt.days_since_J2000 + tt.seconds_of_day / 86400.0;
    return days / 36525.0;
}

/* ---- Precession (IAU 2006, Capitaine et al. 2003) -------- *
 *
 * Computes the bias-precession quantities X, Y, s' that define the
 * CIO (Celestial Intermediate Origin) location in GCRS.
 *
 * Series (IERS Conv. eq. 5.16, truncated to leading polynomial terms
 * — full series + sin/cos table provides 0.01 mas accuracy; the
 * truncated polynomial gives ~1 mas at present epoch which is
 * consistent with the IAU 2000B nutation truncation):
 *
 *   X = -0.016617 + 2004.191898 t - 0.4297829 t² - 0.19861834 t³
 *                 + 7.578e-6   t⁴ + 5.9285e-6 t⁵                  arcsec
 *   Y = -0.006951 - 0.025896 t  - 22.4072747 t² + 0.00190059 t³
 *                 + 0.001112526 t⁴ + 0.0000001358 t⁵              arcsec
 *
 * Full Capitaine series adds ~1377 trig terms in X and Y. For K26
 * v0.1 the polynomial is sufficient; high-fidelity ground-station
 * tracking applications can switch to the full series via the
 * load-table path. */
static void precession_xy_(double T, double *out_X_rad, double *out_Y_rad)
{
    const double as2r = K26A_RAD_PER_ARCSEC;
    double X_as =  -0.016617
                +  2004.191898  * T
                -  0.4297829    * T * T
                -  0.19861834   * T * T * T
                +  7.578e-6     * T * T * T * T
                +  5.9285e-6    * T * T * T * T * T;
    double Y_as =  -0.006951
                -  0.025896     * T
                - 22.4072747    * T * T
                +  0.00190059   * T * T * T
                +  0.001112526  * T * T * T * T
                +  0.0000001358 * T * T * T * T * T;
    if (out_X_rad) *out_X_rad = X_as * as2r;
    if (out_Y_rad) *out_Y_rad = Y_as * as2r;
}

/* CIO locator s'(t) — sub-microarcsecond term that completes the
 * CIRS frame definition. IERS Conv. eq. 5.18.
 *
 *   s'(t) = -47.0 µas * t              (linear approximation)
 *
 * The full series adds short-period terms below microarcsec. */
static double cio_locator_s_(double T)
{
    return -47.0e-6 * K26A_RAD_PER_ARCSEC * T;
}

/* CIO-based bias-precession-nutation matrix Q.
 *
 *   Q = R3(-s) · R3(-E) · R2(-d) · R3(E)
 *
 * where (X, Y) define the position of the CIP in GCRS and
 *   E = atan2(Y, X), d = asin(sqrt(X² + Y²)).
 *
 * After applying nutation (X, Y get the Δψ, Δε contributions via
 * the standard CIO-based formula), Q maps GCRS → CIRS. */
static void build_Q_matrix_(double T,
                            double dpsi, double deps,
                            double out_Q[9])
{
    /* Precession part of (X, Y). */
    double X_p = 0.0, Y_p = 0.0;
    precession_xy_(T, &X_p, &Y_p);

    /* Mean obliquity of the ecliptic at J2000 (IAU 2006). */
    const double eps0 = 84381.406 * K26A_RAD_PER_ARCSEC;
    /* Nutation contribution to X, Y (IERS Conv. eq. 5.20 simplified). */
    double X = X_p + dpsi * sin(eps0);
    double Y = Y_p + deps;

    double s = cio_locator_s_(T);
    double a = 1.0 / (1.0 + sqrt(1.0 - X * X - Y * Y));
    /* Q matrix per IERS Conv. eq. 5.10 (CIO-based, no equinox). */
    /* Q = (
     *   [ 1 - aX²,   -aXY,           X           ]
     *   [ -aXY,       1 - aY²,        Y           ]
     *   [ -X,        -Y,             1 - a(X²+Y²) ] ) * R3(s)
     */
    double M[9];
    M[0] =  1.0 - a * X * X;
    M[3] = -a * X * Y;
    M[6] =  X;
    M[1] = -a * X * Y;
    M[4] =  1.0 - a * Y * Y;
    M[7] =  Y;
    M[2] = -X;
    M[5] = -Y;
    M[8] =  1.0 - a * (X * X + Y * Y);

    /* Apply R3(s) on the right: column-major, R3(s) is
     *   [ cos s, sin s, 0 ]
     *   [-sin s, cos s, 0 ]
     *   [   0,    0,    1 ]
     * In column-major form, R3 on the right combines as
     *   Q[0] = M[0]*cos(s) - M[3]*sin(s)
     *   Q[3] = M[0]*sin(s) + M[3]*cos(s)
     *   Q[1] = M[1]*cos(s) - M[4]*sin(s)
     *   Q[4] = M[1]*sin(s) + M[4]*cos(s)
     *   Q[2] = M[2]*cos(s) - M[5]*sin(s)
     *   Q[5] = M[2]*sin(s) + M[5]*cos(s)
     *   Q[6] = M[6]
     *   Q[7] = M[7]
     *   Q[8] = M[8]
     */
    double cs = cos(s), ss = sin(s);
    out_Q[0] = M[0] * cs - M[3] * ss;
    out_Q[3] = M[0] * ss + M[3] * cs;
    out_Q[1] = M[1] * cs - M[4] * ss;
    out_Q[4] = M[1] * ss + M[4] * cs;
    out_Q[2] = M[2] * cs - M[5] * ss;
    out_Q[5] = M[2] * ss + M[5] * cs;
    out_Q[6] = M[6];
    out_Q[7] = M[7];
    out_Q[8] = M[8];
}

/* Earth Rotation Angle — IAU 2000.
 *
 *   ERA(UT1) = 2π * (0.7790572732640 + 1.00273781191135448 * D_UT1)
 *
 * where D_UT1 = (UT1 - 2451545.0) in days. */
static double compute_era_(const K26AstroEpoch *t)
{
    K26AstroEpoch ut1 = *t;
    if (ut1.scale != K26A_TS_UT1) k26astro_epoch_convert(&ut1, K26A_TS_UT1);
    double D = (double)ut1.days_since_J2000 + ut1.seconds_of_day / 86400.0;
    /* IERS Conv. eq. 5.15. */
    double era = K26A_TWO_PI * (0.7790572732640 + 1.00273781191135448 * D);
    /* Wrap to [0, 2π). */
    era = fmod(era, K26A_TWO_PI);
    if (era < 0.0) era += K26A_TWO_PI;
    return era;
}

/* Polar motion matrix W = R1(-yp) · R2(-xp) · R3(s')
 *
 * Small angles; s' (the TIO locator) is sub-microarcsec and
 * approximated linearly in T. */
static void build_W_matrix_(double T, double out_W[9])
{
    double xp = g_xp_arcsec * K26A_RAD_PER_ARCSEC;
    double yp = g_yp_arcsec * K26A_RAD_PER_ARCSEC;
    /* TIO locator s'(t) ≈ -47 µas * T — same magnitude as CIO s, in
     * the opposite direction. */
    double sp = -47.0e-6 * K26A_RAD_PER_ARCSEC * T;

    /* R3(s') first: small rotation about z. */
    double css = cos(sp), sss = sin(sp);
    /* R2(-xp): small rotation about y by -xp. */
    double cx = cos(-xp), sx = sin(-xp);
    /* R1(-yp): small rotation about x by -yp. */
    double cy = cos(-yp), sy = sin(-yp);

    /* Compose: W = R1(-yp) * R2(-xp) * R3(s'). For small angles
     * (~arcsecond), the linearised form is accurate to (arcsec)²
     * ≈ 1e-12 — well below working precision. We use the exact
     * trig form here for cleanliness; cost is negligible. */
    /* Step 1: A = R2(-xp) * R3(s')  (3x3) */
    double A[9];
    A[0] =  cx * css;       A[3] =  cx * sss;       A[6] = -sx;
    A[1] =  -sss;           A[4] =   css;           A[7] =  0.0;
    A[2] =  sx * css;       A[5] =  sx * sss;       A[8] =  cx;
    /* Step 2: W = R1(-yp) * A. */
    out_W[0] =  A[0];
    out_W[3] =  A[3];
    out_W[6] =  A[6];
    out_W[1] =  cy * A[1] + sy * A[2];
    out_W[4] =  cy * A[4] + sy * A[5];
    out_W[7] =  cy * A[7] + sy * A[8];
    out_W[2] = -sy * A[1] + cy * A[2];
    out_W[5] = -sy * A[4] + cy * A[5];
    out_W[8] = -sy * A[7] + cy * A[8];
}

/* 3x3 column-major matrix multiplication: C = A * B */
static void mat3_mul_(const double A[9], const double B[9], double C[9])
{
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < 3; row++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) {
                s += A[k * 3 + row] * B[col * 3 + k];
            }
            C[col * 3 + row] = s;
        }
    }
}

void k26astro_earth_orientation_matrix(const K26AstroEpoch *t,
                                       int mode,
                                       double out_R[9])
{
    if (!t || !out_R) return;

    double T = epoch_to_T_tt_(t);

    /* Nutation. */
    double dpsi = 0.0, deps = 0.0;
    if (mode == 1) {
        k26astro_nutation_iau2000a(t, &dpsi, &deps);
    } else {
        k26astro_nutation_iau2000b(t, &dpsi, &deps);
    }

    /* Q matrix (bias-precession-nutation). */
    double Q[9];
    build_Q_matrix_(T, dpsi, deps, Q);

    /* ERA rotation about z. */
    double era = compute_era_(t);
    double cE = cos(era), sE = sin(era);
    double Rerar[9] = {
         cE,  sE, 0.0,
        -sE,  cE, 0.0,
         0.0, 0.0, 1.0
    };

    /* Polar motion matrix W. */
    double W[9];
    build_W_matrix_(T, W);

    /* Final: out_R = W * R3(ERA) * Q. */
    double tmp[9];
    mat3_mul_(Rerar, Q,    tmp);   /* R3(ERA) * Q */
    mat3_mul_(W,    tmp,  out_R);  /* W * (R3 * Q) */
}

void k26astro_earth_orientation_apply(const K26AstroEpoch *t,
                                      int mode,
                                      const double in_xyz[3],
                                      double out_xyz[3])
{
    if (!t || !in_xyz || !out_xyz) return;
    double R[9];
    k26astro_earth_orientation_matrix(t, mode, R);
    out_xyz[0] = R[0] * in_xyz[0] + R[3] * in_xyz[1] + R[6] * in_xyz[2];
    out_xyz[1] = R[1] * in_xyz[0] + R[4] * in_xyz[1] + R[7] * in_xyz[2];
    out_xyz[2] = R[2] * in_xyz[0] + R[5] * in_xyz[1] + R[8] * in_xyz[2];
}

void k26astro_earth_orientation_components(const K26AstroEpoch *t,
                                           int mode,
                                           double out_Q[9],
                                           double *out_era,
                                           double out_W[9])
{
    if (!t) return;
    double T = epoch_to_T_tt_(t);
    double dpsi = 0.0, deps = 0.0;
    if (mode == 1) {
        k26astro_nutation_iau2000a(t, &dpsi, &deps);
    } else {
        k26astro_nutation_iau2000b(t, &dpsi, &deps);
    }
    if (out_Q)   build_Q_matrix_(T, dpsi, deps, out_Q);
    if (out_era) *out_era = compute_era_(t);
    if (out_W)   build_W_matrix_(T, out_W);
}
