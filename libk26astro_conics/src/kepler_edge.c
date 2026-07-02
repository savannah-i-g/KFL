/* kepler_edge.c — Kepler propagator edge cases.
 *
 * Implementation notes:
 *
 *  Parabolic (Barker): the universal-variable formulation degenerates
 *  to a cubic in D = tan(ν/2) at α = 0. Barker's equation states
 *
 *    M_p = √(μ/2p³)·(t - t_peri) = D/2 + D³/6
 *
 *  Solve for D via Cardano's formula (one real root); recover r, ν,
 *  then transform back to inertial pos/vel via Lagrange f/g.
 *
 *  Rectilinear: treat as 1-D Kepler. The "orbit" is a radial line;
 *  position is r(t) along that line, velocity along the same axis.
 *  For elliptic-radial (E < 0): standard Kepler's equation
 *  E - sin E = M with eccentricity 1; for hyperbolic-radial,
 *  hyperbolic Kepler form.
 *
 *  Auto-dispatch: bracket α and h to choose the right path; default
 *  to the baseline propagator and only invoke edge cases when the
 *  basic-orbit detector flags the geometry. */
#include "k26astro_conics/kepler_edge.h"
#include "k26astro_conics/kepler.h"

#include <math.h>

/* ---- V3 helpers (shared with kepler.c by intent; both files keep
 *      their own copies to keep the linker happy without re-exposing
 *      via libk26m3d) ---------------------------------------------- */
static double v_dot_(K26V3 a, K26V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static K26V3  v_cross_(K26V3 a, K26V3 b) {
    K26V3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}
static double v_norm_(K26V3 a) { return sqrt(v_dot_(a, a)); }
static K26V3  v_scale_(K26V3 a, double s) {
    K26V3 r = { a.x * s, a.y * s, a.z * s };
    return r;
}
static K26V3  v_add_(K26V3 a, K26V3 b) {
    K26V3 r = { a.x + b.x, a.y + b.y, a.z + b.z };
    return r;
}

/* ---- Parabolic via Barker ---------------------------------------- *
 *
 * For a parabolic orbit, p = h²/μ, r = p/(1 + cos ν). At periapsis
 * (ν = 0), r = p/2. We need to:
 *   1. Determine the time of periapsis relative to t=0.
 *   2. Add dt, get the new M_p, solve Barker for D = tan(ν/2).
 *   3. Reconstruct (r, ν), build f/g, project back to 3D.
 */
int k26astro_kepler_propagate_parabolic(K26V3 *out_pos, K26V3 *out_vel,
                                        K26V3 pos0, K26V3 vel0,
                                        double mu, double dt)
{
    if (!out_pos || !out_vel) return 1;
    if (mu <= 0.0) return 2;

    K26V3 h_vec = v_cross_(pos0, vel0);
    double h = v_norm_(h_vec);
    if (h <= 0.0) return 5;   /* degenerate; caller should use radial path */

    double p = h * h / mu;
    double r0 = v_norm_(pos0);
    if (r0 <= 0.0) return 3;

    /* True anomaly at t=0 from r0 = p/(1 + cos ν0) */
    double cos_nu0 = p / r0 - 1.0;
    if (cos_nu0 >  1.0) cos_nu0 =  1.0;
    if (cos_nu0 < -1.0) cos_nu0 = -1.0;
    double nu0 = acos(cos_nu0);
    /* Sign from r·v: positive r-dot → past periapsis */
    if (v_dot_(pos0, vel0) < 0.0) nu0 = -nu0;

    double D0 = tan(nu0 / 2.0);
    /* M_p(t) = D/2 + D³/6, so M_p(t=0) = D0/2 + D0³/6 */
    double Mp0 = D0 / 2.0 + D0 * D0 * D0 / 6.0;

    /* M_p advances by (μ / (2p³))^(1/2) · dt */
    double n_p = sqrt(mu / (2.0 * p * p * p));
    double Mp1 = Mp0 + n_p * dt;

    /* Solve D/2 + D³/6 = Mp1 via Cardano:
     *   D³ + 3D - 6 Mp1 = 0
     *   discriminant = 9 + 9 Mp1² (always positive)
     *   D = ³√(3 Mp1 + √(9 Mp1² + 9)) + ³√(3 Mp1 - √(9 Mp1² + 9)) */
    double disc = sqrt(9.0 * Mp1 * Mp1 + 9.0);
    double t1 = 3.0 * Mp1 + disc;
    double t2 = 3.0 * Mp1 - disc;
    double D1 = cbrt(t1) + cbrt(t2);

    double nu1 = 2.0 * atan(D1);
    double r1  = p / (1.0 + cos(nu1));

    /* Lagrange f/g from Vallado §2.3 (parabolic) */
    double dnu = nu1 - nu0;
    double f = 1.0 - r1 / p * (1.0 - cos(dnu));
    double g = r1 * r0 / sqrt(mu * p) * sin(dnu);
    double fdot = sqrt(mu / p) * tan(dnu / 2.0)
                * ((1.0 - cos(dnu)) / p - 1.0 / r0 - 1.0 / r1);
    double gdot = 1.0 - r0 / p * (1.0 - cos(dnu));

    *out_pos = v_add_(v_scale_(pos0, f),    v_scale_(vel0, g));
    *out_vel = v_add_(v_scale_(pos0, fdot), v_scale_(vel0, gdot));
    return 0;
}

/* ---- Rectilinear orbit ------------------------------------------- *
 *
 * If h = r × v ≈ 0, position and velocity are colinear: the body is on
 * a radial trajectory. Reduce to 1-D: r(t) along the line direction
 * with v(t) ditto.
 *
 * Energy E = v²/2 - μ/r determines whether the orbit is elliptic-
 * radial (E < 0; falls back), parabolic-radial (E = 0), or hyperbolic-
 * radial (E > 0; escapes). For elliptic: a = -μ/(2E), period T;
 * solve M = n(t - tp) = E - sin E with e = 1 in the line-frame
 * formulation. For hyperbolic: hyperbolic Kepler. */
int k26astro_kepler_propagate_radial(K26V3 *out_pos, K26V3 *out_vel,
                                     K26V3 pos0, K26V3 vel0,
                                     double mu, double dt)
{
    if (!out_pos || !out_vel) return 1;
    if (mu <= 0.0) return 2;

    double r0 = v_norm_(pos0);
    if (r0 <= 0.0) return 3;

    /* Line direction (unit vector). */
    K26V3 dir = v_scale_(pos0, 1.0 / r0);

    /* Sign of v0 along dir. */
    double v0_signed = v_dot_(vel0, dir);

    /* Specific energy. */
    double v2 = v_dot_(vel0, vel0);
    double E  = 0.5 * v2 - mu / r0;

    double r1, v1_signed;

    if (E < -1.0e-12 * (mu / r0)) {
        /* Elliptic-radial. a = -μ/(2E). Kepler eq: M = E_anom - sin E_anom
         * with e = 1 (rectilinear); but use universal-variable approach
         * recast to 1-D for stability. The cleanest form: integrate
         *   dr/dt = ±√(2μ/r - μ/a)
         * via a Newton step on the analytic inverse. For mid-range
         * (~10²−10⁵ s, inner-SS scale) we accept an iterative impl. */
        double a = -mu / (2.0 * E);
        double n = sqrt(mu / (a * a * a));
        /* Eccentric anomaly form: r = a(1 - cos E_anom),
         * t = (E_anom - sin E_anom)/n. Initial E_anom from r0. */
        double cosE0 = 1.0 - r0 / a;
        if (cosE0 >  1.0) cosE0 =  1.0;
        if (cosE0 < -1.0) cosE0 = -1.0;
        double E0 = acos(cosE0);
        if (v0_signed < 0.0) E0 = -E0;    /* falling: negative branch */
        double M0 = E0 - sin(E0);
        double M1 = M0 + n * dt;
        /* Solve E1 - sin E1 = M1 (Kepler with e = 1). */
        double E1 = M1;
        for (int it = 0; it < 64; it++) {
            double f  = E1 - sin(E1) - M1;
            double fp = 1.0 - cos(E1);
            if (fabs(fp) < 1e-16) break;
            double dE = f / fp;
            E1 -= dE;
            if (fabs(dE) < 1e-14) break;
        }
        r1 = a * (1.0 - cos(E1));
        v1_signed = sqrt(2.0 * mu / r1 - mu / a);
        if (sin(E1) < 0.0) v1_signed = -v1_signed;
    } else if (E > 1.0e-12 * (mu / r0)) {
        /* Hyperbolic-radial. */
        double a = -mu / (2.0 * E);     /* negative */
        double n = sqrt(mu / (-a * a * a));
        /* H: r = -a(cosh H - 1); t = (sinh H - H) / n */
        double coshH0 = 1.0 - r0 / a;    /* -a > 0; r0/a < 0 */
        if (coshH0 < 1.0) coshH0 = 1.0;
        double H0 = acosh(coshH0);
        if (v0_signed < 0.0) H0 = -H0;
        double M0 = sinh(H0) - H0;
        double M1 = M0 + n * dt;
        double H1 = M1;
        for (int it = 0; it < 64; it++) {
            double f  = sinh(H1) - H1 - M1;
            double fp = cosh(H1) - 1.0;
            if (fabs(fp) < 1e-16) break;
            double dH = f / fp;
            H1 -= dH;
            if (fabs(dH) < 1e-14) break;
        }
        r1 = -a * (cosh(H1) - 1.0);
        v1_signed = sqrt(2.0 * mu / r1 + 2.0 * E);
        if (sinh(H1) < 0.0) v1_signed = -v1_signed;
    } else {
        /* Parabolic-radial. dr/dt = ±√(2μ/r); integrate analytically
         *  ∫ √r dr = ±√(2μ) dt → (2/3) r^(3/2) = ±√(2μ) t + (2/3) r0^(3/2) */
        double sgn = (v0_signed >= 0.0) ? 1.0 : -1.0;
        double rhs = sgn * sqrt(2.0 * mu) * dt + (2.0 / 3.0) * pow(r0, 1.5);
        if (rhs < 0.0) rhs = 0.0;
        r1 = pow(1.5 * rhs, 2.0 / 3.0);
        v1_signed = sgn * sqrt(2.0 * mu / r1);
    }

    *out_pos = v_scale_(dir, r1);
    *out_vel = v_scale_(dir, v1_signed);
    return 0;
}

/* ---- Dispatch ---------------------------------------------------- */
int k26astro_kepler_propagate_any(K26V3 *out_pos, K26V3 *out_vel,
                                  K26V3 pos0, K26V3 vel0,
                                  double mu, double dt,
                                  int max_iter)
{
    K26V3 h_vec = v_cross_(pos0, vel0);
    double h = v_norm_(h_vec);
    double r0 = v_norm_(pos0);
    double v0 = v_norm_(vel0);

    /* Radial: |h| relative to r0·v0 is tiny. */
    if (r0 > 0.0 && v0 > 0.0
        && h < K26A_KEPLER_RADIAL_EPSILON * r0 * v0) {
        return k26astro_kepler_propagate_radial(out_pos, out_vel,
                                                 pos0, vel0, mu, dt);
    }

    /* Parabolic: |α| ≈ 0, i.e. v² ≈ 2μ/r0. */
    double alpha = 2.0 / r0 - v0 * v0 / mu;
    if (fabs(alpha) < K26A_KEPLER_PARABOLIC_EPSILON / r0) {
        return k26astro_kepler_propagate_parabolic(out_pos, out_vel,
                                                    pos0, vel0, mu, dt);
    }

    return k26astro_kepler_propagate(out_pos, out_vel,
                                      pos0, vel0, mu, dt, max_iter);
}
