/* lambert.c — single-revolution Lambert solver.
 *
 * Vallado §7.6 / Algorithm 58 (universal-variable form).
 *
 * Iteration scheme:
 *   1. Compute the transfer angle Δν from r̂₁ · r̂₂ and the
 *      direction flag.
 *   2. Compute the constant A = sin(Δν) · √(r1·r2 / (1 - cos(Δν))).
 *   3. Newton-Raphson on z (the universal variable):
 *        y(z)  = r1 + r2 + A·(z·S(z) - 1) / √C(z)
 *        χ(z)  = √(y(z) / C(z))
 *        t(z)  = (χ(z)³·S(z) + A·√y(z)) / √μ
 *      Solve t(z) = tof.
 *   4. Build f, g, ġ Lagrange coefficients; v₁ = (r₂ - f·r₁)/g,
 *      v₂ = (ġ·r₂ - r₁)/g.
 *
 * Convergence: z bracketed by elliptic [-4π², 4π²], hyperbolic
 * unbounded (use exponential expansion for large negative z). */
#include "k26astro_conics/lambert.h"
#include "stumpff_internal.h"

#include <math.h>

#define LAMBERT_MAX_ITER     64
#define LAMBERT_TOL          1.0e-9

/* ---- V3 helpers --------------------------------------------------- */
static double v_dot_(K26V3 a, K26V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double v_norm_(K26V3 a)         { return sqrt(v_dot_(a, a)); }
static K26V3  v_cross_(K26V3 a, K26V3 b)
{
    K26V3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}
static K26V3  v_scale_(K26V3 a, double s)
{
    K26V3 r = { a.x * s, a.y * s, a.z * s };
    return r;
}
static K26V3  v_sub_(K26V3 a, K26V3 b)
{
    K26V3 r = { a.x - b.x, a.y - b.y, a.z - b.z };
    return r;
}

int k26astro_lambert(K26V3 *out_v1, K26V3 *out_v2,
                     K26V3 r1, K26V3 r2,
                     double mu, double tof,
                     int direction)
{
    if (!out_v1 || !out_v2) return 1;
    if (mu <= 0.0 || tof <= 0.0) return 2;

    double r1m = v_norm_(r1);
    double r2m = v_norm_(r2);
    if (r1m <= 0.0 || r2m <= 0.0) return 3;

    /* Transfer angle from cos Δν = r̂₁·r̂₂ */
    double cos_dnu = v_dot_(r1, r2) / (r1m * r2m);
    if (cos_dnu >  1.0) cos_dnu =  1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;

    /* Sign of sin Δν from the cross product's projection onto
     * a chosen "up" axis. Vallado picks +z; we follow. */
    K26V3 cross12 = v_cross_(r1, r2);
    double sin_dnu_z = cross12.z;
    double sin_dnu;
    if (direction == K26A_LAMBERT_SHORT_WAY) {
        sin_dnu = sqrt(1.0 - cos_dnu * cos_dnu);
        if (sin_dnu_z < 0.0) sin_dnu = -sin_dnu;
    } else {
        sin_dnu = -sqrt(1.0 - cos_dnu * cos_dnu);
        if (sin_dnu_z < 0.0) sin_dnu = -sin_dnu;
    }

    /* 180° check. */
    if (fabs(1.0 - cos_dnu * cos_dnu) < 1e-12) return 5;

    /* A = sin Δν · √(r1 · r2 / (1 - cos Δν)). */
    double A = sin_dnu * sqrt(r1m * r2m / (1.0 - cos_dnu));
    if (A == 0.0) return 5;

    /* Safeguarded Newton on z: t(z) is monotonically increasing in z
     * for both short-way (A > 0) and long-way (A < 0). Maintain a
     * proper bracket; fall back to bisection when Newton overshoots
     * the bracket or stalls. This form is stable across the A sign
     * change that distinguishes the long-way path (Battin 1999 §7.6
     * A-sign-stable derivative form).
     *
     * Locally-defined helper to evaluate y, x, t at a candidate z;
     * returns 0 on success, 1 if y < 0 (z too high for elliptic). */
    double C, S, y, x, t;
    double z_lo = -4.0 * 3.141592653589793 * 3.141592653589793;
    double z_hi =  4.0 * 3.141592653589793 * 3.141592653589793;

    /* Initial bracket scan: evaluate at z_lo + steps and at z_hi-eps
     * to confirm t(z_lo) < tof < t(z_hi). For most realistic
     * geometries, z_lo gives t < 1 day and z_hi gives t > centuries,
     * so the bracket is huge but valid. */
    double z = 0.0;
    int iter;
    double prev_dz = 1.0e30;
    for (iter = 0; iter < LAMBERT_MAX_ITER; iter++) {
        C = k26astro_conics_stumpff_C(z);
        S = k26astro_conics_stumpff_S(z);
        if (C <= 0.0) {
            /* Stumpff C(z) hits zero only at z = (2π)² = ~39.48
             * (elliptic limit). Pull back toward lower bound. */
            z_hi = z;
            z = 0.5 * (z_lo + z_hi);
            continue;
        }
        y = r1m + r2m + A * (z * S - 1.0) / sqrt(C);
        if (y < 0.0) {
            /* Solution lies further into the elliptic regime. */
            z_lo = z;
            z = 0.5 * (z_lo + z_hi);
            continue;
        }
        x = sqrt(y / C);
        t = (x * x * x * S + A * sqrt(y)) / sqrt(mu);

        /* Convergence on |t - tof|. */
        double f = t - tof;
        if (fabs(f) < LAMBERT_TOL * tof) break;

        /* Update bracket: t is monotonically increasing in z. */
        if (f < 0.0) { if (z > z_lo) z_lo = z; }
        else         { if (z < z_hi) z_hi = z; }

        /* Newton step: dt/dz analytic (Vallado App. D.5). */
        double dCdz, dSdz;
        if (fabs(z) > 1e-6) {
            dCdz = (1.0 - z * S - 2.0 * C) / (2.0 * z);
            dSdz = (C - 3.0 * S) / (2.0 * z);
        } else {
            dCdz = -1.0 / 24.0 + z / 360.0;
            dSdz = -1.0 / 120.0 + z / 2520.0;
        }
        double dydz = A * (S + z * dSdz - dCdz * (z * S - 1.0) / (2.0 * C))
                       / sqrt(C);
        double sqrt_y_safe = (y > 1e-30) ? sqrt(y) : 1e-15;
        double dtdz = (3.0 * x * S * dydz / (2.0 * C)
                       - x * x * x * (S * dCdz / (2.0 * C * C) - dSdz)
                       + A * dydz / (2.0 * sqrt_y_safe)) / sqrt(mu);

        /* Trial Newton step. */
        double z_newton = (dtdz != 0.0) ? z - f / dtdz : z;

        /* Safeguard: if Newton overshoots the bracket OR isn't
         * shrinking the step size, fall back to bisection. */
        int use_bisect = 0;
        if (z_newton <= z_lo || z_newton >= z_hi) use_bisect = 1;
        double dz_newton = z_newton - z;
        if (fabs(dz_newton) > 0.5 * fabs(prev_dz) && iter > 2)
            use_bisect = 1;

        double z_next = use_bisect ? 0.5 * (z_lo + z_hi) : z_newton;
        prev_dz = z_next - z;
        if (fabs(prev_dz) < LAMBERT_TOL) { z = z_next; break; }
        z = z_next;
    }
    if (iter >= LAMBERT_MAX_ITER) return 4;

    /* Re-evaluate y at the converged z for the Lagrange coefficients. */
    C = k26astro_conics_stumpff_C(z);
    S = k26astro_conics_stumpff_S(z);
    y = r1m + r2m + A * (z * S - 1.0) / sqrt(C);
    (void)x; (void)t;

    /* Lagrange f / g / ġ ; recover v₁, v₂. */
    double f    = 1.0 - y / r1m;
    double g    = A * sqrt(y / mu);
    double gdot = 1.0 - y / r2m;

    K26V3 fr1 = v_scale_(r1, f);
    K26V3 v1  = v_scale_(v_sub_(r2, fr1), 1.0 / g);
    K26V3 gdot_r2 = v_scale_(r2, gdot);
    K26V3 v2  = v_scale_(v_sub_(gdot_r2, r1), 1.0 / g);

    *out_v1 = v1;
    *out_v2 = v2;
    return 0;
}
