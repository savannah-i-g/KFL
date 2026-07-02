/* libk26compute — optimisation: 1-D Brent + N-D Nelder-Mead. */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- Brent's method ----
 * Reference: Brent (1973), "Algorithms for Minimization without
 * Derivatives". Standard textbook formulation: parabolic
 * interpolation when safe, golden section as fallback. */

#define K26C_BRENT_GOLDEN_RATIO 0.3819660112501051   /* 1 - 1/φ */
#define K26C_BRENT_TINY         1e-20
#define K26C_BRENT_MAX_ITERS    100

K26CStatus k26c_opt_brent(K26CObj1Fn f, void *user,
                          double lo, double hi, double tol,
                          double *x_min, double *f_min)
{
    if (!f || !x_min || !f_min || !(hi > lo)) return K26C_ERR_INVAL;
    if (tol <= 0.0) tol = 1e-8;

    double a = lo, b = hi;
    double x = a + K26C_BRENT_GOLDEN_RATIO * (b - a);
    double w = x, v = x;
    double fx = f(x, user);
    double fw = fx, fv = fx;
    double d = 0.0, e = 0.0;        /* step sizes */

    for (int iter = 0; iter < K26C_BRENT_MAX_ITERS; iter++) {
        double xm  = 0.5 * (a + b);
        double tol1 = tol * fabs(x) + K26C_BRENT_TINY;
        double tol2 = 2.0 * tol1;

        if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
            *x_min = x;
            *f_min = fx;
            return K26C_OK;
        }

        int parabolic_ok = 0;
        double u, fu;

        if (fabs(e) > tol1) {
            /* Try parabolic interpolation. */
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            double e_old = e;
            e = d;
            if (fabs(p) < fabs(0.5 * q * e_old) &&
                p > q * (a - x) && p < q * (b - x)) {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) {
                    d = (xm - x >= 0) ? tol1 : -tol1;
                }
                parabolic_ok = 1;
            }
        }
        if (!parabolic_ok) {
            /* Golden-section step. */
            e = (x >= xm) ? a - x : b - x;
            d = K26C_BRENT_GOLDEN_RATIO * e;
        }
        u = (fabs(d) >= tol1) ? x + d : x + ((d >= 0) ? tol1 : -tol1);
        fu = f(u, user);

        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v  = w;  fv = fw;
            w  = x;  fw = fx;
            x  = u;  fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v  = w;  fv = fw;
                w  = u;  fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v  = u;  fv = fu;
            }
        }
    }

    *x_min = x;
    *f_min = fx;
    return K26C_ERR_CONVERGE;
}

/* ---- Nelder-Mead ----
 * Standard reflection / expansion / contraction / shrink simplex
 * with the canonical coefficients (α=1, γ=2, ρ=0.5, σ=0.5). Stops
 * when (max_f - min_f) over the simplex < tol or after max_iters. */

#define NM_ALPHA 1.0
#define NM_GAMMA 2.0
#define NM_RHO   0.5
#define NM_SIGMA 0.5

static double obj_at(K26CObjNFn f, void *user, const double *x, size_t n)
{
    K26CVector v = { (double *)x, n };
    return f(&v, user);
}

K26CStatus k26c_opt_nelder_mead(K26CObjNFn f, void *user,
                                const K26CVector *x0, double simplex_step,
                                double tol, size_t max_iters,
                                K26CVector *x_min, double *f_min)
{
    if (!f || !x0 || !x0->data || x0->n == 0 || !x_min || !f_min)
        return K26C_ERR_INVAL;
    if (tol <= 0.0) tol = 1e-8;
    if (max_iters == 0) max_iters = 1000;
    if (simplex_step == 0.0) simplex_step = 0.05;

    size_t n = x0->n;
    /* Simplex: (n+1) vertices, each n-dim. Stored row-major. */
    double *S    = (double *)malloc((n + 1) * n * sizeof(double));
    double *fS   = (double *)malloc((n + 1)     * sizeof(double));
    double *xc   = (double *)malloc(n * sizeof(double));
    double *xr   = (double *)malloc(n * sizeof(double));
    double *xe   = (double *)malloc(n * sizeof(double));
    double *xcc  = (double *)malloc(n * sizeof(double));
    if (!S || !fS || !xc || !xr || !xe || !xcc) {
        free(S); free(fS); free(xc); free(xr); free(xe); free(xcc);
        return K26C_ERR_OOM;
    }

    /* Build initial simplex around x0. */
    for (size_t i = 0; i < n; i++) S[i] = x0->data[i];
    fS[0] = obj_at(f, user, S, n);
    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 0; j < n; j++) S[i*n + j] = x0->data[j];
        S[i*n + (i - 1)] += simplex_step;
        fS[i] = obj_at(f, user, &S[i*n], n);
    }

    K26CStatus rc = K26C_OK;
    size_t iter = 0;
    for (iter = 0; iter < max_iters; iter++) {
        /* Find best (lo), worst (hi), 2nd-worst (hi2) indices. */
        size_t lo = 0, hi = 0, hi2 = 0;
        for (size_t i = 1; i <= n; i++) {
            if (fS[i] < fS[lo]) lo  = i;
            if (fS[i] > fS[hi]) hi  = i;
        }
        for (size_t i = 0; i <= n; i++) {
            if (i == hi) continue;
            if (fS[i] > fS[hi2] || hi2 == hi) hi2 = i;
        }

        /* Convergence: spread of f-values over the simplex. */
        double spread = fS[hi] - fS[lo];
        if (spread < tol) break;

        /* Centroid of all vertices except `hi`. */
        for (size_t j = 0; j < n; j++) xc[j] = 0.0;
        for (size_t i = 0; i <= n; i++) {
            if (i == hi) continue;
            for (size_t j = 0; j < n; j++) xc[j] += S[i*n + j];
        }
        for (size_t j = 0; j < n; j++) xc[j] /= (double)n;

        /* Reflection. */
        for (size_t j = 0; j < n; j++) xr[j] = xc[j] + NM_ALPHA * (xc[j] - S[hi*n + j]);
        double fr = obj_at(f, user, xr, n);

        if (fr < fS[lo]) {
            /* Expansion. */
            for (size_t j = 0; j < n; j++) xe[j] = xc[j] + NM_GAMMA * (xr[j] - xc[j]);
            double fe = obj_at(f, user, xe, n);
            if (fe < fr) {
                memcpy(&S[hi*n], xe, n * sizeof(double));
                fS[hi] = fe;
            } else {
                memcpy(&S[hi*n], xr, n * sizeof(double));
                fS[hi] = fr;
            }
        } else if (fr < fS[hi2]) {
            /* Accept reflection. */
            memcpy(&S[hi*n], xr, n * sizeof(double));
            fS[hi] = fr;
        } else {
            /* Contraction. */
            int outside = (fr < fS[hi]);
            const double *xref = outside ? xr : &S[hi*n];
            double f_ref       = outside ? fr : fS[hi];
            for (size_t j = 0; j < n; j++) xcc[j] = xc[j] + NM_RHO * (xref[j] - xc[j]);
            double fcc = obj_at(f, user, xcc, n);
            if (fcc < f_ref) {
                memcpy(&S[hi*n], xcc, n * sizeof(double));
                fS[hi] = fcc;
            } else {
                /* Shrink toward best. */
                for (size_t i = 0; i <= n; i++) {
                    if (i == lo) continue;
                    for (size_t j = 0; j < n; j++)
                        S[i*n + j] = S[lo*n + j] + NM_SIGMA * (S[i*n + j] - S[lo*n + j]);
                    fS[i] = obj_at(f, user, &S[i*n], n);
                }
            }
        }
    }

    /* Locate best vertex; emit. */
    size_t lo = 0;
    for (size_t i = 1; i <= n; i++) if (fS[i] < fS[lo]) lo = i;

    K26CVector out;
    K26CStatus alloc_rc = k26c_vec_alloc(&out, n);
    if (alloc_rc != K26C_OK) {
        free(S); free(fS); free(xc); free(xr); free(xe); free(xcc);
        return alloc_rc;
    }
    memcpy(out.data, &S[lo*n], n * sizeof(double));
    k26c_vec_free(x_min);
    *x_min = out;
    *f_min = fS[lo];

    if (iter >= max_iters) rc = K26C_ERR_CONVERGE;

    free(S); free(fS); free(xc); free(xr); free(xe); free(xcc);
    return rc;
}
