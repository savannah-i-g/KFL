/* libk26compute — descriptive + inferential statistics. */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double k26c_stats_mean(const double *xs, size_t n)
{
    if (!xs || n == 0) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += xs[i];
    return s / (double)n;
}

double k26c_stats_var(const double *xs, size_t n, int sample)
{
    if (!xs || n == 0) return 0.0;
    if (sample && n < 2) return 0.0;
    double m  = k26c_stats_mean(xs, n);
    double ss = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = xs[i] - m;
        ss += d * d;
    }
    double denom = sample ? (double)(n - 1) : (double)n;
    return ss / denom;
}

double k26c_stats_std(const double *xs, size_t n, int sample)
{
    return sqrt(k26c_stats_var(xs, n, sample));
}

double k26c_stats_min(const double *xs, size_t n)
{
    if (!xs || n == 0) return 0.0;
    double m = xs[0];
    for (size_t i = 1; i < n; i++) if (xs[i] < m) m = xs[i];
    return m;
}

double k26c_stats_max(const double *xs, size_t n)
{
    if (!xs || n == 0) return 0.0;
    double m = xs[0];
    for (size_t i = 1; i < n; i++) if (xs[i] > m) m = xs[i];
    return m;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

double k26c_stats_quantile(const double *xs, size_t n, double q, double *scratch)
{
    if (!xs || !scratch || n == 0) return 0.0;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    memcpy(scratch, xs, n * sizeof(double));
    qsort(scratch, n, sizeof(double), cmp_double);
    /* Linear interpolation between the two nearest order statistics
     * (matches NumPy's default 'linear' interpolation method). */
    double pos = q * (double)(n - 1);
    size_t lo  = (size_t)floor(pos);
    size_t hi  = (size_t)ceil (pos);
    double frac = pos - (double)lo;
    return scratch[lo] + (scratch[hi] - scratch[lo]) * frac;
}

double k26c_stats_median(const double *xs, size_t n, double *scratch)
{
    return k26c_stats_quantile(xs, n, 0.5, scratch);
}

double k26c_stats_pearson(const double *xs, const double *ys, size_t n)
{
    if (!xs || !ys || n < 2) return 0.0;
    double mx = k26c_stats_mean(xs, n);
    double my = k26c_stats_mean(ys, n);
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dx = xs[i] - mx, dy = ys[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    double denom = sqrt(sxx * syy);
    if (denom == 0.0) return 0.0;
    return sxy / denom;
}

K26CStatus k26c_stats_linreg(K26CLinReg *out, const double *xs,
                             const double *ys, size_t n)
{
    if (!out || !xs || !ys || n < 2) return K26C_ERR_INVAL;
    /* Closed-form ordinary least squares. */
    double mx = k26c_stats_mean(xs, n);
    double my = k26c_stats_mean(ys, n);
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dx = xs[i] - mx, dy = ys[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (sxx == 0.0) return K26C_ERR_INVAL;          /* zero-variance x */

    double b = sxy / sxx;
    double a = my  - b * mx;
    /* Residual variance with n-2 dof. Guard against n == 2 (residual
     * variance is 0 by construction). */
    double resid_ss = syy - b * sxy;
    if (resid_ss < 0.0) resid_ss = 0.0;             /* numeric noise */
    double s2  = (n > 2) ? resid_ss / (double)(n - 2) : 0.0;
    double b_se = (sxx > 0.0) ? sqrt(s2 / sxx) : 0.0;
    double a_se = sqrt(s2 * (1.0 / (double)n + mx * mx / sxx));
    double r2   = (syy > 0.0) ? 1.0 - resid_ss / syy : 1.0;

    out->a    = a;
    out->b    = b;
    out->a_se = a_se;
    out->b_se = b_se;
    out->r2   = r2;
    out->n    = n;
    return K26C_OK;
}

double k26c_stats_chi2(const double *observed, const double *expected, size_t n)
{
    if (!observed || !expected || n == 0) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (expected[i] == 0.0) continue;        /* skip zero-expected bins */
        double d = observed[i] - expected[i];
        s += d * d / expected[i];
    }
    return s;
}
