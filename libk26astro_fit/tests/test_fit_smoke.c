/* test_fit_smoke.c — gate for libk26astro_fit's two surfaces.
 *
 * Exercises:
 *   1. C-direct API: Rosenbrock minimum recovery from offset start
 *   2. KFL-callable opaque API: linear/quadratic/power/exp/gaussian
 *      fits against synthetic data with known params
 *   3. Determinism: two runs from the same initial guess produce
 *      bit-identical final params */

#include "k26astro_fit/fit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures_ = 0;

static void expect_close_rel_(const char *tag, double got, double want, double rel_tol)
{
    const double err = fabs(got - want);
    const double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
    if (err / scale > rel_tol) {
        fprintf(stderr,
                "FAIL %s: got %.12g, want %.12g (|err|/scale %.3g > rel_tol %.3g)\n",
                tag, got, want, err / scale, rel_tol);
        failures_++;
    } else {
        fprintf(stderr, "PASS %s: got %.12g, |err|/scale %.3g <= rel_tol %.3g\n",
                tag, got, err / scale, rel_tol);
    }
}

/* C-direct residual: Rosenbrock function.
 *   f(x, y) = (1 - x)² + 100·(y - x²)²
 * The minimum is at (1, 1) with f = 0.
 * Express as 2 residuals: r0 = 10·(y - x²), r1 = 1 - x.
 * Sum of squares is exactly the Rosenbrock function. */
static int rhs_rosenbrock_(int m, int n, const double *params,
                           double *fvec, void *user)
{
    (void)m; (void)n; (void)user;
    const double x = params[0];
    const double y = params[1];
    fvec[0] = 10.0 * (y - x * x);
    fvec[1] = 1.0 - x;
    return 0;
}

static void test_c_direct_rosenbrock_(void)
{
    /* Start from (-1.2, 1) — Rosenbrock's classic test starting point. */
    double x[2] = {-1.2, 1.0};
    int rc = k26astro_fit_lmdif(rhs_rosenbrock_, NULL,
                                 2, 2,    /* m=2 residuals, n=2 params */
                                 x, NULL, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL c_direct_rosenbrock: rc=%d\n", rc);
        failures_++;
        return;
    }
    expect_close_rel_("c_direct_rosenbrock x", x[0], 1.0, 1.0e-6);
    expect_close_rel_("c_direct_rosenbrock y", x[1], 1.0, 1.0e-6);
}

static void test_kfl_opaque_linear_(void)
{
    /* Synthetic data: y = 2 + 3*x at x = 0..9. Recover [2, 3]. */
    const int m = 10;
    double xs[10], ys[10];
    for (int i = 0; i < m; ++i) {
        xs[i] = (double)i;
        ys[i] = 2.0 + 3.0 * (double)i;
    }
    K26AstroFitProblem *p = k26astro_fit_problem_linear(m, xs, ys);
    if (!p) {
        fprintf(stderr, "FAIL kfl_linear: alloc\n");
        failures_++;
        return;
    }
    double params[2] = {0.0, 0.0};
    int rc = k26astro_fit_lmdif_h(p, params, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL kfl_linear: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_linear c0", params[0], 2.0, 1.0e-6);
        expect_close_rel_("kfl_linear c1", params[1], 3.0, 1.0e-6);
    }
    k26astro_fit_problem_destroy(p);
}

static void test_kfl_opaque_quadratic_(void)
{
    /* Synthetic data: y = 1 + 2*x + 3*x² at x = 0..9. Recover [1, 2, 3]. */
    const int m = 10;
    double xs[10], ys[10];
    for (int i = 0; i < m; ++i) {
        xs[i] = (double)i;
        ys[i] = 1.0 + 2.0 * xs[i] + 3.0 * xs[i] * xs[i];
    }
    K26AstroFitProblem *p = k26astro_fit_problem_quadratic(m, xs, ys);
    if (!p) {
        fprintf(stderr, "FAIL kfl_quadratic: alloc\n");
        failures_++;
        return;
    }
    double params[3] = {0.0, 0.0, 0.0};
    int rc = k26astro_fit_lmdif_h(p, params, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL kfl_quadratic: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_quadratic c0", params[0], 1.0, 1.0e-6);
        expect_close_rel_("kfl_quadratic c1", params[1], 2.0, 1.0e-6);
        expect_close_rel_("kfl_quadratic c2", params[2], 3.0, 1.0e-6);
    }
    k26astro_fit_problem_destroy(p);
}

static void test_kfl_opaque_power_(void)
{
    /* Synthetic data: y = 2.5 · x^1.7 at x = 1..10. Recover [2.5, 1.7]. */
    const int m = 10;
    double xs[10], ys[10];
    for (int i = 0; i < m; ++i) {
        xs[i] = (double)(i + 1);
        ys[i] = 2.5 * pow(xs[i], 1.7);
    }
    K26AstroFitProblem *p = k26astro_fit_problem_power(m, xs, ys);
    if (!p) {
        fprintf(stderr, "FAIL kfl_power: alloc\n");
        failures_++;
        return;
    }
    double params[2] = {1.0, 1.0};
    int rc = k26astro_fit_lmdif_h(p, params, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL kfl_power: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_power c", params[0], 2.5, 1.0e-5);
        expect_close_rel_("kfl_power n", params[1], 1.7, 1.0e-5);
    }
    k26astro_fit_problem_destroy(p);
}

static void test_kfl_opaque_exp_(void)
{
    /* Synthetic data: y = 5.0 · exp(-0.3·x) at x = 0..9. Recover [5, 0.3]. */
    const int m = 10;
    double xs[10], ys[10];
    for (int i = 0; i < m; ++i) {
        xs[i] = (double)i;
        ys[i] = 5.0 * exp(-0.3 * xs[i]);
    }
    K26AstroFitProblem *p = k26astro_fit_problem_exp(m, xs, ys);
    if (!p) {
        fprintf(stderr, "FAIL kfl_exp: alloc\n");
        failures_++;
        return;
    }
    double params[2] = {1.0, 0.1};
    int rc = k26astro_fit_lmdif_h(p, params, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL kfl_exp: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_exp c",     params[0], 5.0, 1.0e-5);
        expect_close_rel_("kfl_exp alpha", params[1], 0.3, 1.0e-5);
    }
    k26astro_fit_problem_destroy(p);
}

static void test_kfl_opaque_gaussian_(void)
{
    /* Synthetic data: y = 2.0 · exp(-(x-5)²/(2·1.5²)) at x = 0..10.
     * Recover [μ=5, σ=1.5, amp=2.0]. */
    const int m = 11;
    double xs[11], ys[11];
    const double mu_true = 5.0;
    const double sigma_true = 1.5;
    const double amp_true = 2.0;
    for (int i = 0; i < m; ++i) {
        xs[i] = (double)i;
        const double dx = xs[i] - mu_true;
        ys[i] = amp_true * exp(-(dx * dx) / (2.0 * sigma_true * sigma_true));
    }
    K26AstroFitProblem *p = k26astro_fit_problem_gaussian(m, xs, ys);
    if (!p) {
        fprintf(stderr, "FAIL kfl_gaussian: alloc\n");
        failures_++;
        return;
    }
    double params[3] = {4.0, 1.0, 1.0};  /* offset start */
    int rc = k26astro_fit_lmdif_h(p, params, 1.0e-10);
    if (rc != K26ASTRO_FIT_OK) {
        fprintf(stderr, "FAIL kfl_gaussian: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_gaussian mu",    params[0], mu_true,    1.0e-5);
        expect_close_rel_("kfl_gaussian sigma", params[1], sigma_true, 1.0e-5);
        expect_close_rel_("kfl_gaussian amp",   params[2], amp_true,   1.0e-5);
    }
    k26astro_fit_problem_destroy(p);
}

static void test_constructors_alloc_(void)
{
    double xs[4] = {0.0, 1.0, 2.0, 3.0};
    double ys[4] = {0.0, 1.0, 4.0, 9.0};
    K26AstroFitProblem *handles[5] = {
        k26astro_fit_problem_linear   (4, xs, ys),
        k26astro_fit_problem_quadratic(4, xs, ys),
        k26astro_fit_problem_power    (4, xs, ys),
        k26astro_fit_problem_exp      (4, xs, ys),
        k26astro_fit_problem_gaussian (4, xs, ys),
    };
    int ok = 1;
    for (int i = 0; i < 5; ++i) {
        if (!handles[i]) {
            fprintf(stderr, "FAIL constructors_alloc: handle[%d] NULL\n", i);
            failures_++;
            ok = 0;
        }
    }
    if (ok) {
        fprintf(stderr, "PASS constructors_alloc: all 5 handles non-NULL\n");
    }
    for (int i = 0; i < 5; ++i) k26astro_fit_problem_destroy(handles[i]);
    k26astro_fit_problem_destroy(NULL);
}

static void test_determinism_(void)
{
    double x1[2] = {-1.2, 1.0};
    double x2[2] = {-1.2, 1.0};
    k26astro_fit_lmdif(rhs_rosenbrock_, NULL, 2, 2, x1, NULL, 1.0e-10);
    k26astro_fit_lmdif(rhs_rosenbrock_, NULL, 2, 2, x2, NULL, 1.0e-10);
    if (x1[0] != x2[0] || x1[1] != x2[1]) {
        fprintf(stderr,
                "FAIL determinism: run1 [%.17g, %.17g] != run2 [%.17g, %.17g]\n",
                x1[0], x1[1], x2[0], x2[1]);
        failures_++;
    } else {
        fprintf(stderr, "PASS determinism: byte-identical across two runs\n");
    }
}

int main(void)
{
    fprintf(stderr, "libk26astro_fit smoke test\n");

    test_c_direct_rosenbrock_();
    test_kfl_opaque_linear_();
    test_kfl_opaque_quadratic_();
    test_kfl_opaque_power_();
    test_kfl_opaque_exp_();
    test_kfl_opaque_gaussian_();
    test_constructors_alloc_();
    test_determinism_();

    if (failures_) {
        fprintf(stderr, "TOTAL: %d FAIL\n", failures_);
        return 1;
    }
    fprintf(stderr, "TOTAL: all PASS\n");
    return 0;
}
