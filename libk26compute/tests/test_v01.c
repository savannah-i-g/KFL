/* libk26compute smoke test against known-answer problems.
 *
 * Each test prints PASS / FAIL and the value / error margin. Exits 0 if
 * everything passes, 1 otherwise. Tolerances are deliberately loose
 * (1e-6 to 1e-3 depending on the operation); these are smoke tests,
 * not bit-exact reproducibility checks. The reproducibility contract
 * is exercised separately by repeated calls with same seed. */

#include "k26compute.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failures = 0;

static void check(const char *name, int ok, const char *detail)
{
    if (ok) {
        printf("PASS  %s  (%s)\n", name, detail);
    } else {
        printf("FAIL  %s  (%s)\n", name, detail);
        g_failures++;
    }
}

/* ---- 1. Linear regression: Hooke's law F = k*x with k=2.5.
 *       Synthetic data with no noise → b ≈ 2.5, a ≈ 0, r² ≈ 1.0. */
static void test_linreg(void)
{
    double xs[] = { 0, 1, 2, 3, 4, 5 };
    double ys[6];
    for (int i = 0; i < 6; i++) ys[i] = 2.5 * xs[i];

    K26CLinReg fit;
    K26CStatus rc = k26c_stats_linreg(&fit, xs, ys, 6);
    char detail[160];
    snprintf(detail, sizeof detail,
             "rc=%s a=%.6f b=%.6f r2=%.6f",
             k26compute_status_str(rc), fit.a, fit.b, fit.r2);
    int ok = (rc == K26C_OK
              && fabs(fit.a - 0.0) < 1e-9
              && fabs(fit.b - 2.5) < 1e-9
              && fit.r2 > 0.999999);
    check("linreg: Hooke's law", ok, detail);
}

/* ---- 2. Matrix inverse + multiply: A * A^-1 == I (within tol). */
static void test_mat_inverse(void)
{
    double A_init[] = {  4, 7,
                          2, 6 };          /* det = 4*6 - 7*2 = 10 */
    K26CMatrix A, Ainv, P;
    k26c_mat_from(&A, A_init, 2, 2);
    k26c_mat_alloc(&Ainv, 0, 0);
    k26c_mat_alloc(&P,    0, 0);

    K26CStatus rc = k26c_mat_inverse(&Ainv, &A);
    int ok = (rc == K26C_OK);
    if (ok) {
        rc = k26c_mat_mul(&P, &A, &Ainv);
        if (rc == K26C_OK) {
            /* Compare P to identity. */
            double err = 0.0;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++) {
                    double exp = (i == j) ? 1.0 : 0.0;
                    err += fabs(k26c_mat_get(&P, i, j) - exp);
                }
            char detail[120];
            snprintf(detail, sizeof detail,
                     "rc=%s ||A·A^-1 - I||_1 = %.3e",
                     k26compute_status_str(rc), err);
            check("mat_inverse: 2x2 round-trip", err < 1e-12, detail);
        }
    }
    k26c_mat_free(&A);  k26c_mat_free(&Ainv);  k26c_mat_free(&P);
}

/* ---- 3. Determinant: 2x2 known.  det([[4,7],[2,6]]) = 10. */
static void test_mat_det(void)
{
    double A_init[] = { 4, 7,
                         2, 6 };
    K26CMatrix A;
    k26c_mat_from(&A, A_init, 2, 2);
    double d;
    K26CStatus rc = k26c_mat_det(&d, &A);
    char detail[80];
    snprintf(detail, sizeof detail,
             "rc=%s det=%.10g (expected 10)",
             k26compute_status_str(rc), d);
    check("mat_det: 2x2", rc == K26C_OK && fabs(d - 10.0) < 1e-12, detail);
    k26c_mat_free(&A);
}

/* ---- 4. ODE: y' = -y, y(0) = 1 → y(1) = 1/e ≈ 0.36788. */
static int rhs_decay(double t, const K26CVector *y, K26CVector *dy, void *user)
{
    (void)t; (void)user;
    dy->data[0] = -y->data[0];
    return 0;
}
static void test_ode_rk45(void)
{
    K26CVector y;
    k26c_vec_alloc(&y, 1);
    y.data[0] = 1.0;
    K26CStatus rc = k26c_ode_rk45(rhs_decay, NULL, 0.0, 1.0, 1e-9, 1e-12, &y);
    double expected = 1.0 / M_E;
    char detail[120];
    snprintf(detail, sizeof detail,
             "rc=%s y(1)=%.12f  expected 1/e=%.12f  err=%.3e",
             k26compute_status_str(rc), y.data[0], expected,
             fabs(y.data[0] - expected));
    check("ode_rk45: y'=-y", rc == K26C_OK && fabs(y.data[0] - expected) < 1e-7, detail);
    k26c_vec_free(&y);
}

/* ---- 5. Brent: minimum of (x-3)^2 + 1 is at x=3, f=1. */
static double f_parabola(double x, void *user)
{
    (void)user;
    return (x - 3.0) * (x - 3.0) + 1.0;
}
static void test_brent(void)
{
    double xm, fm;
    K26CStatus rc = k26c_opt_brent(f_parabola, NULL, -10, 10, 1e-9, &xm, &fm);
    char detail[120];
    snprintf(detail, sizeof detail,
             "rc=%s x=%.10f f=%.10f", k26compute_status_str(rc), xm, fm);
    check("opt_brent: (x-3)^2+1",
          rc == K26C_OK && fabs(xm - 3.0) < 1e-6 && fabs(fm - 1.0) < 1e-12,
          detail);
}

/* ---- 6. Nelder-Mead: Rosenbrock function — minimum at (1, 1) = 0. */
static double f_rosenbrock(const K26CVector *x, void *user)
{
    (void)user;
    double a = 1.0 - x->data[0];
    double b = x->data[1] - x->data[0] * x->data[0];
    return a*a + 100.0 * b * b;
}
static void test_nelder_mead(void)
{
    K26CVector x0;
    k26c_vec_from(&x0, (double[]){-1.2, 1.0}, 2);
    K26CVector xm;
    k26c_vec_alloc(&xm, 0);
    double fm;
    K26CStatus rc = k26c_opt_nelder_mead(f_rosenbrock, NULL, &x0, 0.5, 1e-10,
                                          5000, &xm, &fm);
    char detail[160];
    snprintf(detail, sizeof detail,
             "rc=%s x=(%.6f,%.6f) f=%.3e",
             k26compute_status_str(rc),
             xm.n >= 1 ? xm.data[0] : 0.0,
             xm.n >= 2 ? xm.data[1] : 0.0, fm);
    int ok = (rc == K26C_OK
              && fabs(xm.data[0] - 1.0) < 1e-3
              && fabs(xm.data[1] - 1.0) < 1e-3
              && fm < 1e-6);
    check("opt_nelder_mead: Rosenbrock", ok, detail);
    k26c_vec_free(&x0); k26c_vec_free(&xm);
}

/* ---- 7. RNG determinism: same seed → same first three uniforms. */
static void test_rng_determinism(void)
{
    K26CRng a, b;
    k26c_rng_init(&a, 42);
    k26c_rng_init(&b, 42);
    double a1 = k26c_rng_uniform(&a), a2 = k26c_rng_uniform(&a), a3 = k26c_rng_uniform(&a);
    double b1 = k26c_rng_uniform(&b), b2 = k26c_rng_uniform(&b), b3 = k26c_rng_uniform(&b);
    int ok = (a1 == b1 && a2 == b2 && a3 == b3);
    char detail[160];
    snprintf(detail, sizeof detail,
             "(%.6f,%.6f,%.6f) vs (%.6f,%.6f,%.6f)",
             a1, a2, a3, b1, b2, b3);
    check("rng: same seed → same outputs", ok, detail);
}

int main(void)
{
    printf("libk26compute v%s smoke test\n\n", K26COMPUTE_LIB_VERSION);

    test_linreg();
    test_mat_inverse();
    test_mat_det();
    test_ode_rk45();
    test_brent();
    test_nelder_mead();
    test_rng_determinism();

    printf("\n%s — %d failure(s)\n", g_failures == 0 ? "OK" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
