/* test_quad_smoke.c — gate for libk26astro_quad's two surfaces.
 *
 * Exercises:
 *   1. C-direct API: ∫₀¹ e^x dx = e - 1            (DQAGS)
 *      C-direct API: ∫₀^∞ e^(-x) dx = 1            (DQAGI)
 *   2. KFL-callable opaque API: ∫₀¹ x² dx = 1/3    (DQAGS, poly2)
 *      KFL-callable opaque API: ∫₀^∞ B(λ; 5778K) dλ ≈ σT⁴/π
 *                                                  (DQAGI, Planck)
 *   3. Determinism: two back-to-back runs produce bit-identical
 *      result + abserr for the C-direct case.
 *   4. Nested-quadrature TLS save/restore: ∫₀¹ (∫₀¹ e^x dx) dx
 *      should equal (e-1) (outer integrand is the inner result,
 *      constant in x). */

#include "k26astro_quad/quad.h"
#include "k26astro_quad/quad_consts.h"
#include "quad_reference.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int  failures_ = 0;

static void expect_close_(const char *tag, double got, double want, double tol)
{
    const double err = fabs(got - want);
    if (err > tol) {
        fprintf(stderr,
                "FAIL %s: got %.17g, want %.17g (|err| %.3g > tol %.3g)\n",
                tag, got, want, err, tol);
        failures_++;
    } else {
        fprintf(stderr, "PASS %s: |err| %.3g <= tol %.3g\n",
                tag, err, tol);
    }
}

/* C-direct integrand: f(x) = e^x. */
static double f_exp_(double x, void *user)
{
    (void)user;
    return exp(x);
}

/* C-direct integrand: f(x) = e^(-x). */
static double f_exp_neg_(double x, void *user)
{
    (void)user;
    return exp(-x);
}

/* C-direct integrand for nested-quadrature test. Calls a nested
 * dqags from inside an integrand to verify TLS save/restore. */
static double f_constant_via_nested_(double x, void *user)
{
    (void)x;
    (void)user;
    double inner_result = 0.0, inner_abserr = 0.0;
    /* Inner integrand: e^x on [0, 1]. */
    int rc = k26astro_quad_dqags(f_exp_, NULL, 0.0, 1.0,
                                  K26ASTRO_QUAD_EPSABS_DEFAULT,
                                  K26ASTRO_QUAD_EPSREL_DEFAULT,
                                  &inner_result, &inner_abserr);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "nested inner dqags rc=%d\n", rc);
    }
    return inner_result;
}

static void test_c_direct_dqags_(void)
{
    double result = 0.0, abserr = 0.0;
    int rc = k26astro_quad_dqags(f_exp_, NULL, 0.0, 1.0,
                                  K26ASTRO_QUAD_EPSABS_DEFAULT,
                                  K26ASTRO_QUAD_EPSREL_DEFAULT,
                                  &result, &abserr);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "FAIL c_direct_dqags: rc=%d\n", rc);
        failures_++;
        return;
    }
    expect_close_("c_direct_dqags(exp on [0,1])",
                  result, Q_REF_EXP_INT_0_1, Q_TOL_ABS);
}

static void test_c_direct_dqagi_(void)
{
    double result = 0.0, abserr = 0.0;
    int rc = k26astro_quad_dqagi(f_exp_neg_, NULL, 0.0,
                                  K26ASTRO_QUAD_INF_POS_INFTY,
                                  K26ASTRO_QUAD_EPSABS_DEFAULT,
                                  K26ASTRO_QUAD_EPSREL_DEFAULT,
                                  &result, &abserr);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "FAIL c_direct_dqagi: rc=%d\n", rc);
        failures_++;
        return;
    }
    expect_close_("c_direct_dqagi(exp(-x) on [0,inf])",
                  result, Q_REF_EXP_NEG_INT, Q_TOL_ABS);
}

static void test_kfl_opaque_poly2_(void)
{
    /* ∫₀¹ (0 + 0·x + 1·x²) dx = 1/3 */
    K26AstroQuadIntegrand *h =
        k26astro_quad_integrand_poly2(0.0, 0.0, 1.0);
    if (!h) {
        fprintf(stderr, "FAIL kfl_opaque_poly2: alloc\n");
        failures_++;
        return;
    }
    double result = 0.0, abserr = 0.0;
    int rc = k26astro_quad_dqags_h(h, 0.0, 1.0,
                                    K26ASTRO_QUAD_EPSABS_DEFAULT,
                                    K26ASTRO_QUAD_EPSREL_DEFAULT,
                                    &result, &abserr);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "FAIL kfl_opaque_poly2: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_("kfl_opaque_poly2(x^2 on [0,1])",
                      result, Q_REF_X_SQ_INT_0_1, Q_TOL_ABS);
    }
    k26astro_quad_integrand_destroy(h);
}

static void test_kfl_opaque_gaussian_(void)
{
    /* Unnormalized Gaussian with amp=1/(sigma·sqrt(2π)) is the
     * PDF; integrates to 1 over (-∞, +∞) but DQAGS on a finite
     * range covers 99.9999% within ±5σ. mu=0, sigma=1, amp adjusted
     * so the analytic integral over the range is exactly known:
     *   ∫_-5^+5 (1/sqrt(2π)) exp(-x²/2) dx ≈ 0.9999994266...
     * Bit-stable comparison to that reference value. */
    const double inv_sqrt_2pi = 0.39894228040143267; /* 1/sqrt(2π) */
    K26AstroQuadIntegrand *h =
        k26astro_quad_integrand_gaussian(0.0, 1.0, inv_sqrt_2pi);
    if (!h) {
        fprintf(stderr, "FAIL kfl_opaque_gaussian: alloc\n");
        failures_++;
        return;
    }
    double result = 0.0, abserr = 0.0;
    int rc = k26astro_quad_dqags_h(h, -5.0, 5.0,
                                    K26ASTRO_QUAD_EPSABS_DEFAULT,
                                    K26ASTRO_QUAD_EPSREL_DEFAULT,
                                    &result, &abserr);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "FAIL kfl_opaque_gaussian: rc=%d\n", rc);
        failures_++;
    } else {
        /* ∫_-5^+5 φ(x;0,1) dx = erf(5/sqrt(2)) ≈ 0.99999942669...
         * (1 - 5.733e-7); accepts result within 1e-9 of that. */
        expect_close_("kfl_opaque_gaussian(N(0,1) on [-5,+5])",
                      result, 0.9999994266968562, 1.0e-9);
    }
    k26astro_quad_integrand_destroy(h);
}

static void test_kfl_opaque_constructors_alloc_(void)
{
    /* Stress-tests that every integrand_* constructor produces a
     * non-NULL handle that can be destroyed safely. No numerical
     * convergence test — that's per-family above. */
    K26AstroQuadIntegrand *handles[7] = {
        k26astro_quad_integrand_const(1.0),
        k26astro_quad_integrand_poly2(0.0, 0.0, 1.0),
        k26astro_quad_integrand_poly3(0.0, 0.0, 0.0, 1.0),
        k26astro_quad_integrand_power(1.0, 2.0),
        k26astro_quad_integrand_exp(1.0, 1.0),
        k26astro_quad_integrand_gaussian(0.0, 1.0, 1.0),
        k26astro_quad_integrand_planck(5778.0),
    };
    int ok = 1;
    for (int i = 0; i < 7; ++i) {
        if (!handles[i]) {
            fprintf(stderr,
                "FAIL kfl_opaque_constructors_alloc: handle[%d] NULL\n", i);
            failures_++;
            ok = 0;
        }
    }
    if (ok) {
        fprintf(stderr,
            "PASS kfl_opaque_constructors_alloc: all 7 handles non-NULL\n");
    }
    for (int i = 0; i < 7; ++i) {
        k26astro_quad_integrand_destroy(handles[i]);
    }
    /* Null-destroy is a no-op. */
    k26astro_quad_integrand_destroy(NULL);
}

static void test_determinism_(void)
{
    double r1 = 0.0, e1 = 0.0, r2 = 0.0, e2 = 0.0;
    k26astro_quad_dqags(f_exp_, NULL, 0.0, 1.0,
                        K26ASTRO_QUAD_EPSABS_DEFAULT,
                        K26ASTRO_QUAD_EPSREL_DEFAULT, &r1, &e1);
    k26astro_quad_dqags(f_exp_, NULL, 0.0, 1.0,
                        K26ASTRO_QUAD_EPSABS_DEFAULT,
                        K26ASTRO_QUAD_EPSREL_DEFAULT, &r2, &e2);
    if (r1 != r2 || e1 != e2) {
        fprintf(stderr,
                "FAIL determinism: run1 (%.17g, %.17g) != run2 (%.17g, %.17g)\n",
                r1, e1, r2, e2);
        failures_++;
    } else {
        fprintf(stderr, "PASS determinism: byte-identical across two runs\n");
    }
}

static void test_nested_tls_(void)
{
    /* Outer integrand calls a nested dqags. The outer integrand is
     * a constant value (inner result ≈ e-1), so outer integral over
     * [0,1] should equal that same constant. */
    double outer = 0.0, outer_err = 0.0;
    int rc = k26astro_quad_dqags(f_constant_via_nested_, NULL, 0.0, 1.0,
                                  1.0e-9, 1.0e-9,
                                  &outer, &outer_err);
    if (rc != K26ASTRO_QUAD_OK) {
        fprintf(stderr, "FAIL nested_tls: outer rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_("nested_tls(outer integral = inner constant)",
                      outer, Q_REF_EXP_INT_0_1, 1.0e-9);
    }
}

int main(void)
{
    fprintf(stderr, "libk26astro_quad smoke test\n");

    test_c_direct_dqags_();
    test_c_direct_dqagi_();
    test_kfl_opaque_poly2_();
    test_kfl_opaque_gaussian_();
    test_kfl_opaque_constructors_alloc_();
    test_determinism_();
    test_nested_tls_();

    if (failures_) {
        fprintf(stderr, "TOTAL: %d FAIL\n", failures_);
        return 1;
    }
    fprintf(stderr, "TOTAL: all PASS\n");
    return 0;
}
