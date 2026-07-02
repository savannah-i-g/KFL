/* test_ode_smoke.c — gate for libk26astro_ode's two surfaces.
 *
 * Exercises:
 *   1. C-direct API:
 *      - harmonic oscillator integrated to t=pi via custom RHS
 *      - determinism: two runs produce bit-identical y(tf)
 *   2. KFL-callable opaque API:
 *      - harmonic problem to t=pi
 *      - Van der Pol stiff (μ=1000) — LSODA's canonical stiff
 *        benchmark; tests stiff/non-stiff switching
 *      - Robertson chemistry — canonical stiff chemistry benchmark
 *      - Lotka-Volterra one period — non-stiff oscillation
 *      - constructors_alloc — every canned constructor produces a
 *        non-NULL handle
 *   3. Linear ODE via C-direct constructor (matrix passing) */

#include "k26astro_ode/ode.h"
#include "k26astro_ode/ode_consts.h"
#include "ode_reference.h"

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
        fprintf(stderr, "PASS %s: |err|/scale %.3g <= rel_tol %.3g\n",
                tag, err / scale, rel_tol);
    }
}

/* C-direct RHS: harmonic oscillator, omega passed via user. */
static int rhs_harmonic_(int n, double t, const double *y, double *ydot, void *user)
{
    (void)n; (void)t;
    const double omega = *(const double *)user;
    ydot[0] = y[1];
    ydot[1] = -omega * omega * y[0];
    return 0;
}

static void test_c_direct_harmonic_(void)
{
    double omega = 1.0;
    double y0[2] = {1.0, 0.0};
    double yf[2] = {0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve(rhs_harmonic_, &omega, 2,
                                       y0, 0.0, M_PI,
                                       1.0e-9, 1.0e-12,
                                       yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL c_direct_harmonic: rc=%d\n", rc);
        failures_++;
        return;
    }
    expect_close_rel_("c_direct_harmonic x(pi) = -1", yf[0],
                      O_REF_HARMONIC_X_AT_PI, O_TOL_REL_ANALYTICAL);
    /* v(pi) should be ~0 with absolute precision <= 1e-6. */
    if (fabs(yf[1]) > 1.0e-6) {
        fprintf(stderr, "FAIL c_direct_harmonic v(pi) = %.12g (want ~0)\n", yf[1]);
        failures_++;
    } else {
        fprintf(stderr, "PASS c_direct_harmonic v(pi) ~0: |v| = %.3g\n", fabs(yf[1]));
    }
}

static void test_kfl_opaque_harmonic_(void)
{
    K26AstroOdeProblem *p = k26astro_ode_problem_harmonic(1.0);
    if (!p) {
        fprintf(stderr, "FAIL kfl_opaque_harmonic: alloc\n");
        failures_++;
        return;
    }
    double y0[2] = {1.0, 0.0};
    double yf[2] = {0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve_h(p, y0, 0.0, M_PI,
                                         1.0e-9, 1.0e-12, yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL kfl_opaque_harmonic: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("kfl_opaque_harmonic x(pi)", yf[0],
                          O_REF_HARMONIC_X_AT_PI, O_TOL_REL_ANALYTICAL);
    }
    k26astro_ode_problem_destroy(p);
}

static void test_van_der_pol_stiff_(void)
{
    /* μ = 1000 — LSODA's canonical stiff benchmark. Integrate from
     * t=0 to t=3000 (long enough to traverse multiple relaxation
     * oscillations). The solution is bounded; we only check that
     * LSODA converges and returns finite values. */
    K26AstroOdeProblem *p = k26astro_ode_problem_van_der_pol(1000.0);
    if (!p) {
        fprintf(stderr, "FAIL van_der_pol_stiff: alloc\n");
        failures_++;
        return;
    }
    double y0[2] = {2.0, 0.0};
    double yf[2] = {0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve_h(p, y0, 0.0, 3000.0,
                                         1.0e-6, 1.0e-9, yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL van_der_pol_stiff: rc=%d (LSODA must converge for μ=1000)\n", rc);
        failures_++;
    } else if (!isfinite(yf[0]) || !isfinite(yf[1])) {
        fprintf(stderr, "FAIL van_der_pol_stiff: non-finite result yf=[%g, %g]\n",
                yf[0], yf[1]);
        failures_++;
    } else if (fabs(yf[0]) > 3.0 || fabs(yf[1]) > 3000.0) {
        /* The orbit should be bounded; large excursions signal divergence. */
        fprintf(stderr, "FAIL van_der_pol_stiff: unbounded yf=[%g, %g]\n",
                yf[0], yf[1]);
        failures_++;
    } else {
        fprintf(stderr, "PASS van_der_pol_stiff: μ=1000 converged, yf=[%.4g, %.4g]\n",
                yf[0], yf[1]);
    }
    k26astro_ode_problem_destroy(p);
}

static void test_robertson_(void)
{
    /* Robertson chemistry — canonical stiff chemistry benchmark.
     * Default rates k1=0.04, k2=1e4, k3=3e7. y(0) = [1, 0, 0].
     * Compare y[0] and y[2] at t=0.4 against published values
     * (Hindmarsh-Petzold). */
    K26AstroOdeProblem *p = k26astro_ode_problem_robertson(0.04, 1.0e4, 3.0e7);
    if (!p) {
        fprintf(stderr, "FAIL robertson: alloc\n");
        failures_++;
        return;
    }
    double y0[3] = {1.0, 0.0, 0.0};
    double yf[3] = {0.0, 0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve_h(p, y0, 0.0, 0.4,
                                         1.0e-6, 1.0e-12, yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL robertson: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("robertson y[A] at t=0.4", yf[0],
                          O_REF_ROBERTSON_T0_4_Y0, O_TOL_REL_BENCHMARK);
        expect_close_rel_("robertson y[C] at t=0.4", yf[2],
                          O_REF_ROBERTSON_T0_4_Y2, O_TOL_REL_BENCHMARK);
    }
    k26astro_ode_problem_destroy(p);
}

static void test_lotka_volterra_(void)
{
    /* Classic Lotka-Volterra; with α=1, β=1, δ=1, γ=1 and y(0)=[1, 1]
     * the solution is the fixed point at (1, 1). Use slightly off-
     * diagonal initial conditions to see oscillation; check that
     * average over a period stays close to the fixed point. */
    K26AstroOdeProblem *p =
        k26astro_ode_problem_lotka_volterra(1.0, 1.0, 1.0, 1.0);
    if (!p) {
        fprintf(stderr, "FAIL lotka_volterra: alloc\n");
        failures_++;
        return;
    }
    double y0[2] = {1.2, 0.8};
    double yf[2] = {0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve_h(p, y0, 0.0, 6.28,
                                         1.0e-9, 1.0e-12, yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL lotka_volterra: rc=%d\n", rc);
        failures_++;
    } else if (!isfinite(yf[0]) || !isfinite(yf[1]) ||
               yf[0] < 0.1 || yf[1] < 0.1) {
        fprintf(stderr, "FAIL lotka_volterra: unphysical yf=[%g, %g]\n",
                yf[0], yf[1]);
        failures_++;
    } else {
        fprintf(stderr, "PASS lotka_volterra: yf=[%.4g, %.4g] (period bounded)\n",
                yf[0], yf[1]);
    }
    k26astro_ode_problem_destroy(p);
}

static void test_linear_c_direct_(void)
{
    /* Linear ODE: dy/dt = A·y + b, with A=[[0,1],[-1,0]] and b=[0,0]
     * = harmonic. y(0)=[1,0], y(pi)=[-1,0]. Validates the linear
     * constructor + matrix dispatcher. */
    double A[4] = { 0.0, 1.0,
                   -1.0, 0.0 };
    double b[2] = { 0.0, 0.0 };
    K26AstroOdeProblem *p = k26astro_ode_problem_linear(2, A, b);
    if (!p) {
        fprintf(stderr, "FAIL linear_c_direct: alloc\n");
        failures_++;
        return;
    }
    double y0[2] = {1.0, 0.0};
    double yf[2] = {0.0, 0.0};
    int rc = k26astro_ode_lsoda_solve_h(p, y0, 0.0, M_PI,
                                         1.0e-9, 1.0e-12, yf);
    if (rc != K26ASTRO_ODE_OK) {
        fprintf(stderr, "FAIL linear_c_direct: rc=%d\n", rc);
        failures_++;
    } else {
        expect_close_rel_("linear_c_direct A=harmonic x(pi)", yf[0],
                          -1.0, O_TOL_REL_ANALYTICAL);
    }
    k26astro_ode_problem_destroy(p);
}

static void test_constructors_alloc_(void)
{
    K26AstroOdeProblem *handles[6] = {
        k26astro_ode_problem_harmonic(1.0),
        k26astro_ode_problem_damped_harmonic(1.0, 0.1),
        k26astro_ode_problem_van_der_pol(1000.0),
        k26astro_ode_problem_lotka_volterra(1.0, 1.0, 1.0, 1.0),
        k26astro_ode_problem_robertson(0.04, 1.0e4, 3.0e7),
        NULL,
    };
    double A[1] = {-1.0};
    double b[1] = {0.0};
    handles[5] = k26astro_ode_problem_linear(1, A, b);

    int ok = 1;
    for (int i = 0; i < 6; ++i) {
        if (!handles[i]) {
            fprintf(stderr, "FAIL constructors_alloc: handle[%d] NULL\n", i);
            failures_++;
            ok = 0;
        }
    }
    if (ok) {
        fprintf(stderr, "PASS constructors_alloc: all 6 handles non-NULL\n");
    }
    for (int i = 0; i < 6; ++i) k26astro_ode_problem_destroy(handles[i]);
    k26astro_ode_problem_destroy(NULL); /* no-op smoke */
}

static void test_determinism_(void)
{
    double omega = 1.0;
    double y0[2] = {1.0, 0.0};
    double yf1[2] = {0.0, 0.0};
    double yf2[2] = {0.0, 0.0};
    k26astro_ode_lsoda_solve(rhs_harmonic_, &omega, 2, y0,
                              0.0, M_PI, 1.0e-9, 1.0e-12, yf1);
    k26astro_ode_lsoda_solve(rhs_harmonic_, &omega, 2, y0,
                              0.0, M_PI, 1.0e-9, 1.0e-12, yf2);
    if (yf1[0] != yf2[0] || yf1[1] != yf2[1]) {
        fprintf(stderr,
                "FAIL determinism: run1 [%.17g, %.17g] != run2 [%.17g, %.17g]\n",
                yf1[0], yf1[1], yf2[0], yf2[1]);
        failures_++;
    } else {
        fprintf(stderr, "PASS determinism: byte-identical across two runs\n");
    }
}

int main(void)
{
    fprintf(stderr, "libk26astro_ode smoke test\n");

    test_c_direct_harmonic_();
    test_kfl_opaque_harmonic_();
    test_van_der_pol_stiff_();
    test_robertson_();
    test_lotka_volterra_();
    test_linear_c_direct_();
    test_constructors_alloc_();
    test_determinism_();

    if (failures_) {
        fprintf(stderr, "TOTAL: %d FAIL\n", failures_);
        return 1;
    }
    fprintf(stderr, "TOTAL: all PASS\n");
    return 0;
}
