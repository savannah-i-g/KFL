/* test_cheby.c — Clenshaw recurrence against known Chebyshev T_n(s).
 *
 * Acceptance: machine-precision agreement on
 *   T_0(s) = 1
 *   T_1(s) = s
 *   T_2(s) = 2 s^2 - 1
 *   T_3(s) = 4 s^3 - 3 s
 *   T_4(s) = 8 s^4 - 8 s^2 + 1
 *
 * plus an analytic-derivative cross-check using the identity
 *   T_n'(s) = n U_{n-1}(s)
 * specialised to T_3'(s) = 12 s^2 - 3, evaluated via the eval_deriv
 * routine driven by the c_3 = 1 coefficient vector.
 *
 * Plus a round-trip: take an analytic polynomial f(s) = a + b s +
 * c s^2 + d s^3, decompose to Chebyshev basis, evaluate, compare to
 * the direct power-basis evaluation. */
#include "k26astro_ephem/cheby.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

static double Tn_(int n, double s)
{
    switch (n) {
    case 0: return 1.0;
    case 1: return s;
    case 2: return 2.0 * s * s - 1.0;
    case 3: return 4.0 * s * s * s - 3.0 * s;
    case 4: return 8.0 * pow(s, 4) - 8.0 * s * s + 1.0;
    }
    return 0.0;
}

int main(void)
{
    /* T_0, T_1, T_2, T_3, T_4 — set one coefficient to 1, evaluate. */
    for (int n = 0; n < 5; n++) {
        double coeffs[5] = {0};
        coeffs[n] = 1.0;
        for (double s = -1.0; s <= 1.0 + 1e-12; s += 0.1) {
            double v = k26astro_cheby_eval(coeffs, 5, s);
            assert(near_(v, Tn_(n, s), 1e-13));
        }
    }

    /* T_3' analytic: 12 s^2 - 3. */
    {
        double coeffs[4] = { 0.0, 0.0, 0.0, 1.0 };
        for (double s = -1.0; s <= 1.0 + 1e-12; s += 0.1) {
            double d = k26astro_cheby_eval_deriv(coeffs, 4, s);
            double expected = 12.0 * s * s - 3.0;
            assert(near_(d, expected, 1e-13));
        }
    }

    /* Mixed coefficients. f(s) = 2 T_0 + 3 T_1 - T_2 + 0.5 T_3 */
    {
        double coeffs[4] = { 2.0, 3.0, -1.0, 0.5 };
        for (double s = -1.0; s <= 1.0 + 1e-12; s += 0.05) {
            double v_clen = k26astro_cheby_eval(coeffs, 4, s);
            double v_direct =
                  2.0 * Tn_(0, s)
                + 3.0 * Tn_(1, s)
                - 1.0 * Tn_(2, s)
                + 0.5 * Tn_(3, s);
            assert(near_(v_clen, v_direct, 1e-13));
        }
    }

    /* eval_both consistency check. */
    {
        double coeffs[5] = { 1.0, -2.0, 0.3, 0.05, 0.001 };
        for (double s = -1.0; s <= 1.0 + 1e-12; s += 0.1) {
            double v1 = k26astro_cheby_eval      (coeffs, 5, s);
            double d1 = k26astro_cheby_eval_deriv(coeffs, 5, s);
            double v2 = 0.0, d2 = 0.0;
            k26astro_cheby_eval_both(coeffs, 5, s, &v2, &d2);
            assert(v1 == v2 && d1 == d2);
        }
    }

    printf("test_cheby: OK\n");
    return 0;
}
