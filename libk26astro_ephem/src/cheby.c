/* cheby.c — Clenshaw recurrence for Chebyshev T_k(s).
 *
 * Reference: Numerical Recipes 3rd ed. §5.8 — the Clenshaw form is
 * numerically stable across |s| ≤ 1 and produces bit-for-bit
 * identical output regardless of compiler / FPU mode (in our
 * `-ffp-contract=off`, `-fexcess-precision=standard` configuration).
 *
 * For position evaluation: T_k(s) recurrence, k = 0..n-1.
 *
 *   b_{n}   = 0
 *   b_{n-1} = 0
 *   b_{k}   = 2 s b_{k+1} - b_{k+2} + c_k     for k = n-2 down to 1
 *   result  = s b_1 - b_2 + c_0
 *
 * For derivative w.r.t. s: T_k'(s) = k U_{k-1}(s), where U_n satisfies
 * the same recurrence (different boundary conditions).
 *
 *   d_{n}   = 0
 *   d_{n-1} = 0
 *   d_{k}   = 2 s d_{k+1} - d_{k+2} + (k+1) c_{k+1}    for k = n-2 down to 1
 *   result  = s d_1 - d_2 + c_1
 *
 * (The "(k+1) c_{k+1}" injection accounts for k T_k = T_k', which in
 * U-form becomes k U_{k-1}. For k=0, c_0 has no derivative term so
 * its derivative contribution vanishes.) */
#include "k26astro_ephem/cheby.h"

double k26astro_cheby_eval(const double *coeffs, int n, double s)
{
    if (n <= 0) return 0.0;
    if (n == 1) return coeffs[0];
    double b_kp2 = 0.0;
    double b_kp1 = 0.0;
    double b_k   = 0.0;
    /* Walk k = n-1 down to 1; at each step compute b_k from b_kp1,
     * b_kp2, and c_k. */
    for (int k = n - 1; k >= 1; k--) {
        b_k   = 2.0 * s * b_kp1 - b_kp2 + coeffs[k];
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    /* Final fold-in of c_0. */
    return s * b_kp1 - b_kp2 + coeffs[0];
}

double k26astro_cheby_eval_deriv(const double *coeffs, int n, double s)
{
    if (n <= 1) return 0.0;
    if (n == 2) return coeffs[1];
    /* d/ds Σ c_k T_k(s) = Σ_{k≥1} c_k k U_{k-1}(s), which we
     * rewrite as Σ_{j=0..n-2} a_j U_j(s) with a_j = (j+1) c_{j+1}.
     * For a Chebyshev-U series, Clenshaw uses the boundary
     *   b_{N+1} = b_{N+2} = 0
     *   b_k     = 2 s b_{k+1} - b_{k+2} + a_k    (k = N..0)
     *   result  = b_0
     * since U_0 = 1 and the φ_1 - α φ_0 term vanishes. */
    double b_kp2 = 0.0;
    double b_kp1 = 0.0;
    double b_k   = 0.0;
    for (int k = n - 2; k >= 0; k--) {
        double a_k = (double)(k + 1) * coeffs[k + 1];
        b_k   = 2.0 * s * b_kp1 - b_kp2 + a_k;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    return b_k;
}

void k26astro_cheby_eval_both(const double *coeffs, int n, double s,
                               double *out_val, double *out_deriv)
{
    if (out_val)   *out_val   = k26astro_cheby_eval      (coeffs, n, s);
    if (out_deriv) *out_deriv = k26astro_cheby_eval_deriv(coeffs, n, s);
}
