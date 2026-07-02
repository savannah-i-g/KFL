/* k26astro_ephem/cheby.h — Chebyshev polynomial evaluator.
 *
 * The SPK Type 2 format stores planetary positions as Chebyshev
 * polynomial approximations of low order (typically 14) per fixed-
 * duration block. Evaluating a position at an arbitrary epoch is:
 *
 *   1. Find the block whose [t_start, t_end] contains the epoch.
 *   2. Normalise the epoch into s ∈ [-1, 1]:
 *        s = 2 * (t - t_mid) / t_dur
 *   3. Evaluate Σ c_k * T_k(s) via Clenshaw recurrence.
 *
 * Velocity is the analytic derivative of the same series:
 *   T_k'(s) = k * U_{k-1}(s)
 * evaluated by a parallel Clenshaw-style recurrence on the Chebyshev
 * polynomials of the second kind U_k.
 *
 * Both routines are pure functions; no allocation, no global state.
 */
#ifndef K26ASTRO_EPHEM_CHEBY_H
#define K26ASTRO_EPHEM_CHEBY_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Evaluate Σ_{k=0}^{n-1} coeffs[k] * T_k(s) via Clenshaw recurrence.
 * Numerically stable; identical bit-for-bit across binary64 platforms.
 *
 * Arguments:
 *   coeffs : array of n coefficients (T_0 through T_{n-1})
 *   n      : number of coefficients (≥ 1)
 *   s      : evaluation point in [-1, 1] (no range check)
 *
 * Returns the polynomial value. For n = 1 returns coeffs[0]. */
double k26astro_cheby_eval(const double *coeffs, int n, double s);

/* Evaluate the derivative of the same series at point s, returning
 * the value with respect to s. To get derivative w.r.t. real time
 * (e.g. position → velocity), multiply by ds/dt = 2 / t_dur where
 * t_dur is the block's time span. */
double k26astro_cheby_eval_deriv(const double *coeffs, int n, double s);

/* Evaluate both value and derivative in a single pass. Slightly
 * cheaper than two separate calls (shared Clenshaw pre-roll). */
void   k26astro_cheby_eval_both(const double *coeffs, int n, double s,
                                double *out_val, double *out_deriv);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_EPHEM_CHEBY_H */
