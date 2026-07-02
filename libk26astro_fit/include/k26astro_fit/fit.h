/* libk26astro_fit - nonlinear least-squares via CMINPACK (MINPACK in C).
 *
 * Levenberg-Marquardt minimizer for sum-of-squares cost functions.
 * cminpack (devernay; BSD-3) is a pure-C rewrite of MINPACK 1
 * (Argonne 1980); no Fortran involvement. Pure C, joining the
 * numerics suite alongside the Fortran-backed libs
 * (libk26astro_quad, libk26astro_ode).
 *
 * Two surfaces:
 *   1. C-direct API - `k26astro_fit_lmdif(K26AstroFitFn f, void *user,
 *      ...)`. Fully general; user supplies a residual callback.
 *   2. KFL-callable opaque API - `K26AstroFitProblem` handle
 *      constructed via canned residual families (linear regression,
 *      power-law fit, gaussian fit, exponential decay fit). Arbitrary
 *      user-defined KFL residual callbacks require a kflc `fn_ref`
 *      arg-type that is not yet provided; the pre-baked families
 *      cover the current consumer set (orbit determination linear
 *      systems, drag coefficient power-law fits, spectroscopic
 *      gaussian profiles).
 *
 * User-data threading: K26AstroFitFn takes an opaque `void *user`.
 * cminpack's callback signature already includes a user-data pointer,
 * so no TLS is required (unlike the Fortran-backed wrappers in
 * libk26astro_quad / libk26astro_ode). The K26 → cminpack callback
 * adapter is a thin struct unwrap.
 *
 * Error model: returns the cminpack info code at exit (0..7). info=1
 * is the "primary convergence" success path; info=2,3 are also
 * acceptable (relative-error / gradient-orthogonality convergence);
 * info=0 means bad input; info>=4 means non-convergence (max
 * iterations hit, no further progress, etc.). K26 maps these into
 * the K26ASTRO_FIT_* enum below.
 *
 * Determinism contract: deterministic given the residual callback
 * is deterministic. -O2 -ffp-contract=off -fexcess-precision=standard
 * -frounding-math (matches libk26astro_*). cminpack's LM iteration
 * is deterministic given fixed inputs. */
#ifndef K26ASTRO_FIT_H
#define K26ASTRO_FIT_H

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_FIT_LIB_VERSION "0.1.0"

/* ----- Return codes (mapped from cminpack info) -------------------- */

/* K26ASTRO_FIT_OK covers cminpack info codes 1-4: convergence on
 * residual relative-error (1), parameter relative-error (2), both (3),
 * or Jacobian orthogonal to residuals at machine precision (4 — the
 * standard "gradient-flat success" signal recognized by scipy / GSL /
 * Octave). */
#define K26ASTRO_FIT_OK                  0
#define K26ASTRO_FIT_E_BAD_INPUT        -1   /* invalid args (e.g. m<n) */
#define K26ASTRO_FIT_E_USER_ABORT       -2   /* callback returned negative */
#define K26ASTRO_FIT_E_MAX_ITER         -3   /* maxfev exhausted */
#define K26ASTRO_FIT_E_TOL_TOO_SMALL    -5   /* tolerances at roundoff floor */

/* Max NEQ supported by the canned KFL surface (m residuals, n params). */
#define K26ASTRO_FIT_MAX_N_HANDLE       16    /* parameters */
#define K26ASTRO_FIT_MAX_M_HANDLE      256    /* observations */

/* ----- C-direct API (general residual callback) -------------------- */

/* Residual callback: writes `m` residual values (model(x) - data) at
 * the current `n` parameters `x` into `fvec`. Return 0 on success;
 * non-zero to abort the fit (cminpack will halt with user-abort). */
typedef int (*K26AstroFitFn)(int m, int n,
                             const double *x, double *fvec,
                             void *user);

/* Levenberg-Marquardt with numerical Jacobian (forward differences).
 *
 * Parameters:
 *   f      : residual callback (must be non-NULL)
 *   user   : opaque pointer passed to f
 *   m      : number of residuals (data points); m >= n
 *   n      : number of free parameters; 1 <= n <= m
 *   x      : in: initial parameter guess (length n); out: best-fit
 *            parameters at exit
 *   fvec   : optional output buffer for the final residuals (length m);
 *            may be NULL if the caller doesn't care
 *   tol    : convergence tolerance (rel + abs combined per cminpack
 *            convention; tol = 1e-8 is a sane default for double
 *            precision)
 *
 * Returns K26ASTRO_FIT_OK on convergence; one of K26ASTRO_FIT_E_* on
 * partial completion. `x` is populated regardless. */
int k26astro_fit_lmdif(K26AstroFitFn f, void *user,
                       int m, int n,
                       double *x, double *fvec,
                       double tol);

/* ----- KFL-callable surface (opaque-handle residual problems) ------ */

/* Opaque fit problem handle. Constructs a residual function from a
 * canned family + observation data (xs, ys). Caller owns + destroys. */
typedef struct K26AstroFitProblem K26AstroFitProblem;

/* Number of free parameters of the canned problem. */
int k26astro_fit_problem_n(const K26AstroFitProblem *p);

/* Number of observations (m) registered on the problem. */
int k26astro_fit_problem_m(const K26AstroFitProblem *p);

/* Linear regression: y_i = c0 + c1·x_i.
 * 2 free parameters [c0, c1]. m observations. */
K26AstroFitProblem *k26astro_fit_problem_linear(int m,
                                                 const double *xs,
                                                 const double *ys);

/* Quadratic regression: y_i = c0 + c1·x_i + c2·x_i².
 * 3 free parameters. */
K26AstroFitProblem *k26astro_fit_problem_quadratic(int m,
                                                    const double *xs,
                                                    const double *ys);

/* Power-law fit: y_i = c · x_i^n.
 * 2 free parameters [c, n]. x_i must be positive for the residual
 * computation. */
K26AstroFitProblem *k26astro_fit_problem_power(int m,
                                                const double *xs,
                                                const double *ys);

/* Exponential-decay fit: y_i = c · exp(-α · x_i).
 * 2 free parameters [c, α]. */
K26AstroFitProblem *k26astro_fit_problem_exp(int m,
                                              const double *xs,
                                              const double *ys);

/* Gaussian fit: y_i = amp · exp(-(x_i - μ)² / (2σ²)).
 * 3 free parameters [μ, σ, amp]. */
K26AstroFitProblem *k26astro_fit_problem_gaussian(int m,
                                                   const double *xs,
                                                   const double *ys);

/* Free a problem handle. NULL is a no-op. */
void k26astro_fit_problem_destroy(K26AstroFitProblem *p);

/* Handle-based lmdif. Initial guess + final params in/out at `x`
 * (length k26astro_fit_problem_n(p)). */
int k26astro_fit_lmdif_h(K26AstroFitProblem *p,
                         double *x, double tol);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_FIT_H */
