/* libk26astro_quad - adaptive one-dimensional numerical quadrature.
 *
 * QUADPACK (jacobwilliams modernized fork; BSD-3) wrapped behind a
 * pure-C ABI via ISO_C_BINDING (Fortran 2003). Establishes the K26
 * Fortran-tier library pattern (pure-C ABI surface + ISO_C_BINDING
 * shim + bundled gfortran link).
 *
 * Two surfaces:
 *   1. C-direct API — `k26astro_quad_dqags(K26AstroQuadFn f, void *user,
 *      ...)`. Fully general; any C-side callable can be an integrand.
 *      Used by other K26 libs (e.g. libk26astro_atmos refraction
 *      column integrals, future libk26astro_rad Planck integrals).
 *   2. KFL-callable surface — handle-based via `K26AstroQuadIntegrand`
 *      opaque. KFL programs construct integrands from a fixed set of
 *      canned families (polynomial, exponential, gaussian, Planck,
 *      etc.) and pass the handle to `k26astro_quad_dqags_h`. Until
 *      kflc gains a `fn_ref` argument type, arbitrary user-defined
 *      KFL integrands are not supported; the pre-baked families
 *      cover the typical consumers (refraction, blackbody spectra,
 *      gaussian kernels).
 *
 * User-data threading (C-direct API): K26AstroQuadFn takes an opaque
 * `void *user`. Fortran's `procedure(func)` interface in QUADPACK has
 * no user-data slot; the C-side wrapper threads the (fn, user) pair
 * through thread-local storage during the call and restores any
 * prior pair on return. Nested quadrature is supported.
 *
 * Error model: returns 0 on full convergence; non-zero on incomplete
 * convergence with the QUADPACK `ier` code. Result + abserr are
 * populated even on non-zero return (best-effort estimates).
 *
 * Workspace: each call uses a 500-subdivision stack-allocated
 * workspace (Limit=500, Lenw=2000 doubles). No heap allocation in
 * the integration path; safe inside REFERENCED-mode and KFL `tick`
 * bodies. Integrand-handle constructors heap-allocate the handle;
 * the integrand fn itself runs allocator-free.
 *
 * Determinism contract: deterministic given the integrand is
 * deterministic. -O2 -ffp-contract=off -fexcess-precision=standard
 * -frounding-math (matches the rest of libk26astro_*). */
#ifndef K26ASTRO_QUAD_H
#define K26ASTRO_QUAD_H

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_QUAD_LIB_VERSION "0.1.0"

/* ----- Return codes (passthrough of QUADPACK ier) ------------------- */

#define K26ASTRO_QUAD_OK                 0
#define K26ASTRO_QUAD_E_LIMIT_REACHED    1
#define K26ASTRO_QUAD_E_ROUNDOFF         2
#define K26ASTRO_QUAD_E_BAD_BEHAVIOUR    3
#define K26ASTRO_QUAD_E_NO_CONVERGE      4
#define K26ASTRO_QUAD_E_DIVERGENT        5
#define K26ASTRO_QUAD_E_BAD_INPUT        6

/* Infinite-interval direction codes for DQAGI's `infinite_dir`. */
#define K26ASTRO_QUAD_INF_POS_INFTY      1   /* (bound, +inf) */
#define K26ASTRO_QUAD_INF_NEG_INFTY     -1   /* (-inf, bound) */
#define K26ASTRO_QUAD_INF_BOTH           2   /* (-inf, +inf); `bound` ignored */

/* ----- C-direct API (general integrand callback) -------------------- */

/* Integrand callback. Returns f(x) for the supplied x; user is the
 * opaque pointer passed at quadrature-call time. Must be deterministic
 * for the quadrature to be deterministic. */
typedef double (*K26AstroQuadFn)(double x, void *user);

/* DQAGS — adaptive integration of f(x) over [a, b]. */
int k26astro_quad_dqags(K26AstroQuadFn f, void *user,
                        double a, double b,
                        double epsabs, double epsrel,
                        double *out_result, double *out_abserr);

/* DQAGI — adaptive integration over a (semi-)infinite interval.
 * Direction selected by `infinite_dir`:
 *   K26ASTRO_QUAD_INF_POS_INFTY : (bound, +inf)
 *   K26ASTRO_QUAD_INF_NEG_INFTY : (-inf, bound)
 *   K26ASTRO_QUAD_INF_BOTH      : (-inf, +inf); `bound` ignored */
int k26astro_quad_dqagi(K26AstroQuadFn f, void *user,
                        double bound, int infinite_dir,
                        double epsabs, double epsrel,
                        double *out_result, double *out_abserr);

/* ----- KFL-callable surface (opaque-handle integrands) -------------- */

/* Opaque integrand handle. Constructed by one of the integrand_*
 * factory functions below; consumed by k26astro_quad_dqags_h / _dqagi_h.
 * Caller owns the handle and must call k26astro_quad_integrand_destroy
 * to free it. */
typedef struct K26AstroQuadIntegrand K26AstroQuadIntegrand;

/* Constant integrand: f(x) = c. Useful for testing. */
K26AstroQuadIntegrand *k26astro_quad_integrand_const(double c);

/* Quadratic polynomial: f(x) = c0 + c1·x + c2·x². */
K26AstroQuadIntegrand *k26astro_quad_integrand_poly2(double c0, double c1, double c2);

/* Cubic polynomial: f(x) = c0 + c1·x + c2·x² + c3·x³. */
K26AstroQuadIntegrand *k26astro_quad_integrand_poly3(double c0, double c1,
                                                     double c2, double c3);

/* Power law: f(x) = c · x^n (n may be real-valued, including negative
 * for inverse-power decay; x must stay positive for non-integer n). */
K26AstroQuadIntegrand *k26astro_quad_integrand_power(double c, double n);

/* Exponential decay: f(x) = c · exp(-α·x). For α > 0, integrable on
 * [0, +∞). */
K26AstroQuadIntegrand *k26astro_quad_integrand_exp(double c, double alpha);

/* Gaussian: f(x) = amp · exp(-(x - mu)² / (2·σ²)). Note the standard
 * Gaussian PDF prefactor 1/(σ·√(2π)) is NOT applied — caller chooses
 * `amp` to control normalization. */
K26AstroQuadIntegrand *k26astro_quad_integrand_gaussian(double mu, double sigma, double amp);

/* Planck spectral radiance: f(λ) = (2hc²/λ⁵) · 1/(exp(hc/(λkT)) - 1).
 * λ in metres, T in kelvin, output in W/(m²·sr·m). Used by future
 * libk26astro_rad work; useful here for the Stefan-Boltzmann
 * regression (∫₀^∞ B(λ;T) dλ = σT⁴/π). */
K26AstroQuadIntegrand *k26astro_quad_integrand_planck(double T_K);

/* Free an integrand handle previously returned by any of the
 * constructors above. Calling with NULL is a no-op. */
void k26astro_quad_integrand_destroy(K26AstroQuadIntegrand *h);

/* Handle-based DQAGS — KFL-callable equivalent of dqags. */
int k26astro_quad_dqags_h(K26AstroQuadIntegrand *h,
                          double a, double b,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr);

/* Handle-based DQAGI — KFL-callable equivalent of dqagi. */
int k26astro_quad_dqagi_h(K26AstroQuadIntegrand *h,
                          double bound, int infinite_dir,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_QUAD_H */
