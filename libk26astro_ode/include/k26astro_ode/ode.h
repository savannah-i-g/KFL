/* libk26astro_ode — initial-value ODE solver via ODEPACK / DLSODA.
 *
 * DLSODA (Petzold-Hindmarsh 1983, public domain LLNL release) is the
 * stiff/non-stiff automatic-switching variant of LSODE. It detects
 * stiffness on the fly and switches between an Adams predictor-
 * corrector (non-stiff) and a BDF method with Newton iteration
 * (stiff). Where K26's IAS15 shrinks dt to picoseconds near close
 * encounters, LSODA switches solver and uses macroscopic steps;
 * complementary integrator.
 *
 * jacobwilliams/odepack modern Fortran refactor vendored under
 * src/upstream/odepack/ (35 kLOC; public domain). ISO_C_BINDING
 * wrapper exposes a plain-C ABI.
 *
 * Two surfaces:
 *   1. C-direct API — `k26astro_ode_lsoda_solve(K26AstroOdeRhsFn rhs,
 *      void *user, int n, ...)`. Fully general; any C-side RHS
 *      callback.
 *   2. KFL-callable opaque API - `K26AstroOdeProblem` handle
 *      constructed via canned RHS families (linear, harmonic, damped
 *      harmonic, Van der Pol, Lotka-Volterra, Robertson). Arbitrary
 *      KFL-side RHS callbacks require a kflc `fn_ref` arg-type that
 *      is not yet provided; the pre-baked families cover typical
 *      consumers (stiff chemistry, oscillator dynamics,
 *      predator-prey populations, Van der Pol benchmark).
 *
 * User-data threading (C-direct API): K26AstroOdeRhsFn takes an
 * opaque `void *user`. The Fortran DLSODA entry uses an
 * `external f` argument with no user-data slot; the C-side wrapper
 * threads (rhs, user) through thread-local storage, and the Fortran
 * shim calls back through a C trampoline that reads the TLS. Nested
 * solves are supported via save/restore on entry.
 *
 * Error model: returns 0 on full convergence (DLSODA istate=2 at
 * exit); non-zero passes through DLSODA's `istate` codes. Result
 * vector is populated even on partial completion.
 *
 * Determinism contract: deterministic given the RHS is deterministic.
 * -O2 -ffp-contract=off -fexcess-precision=standard -frounding-math
 * (matches libk26astro_*). DLSODA's stiff/non-stiff switching is
 * deterministic given fixed inputs (switching is decided by error
 * estimates, not thread scheduling). */
#ifndef K26ASTRO_ODE_H
#define K26ASTRO_ODE_H

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_ODE_LIB_VERSION "0.1.0"

/* ----- Return codes (passthrough of DLSODA istate at exit) ----------
 * DLSODA istate convention:
 *   2  — successful return (normal completion)
 *  -1  — excess work done (more than `mxstep` internal steps needed)
 *  -2  — too much accuracy requested (rtol/atol too tight)
 *  -3  — illegal input (verified at C-wrapper entry; should not reach
 *        DLSODA)
 *  -4  — repeated convergence failures (likely bad Jacobian or
 *        singularity)
 *  -5  — repeated error-test failures
 *  -6  — pure absolute-error tolerance failure (atol*y(i)=0 for some i)
 *  -7  — work-array length insufficient (rwork/iwork too short)
 *
 * K26 maps these to non-zero negatives in the same convention; 0 on
 * success (note: DLSODA returns istate=2 on success, we report 0). */
#define K26ASTRO_ODE_OK                  0
#define K26ASTRO_ODE_E_EXCESS_WORK      -1
#define K26ASTRO_ODE_E_TOO_TIGHT        -2
#define K26ASTRO_ODE_E_BAD_INPUT        -3
#define K26ASTRO_ODE_E_CONVERGENCE      -4
#define K26ASTRO_ODE_E_ERROR_TEST       -5
#define K26ASTRO_ODE_E_ATOL_FAILURE     -6
#define K26ASTRO_ODE_E_WORKSPACE        -7

/* Maximum NEQ supported by the canned KFL surface. The C-direct API
 * is unconstrained beyond available memory. */
#define K26ASTRO_ODE_MAX_NEQ_HANDLE     16

/* ----- C-direct API (general RHS callback) -------------------------- */

/* RHS callback: writes dy/dt at (t, y) into ydot.
 * Return 0 on success; non-zero signals an integration abort (DLSODA
 * does not directly consume this signal — the wrapper currently
 * ignores the return value, but consumers can use it for diagnostic
 * logging). */
typedef int (*K26AstroOdeRhsFn)(int n, double t,
                                const double *y, double *ydot,
                                void *user);

/* DLSODA — integrate y' = f(t, y) from t0 to tf with stiff/non-stiff
 * automatic switching.
 *
 * Parameters:
 *   rhs        : RHS callback (must be non-NULL)
 *   user       : opaque pointer passed to rhs on each invocation
 *   n          : dimension of the state vector (n > 0)
 *   y0         : initial state at t0 (length n)
 *   t0, tf     : integration interval. tf may be greater or less than
 *                t0 (DLSODA supports backwards integration when tf < t0)
 *   rtol, atol : scalar relative and absolute tolerances. Convergence
 *                gate is roughly |error_i| <= rtol*|y_i| + atol.
 *                See ode_consts.h for sane defaults.
 *   out_y      : caller-allocated output buffer (length n). On
 *                successful return, contains y(tf).
 *
 * Returns K26ASTRO_ODE_OK on success; one of K26ASTRO_ODE_E_* on
 * partial completion. out_y is populated regardless. */
int k26astro_ode_lsoda_solve(K26AstroOdeRhsFn rhs, void *user,
                             int n,
                             const double *y0,
                             double t0, double tf,
                             double rtol, double atol,
                             double *out_y);

/* ----- KFL-callable surface (opaque-handle RHS) --------------------- */

/* Opaque ODE problem handle. Constructed via integrand_*-style
 * factory functions below; consumed by k26astro_ode_lsoda_solve_h.
 * Caller owns the handle and must call k26astro_ode_problem_destroy. */
typedef struct K26AstroOdeProblem K26AstroOdeProblem;

/* Number of state-vector components of the canned problem. */
int k26astro_ode_problem_n(const K26AstroOdeProblem *p);

/* Linear: dy/dt = A·y + b. A is row-major n*n; b is length n.
 * Used for stiff linear systems, control-theory problems, exponential
 * decay matrices. n must satisfy 1 <= n <= K26ASTRO_ODE_MAX_NEQ_HANDLE. */
K26AstroOdeProblem *k26astro_ode_problem_linear(int n,
                                                const double *A_row_major,
                                                const double *b);

/* Harmonic oscillator: y = [x, v], dy/dt = [v, -ω²·x].
 * Non-stiff. n=2. */
K26AstroOdeProblem *k26astro_ode_problem_harmonic(double omega);

/* Damped harmonic: y = [x, v], dy/dt = [v, -ω²·x - 2γ·v].
 * Non-stiff for normal γ; becomes stiff for ω >> γ. n=2. */
K26AstroOdeProblem *k26astro_ode_problem_damped_harmonic(double omega, double gamma);

/* Van der Pol oscillator: y = [x, v], dy/dt = [v, μ·(1-x²)·v - x].
 * Increasingly stiff as μ grows; μ=1000 is the canonical LSODA
 * stiffness benchmark. n=2. */
K26AstroOdeProblem *k26astro_ode_problem_van_der_pol(double mu);

/* Lotka-Volterra predator-prey: y = [prey, predator],
 *   dy[0]/dt =  α·prey - β·prey·predator
 *   dy[1]/dt =  δ·prey·predator - γ·predator
 * Non-stiff but produces limit cycles useful for regression. n=2. */
K26AstroOdeProblem *k26astro_ode_problem_lotka_volterra(double alpha, double beta,
                                                        double delta, double gamma);

/* Robertson chemistry: y = [A, B, C] concentrations,
 *   dy[0]/dt = -k1·A + k2·B·C
 *   dy[1]/dt =  k1·A - k2·B·C - k3·B²
 *   dy[2]/dt =                  k3·B²
 * Classic stiff chemistry benchmark; rate constants span 7 orders of
 * magnitude. Default k1=0.04, k2=1e4, k3=3e7. n=3. */
K26AstroOdeProblem *k26astro_ode_problem_robertson(double k1, double k2, double k3);

/* Free a problem handle previously returned by any constructor. NULL
 * is a no-op. */
void k26astro_ode_problem_destroy(K26AstroOdeProblem *p);

/* Handle-based DLSODA — KFL-callable equivalent of lsoda_solve.
 * `y0` and `out_y` must each be length k26astro_ode_problem_n(p). */
int k26astro_ode_lsoda_solve_h(K26AstroOdeProblem *p,
                               const double *y0,
                               double t0, double tf,
                               double rtol, double atol,
                               double *out_y);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_ODE_H */
