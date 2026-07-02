/* k26astro_ode.c — C-side wrapper for libk26astro_ode.
 *
 * Threads (rhs, user, n) through thread-local storage to the Fortran
 * rhs_shim in src/k26astro_ode_iface.f90, and implements the
 * KFL-callable opaque-handle ODE-problem surface.
 *
 * Thread-local storage uses C11 _Thread_local. Same pattern as
 * libk26astro_quad. */

#include "k26astro_ode/ode.h"
#include "k26astro_ode/ode_consts.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ----- Fortran bind(C) entry (defined in k26astro_ode_iface.f90) --- */

extern void k26astro_ode_lsoda_call(int n,
                                    const double *y0,
                                    double t0, double tf,
                                    double rtol, double atol,
                                    double *out_y,
                                    int *out_istate);

/* ----- Thread-local (rhs, user, n) for the Fortran rhs shim --------- */

static _Thread_local K26AstroOdeRhsFn k_tls_rhs  = NULL;
static _Thread_local void            *k_tls_user = NULL;
static _Thread_local int              k_tls_n    = 0;

/* C trampoline invoked by the Fortran rhs_shim. Reads the TLS pair
 * set by the C-direct entry below. Signature matches the F-side
 * interface block in k26astro_ode_iface.f90. */
void k26astro_ode_tls_rhs_trampoline(int n, double t,
                                     const double *y, double *ydot)
{
    if (k_tls_rhs) {
        (void)k_tls_rhs(n, t, y, ydot, k_tls_user);
    } else {
        /* No RHS bound — fill ydot with zeros (degenerate). */
        for (int i = 0; i < n; ++i) ydot[i] = 0.0;
    }
}

/* ----- DLSODA istate → K26 return code mapping --------------------- */

static int k_map_istate_(int istate)
{
    if (istate == 2) return K26ASTRO_ODE_OK;
    /* istate values -1 .. -7 already match our K26ASTRO_ODE_E_* enum. */
    return istate;
}

/* ----- C-direct API ------------------------------------------------- */

int k26astro_ode_lsoda_solve(K26AstroOdeRhsFn rhs, void *user,
                             int n,
                             const double *y0,
                             double t0, double tf,
                             double rtol, double atol,
                             double *out_y)
{
    if (!rhs || !y0 || !out_y || n <= 0 || rtol < 0.0 || atol < 0.0) {
        return K26ASTRO_ODE_E_BAD_INPUT;
    }

    /* Save/restore so nested solves compose. */
    K26AstroOdeRhsFn save_rhs  = k_tls_rhs;
    void            *save_user = k_tls_user;
    int              save_n    = k_tls_n;
    k_tls_rhs  = rhs;
    k_tls_user = user;
    k_tls_n    = n;

    int istate = 0;
    k26astro_ode_lsoda_call(n, y0, t0, tf, rtol, atol, out_y, &istate);

    k_tls_rhs  = save_rhs;
    k_tls_user = save_user;
    k_tls_n    = save_n;

    return k_map_istate_(istate);
}

/* ----- KFL-callable opaque-handle problem surface ------------------- */

typedef enum {
    K_ODE_KIND_LINEAR          = 1,
    K_ODE_KIND_HARMONIC        = 2,
    K_ODE_KIND_DAMPED_HARMONIC = 3,
    K_ODE_KIND_VAN_DER_POL     = 4,
    K_ODE_KIND_LOTKA_VOLTERRA  = 5,
    K_ODE_KIND_ROBERTSON       = 6
} K_OdeKind;

struct K26AstroOdeProblem {
    K_OdeKind kind;
    int       n;
    double    params[8];
    double   *matrix; /* linear A row-major, length n*n; NULL for non-linear */
    double   *offset; /* linear b, length n; NULL for non-linear */
};

int k26astro_ode_problem_n(const K26AstroOdeProblem *p)
{
    return p ? p->n : 0;
}

static K26AstroOdeProblem *k_alloc_(K_OdeKind kind, int n)
{
    if (n <= 0 || n > K26ASTRO_ODE_MAX_NEQ_HANDLE) return NULL;
    K26AstroOdeProblem *p = (K26AstroOdeProblem *)calloc(1, sizeof *p);
    if (!p) return NULL;
    p->kind = kind;
    p->n    = n;
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_linear(int n,
                                                const double *A_row_major,
                                                const double *b)
{
    if (!A_row_major || !b) return NULL;
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_LINEAR, n);
    if (!p) return NULL;
    p->matrix = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    p->offset = (double *)calloc((size_t)n, sizeof(double));
    if (!p->matrix || !p->offset) {
        free(p->matrix);
        free(p->offset);
        free(p);
        return NULL;
    }
    memcpy(p->matrix, A_row_major, (size_t)n * (size_t)n * sizeof(double));
    memcpy(p->offset, b,           (size_t)n * sizeof(double));
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_harmonic(double omega)
{
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_HARMONIC, 2);
    if (!p) return NULL;
    p->params[0] = omega;
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_damped_harmonic(double omega, double gamma)
{
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_DAMPED_HARMONIC, 2);
    if (!p) return NULL;
    p->params[0] = omega;
    p->params[1] = gamma;
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_van_der_pol(double mu)
{
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_VAN_DER_POL, 2);
    if (!p) return NULL;
    p->params[0] = mu;
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_lotka_volterra(double alpha, double beta,
                                                        double delta, double gamma)
{
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_LOTKA_VOLTERRA, 2);
    if (!p) return NULL;
    p->params[0] = alpha;
    p->params[1] = beta;
    p->params[2] = delta;
    p->params[3] = gamma;
    return p;
}

K26AstroOdeProblem *k26astro_ode_problem_robertson(double k1, double k2, double k3)
{
    K26AstroOdeProblem *p = k_alloc_(K_ODE_KIND_ROBERTSON, 3);
    if (!p) return NULL;
    p->params[0] = k1;
    p->params[1] = k2;
    p->params[2] = k3;
    return p;
}

void k26astro_ode_problem_destroy(K26AstroOdeProblem *p)
{
    if (!p) return;
    free(p->matrix);
    free(p->offset);
    free(p);
}

/* RHS dispatcher for canned families. Reads from the TLS-bound
 * problem pointer (set by lsoda_solve_h). */
static int k_problem_rhs_(int n, double t,
                          const double *y, double *ydot, void *user)
{
    (void)t;
    const K26AstroOdeProblem *p = (const K26AstroOdeProblem *)user;
    if (!p || n != p->n) {
        for (int i = 0; i < n; ++i) ydot[i] = 0.0;
        return 0;
    }
    switch (p->kind) {
    case K_ODE_KIND_LINEAR: {
        /* ydot = A·y + b. A is row-major n*n. */
        for (int i = 0; i < n; ++i) {
            double s = p->offset[i];
            const double *row = &p->matrix[i * n];
            for (int j = 0; j < n; ++j) s += row[j] * y[j];
            ydot[i] = s;
        }
        return 0;
    }
    case K_ODE_KIND_HARMONIC: {
        /* y = [x, v]; ydot = [v, -ω²·x] */
        const double omega = p->params[0];
        ydot[0] = y[1];
        ydot[1] = -omega * omega * y[0];
        return 0;
    }
    case K_ODE_KIND_DAMPED_HARMONIC: {
        /* y = [x, v]; ydot = [v, -ω²·x - 2γ·v] */
        const double omega = p->params[0];
        const double gamma = p->params[1];
        ydot[0] = y[1];
        ydot[1] = -omega * omega * y[0] - 2.0 * gamma * y[1];
        return 0;
    }
    case K_ODE_KIND_VAN_DER_POL: {
        /* y = [x, v]; ydot = [v, μ·(1-x²)·v - x] */
        const double mu = p->params[0];
        ydot[0] = y[1];
        ydot[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
        return 0;
    }
    case K_ODE_KIND_LOTKA_VOLTERRA: {
        /* y = [prey, predator]
         *   dy[0] =  α·prey - β·prey·predator
         *   dy[1] =  δ·prey·predator - γ·predator */
        const double alpha = p->params[0];
        const double beta  = p->params[1];
        const double delta = p->params[2];
        const double gamma = p->params[3];
        ydot[0] =  alpha * y[0] - beta  * y[0] * y[1];
        ydot[1] =  delta * y[0] * y[1] - gamma * y[1];
        return 0;
    }
    case K_ODE_KIND_ROBERTSON: {
        /* y = [A, B, C]
         *   dy[0] = -k1·A + k2·B·C
         *   dy[1] =  k1·A - k2·B·C - k3·B²
         *   dy[2] =                  k3·B² */
        const double k1 = p->params[0];
        const double k2 = p->params[1];
        const double k3 = p->params[2];
        const double dyB_sq = k3 * y[1] * y[1];
        ydot[0] = -k1 * y[0] + k2 * y[1] * y[2];
        ydot[2] =                                 dyB_sq;
        ydot[1] = -ydot[0] - ydot[2];
        return 0;
    }
    }
    for (int i = 0; i < n; ++i) ydot[i] = 0.0;
    return 0;
}

int k26astro_ode_lsoda_solve_h(K26AstroOdeProblem *p,
                               const double *y0,
                               double t0, double tf,
                               double rtol, double atol,
                               double *out_y)
{
    if (!p) return K26ASTRO_ODE_E_BAD_INPUT;
    return k26astro_ode_lsoda_solve(k_problem_rhs_, p, p->n, y0,
                                    t0, tf, rtol, atol, out_y);
}
