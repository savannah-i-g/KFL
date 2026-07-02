/* k26astro_fit.c - C-side wrapper for libk26astro_fit.
 *
 * cminpack's callback signature already includes a user-data pointer
 * (cminpack's `p` arg), so unlike the Fortran-backed wrappers in
 * libk26astro_quad / libk26astro_ode no thread-local storage is
 * needed. The K26 callback shape (K26AstroFitFn) takes the user
 * pointer at the END of the arg list; cminpack passes it FIRST. A
 * trivial struct + adapter bridges them.
 *
 * KFL-callable opaque problem handles bundle a residual family + the
 * observation arrays. The adapter computes residuals on each fit
 * iteration. */

#include "k26astro_fit/fit.h"
#include "../src/upstream/cminpack/cminpack.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ----- K26 → cminpack callback adapter ------------------------------ */

typedef struct {
    K26AstroFitFn user_fn;
    void         *user_data;
} K_FitAdapter;

static int k_cm_adapter_(void *p, int m, int n,
                         const double *x, double *fvec, int iflag)
{
    (void)iflag;  /* K26 callback doesn't distinguish residual vs Jacobian
                   * because we always use numerical differencing
                   * (fdjac2 inside lmdif). */
    const K_FitAdapter *a = (const K_FitAdapter *)p;
    return a->user_fn(m, n, x, fvec, a->user_data);
}

/* ----- info → K26 code mapping ------------------------------------- */

static int k_map_info_(int info)
{
    /* info 1/2/3 = convergence on residual / parameter / both relative
     * tolerances. info 4 = Jacobian orthogonal to residuals at machine
     * precision — the standard "gradient-flat" success signal (scipy,
     * GSL, Octave all treat as success). Only 5/6/7 are failures. */
    if (info == 1 || info == 2 || info == 3 || info == 4) {
        return K26ASTRO_FIT_OK;
    }
    if (info == 0) return K26ASTRO_FIT_E_BAD_INPUT;
    if (info == 5) return K26ASTRO_FIT_E_MAX_ITER;
    if (info == 6 || info == 7) return K26ASTRO_FIT_E_TOL_TOO_SMALL;
    if (info < 0) return K26ASTRO_FIT_E_USER_ABORT;
    return K26ASTRO_FIT_E_BAD_INPUT;
}

/* ----- C-direct API ------------------------------------------------- */

int k26astro_fit_lmdif(K26AstroFitFn f, void *user,
                       int m, int n,
                       double *x, double *fvec,
                       double tol)
{
    if (!f || !x || n <= 0 || m < n || tol < 0.0) {
        return K26ASTRO_FIT_E_BAD_INPUT;
    }

    /* Workspace: cminpack lmdif1 wants lwa >= m*n + 5*n + m. */
    const int lwa = m * n + 5 * n + m;
    int    *iwa = (int *)calloc((size_t)n,   sizeof(int));
    double *wa  = (double *)calloc((size_t)lwa, sizeof(double));
    double *fvec_local = NULL;
    double *fvec_use   = fvec;
    if (!fvec_use) {
        fvec_local = (double *)calloc((size_t)m, sizeof(double));
        fvec_use   = fvec_local;
    }
    if (!iwa || !wa || !fvec_use) {
        free(iwa); free(wa); free(fvec_local);
        return K26ASTRO_FIT_E_BAD_INPUT;
    }

    K_FitAdapter adapter = { f, user };
    int info = lmdif1(k_cm_adapter_, &adapter, m, n, x, fvec_use,
                      tol, iwa, wa, lwa);

    free(iwa);
    free(wa);
    free(fvec_local);
    return k_map_info_(info);
}

/* ----- KFL-callable opaque problem surface -------------------------- */

typedef enum {
    K_FIT_KIND_LINEAR    = 1,
    K_FIT_KIND_QUADRATIC = 2,
    K_FIT_KIND_POWER     = 3,
    K_FIT_KIND_EXP       = 4,
    K_FIT_KIND_GAUSSIAN  = 5
} K_FitKind;

struct K26AstroFitProblem {
    K_FitKind kind;
    int       m;     /* observations */
    int       n;     /* parameters */
    double   *xs;    /* heap-allocated, length m */
    double   *ys;    /* heap-allocated, length m */
};

int k26astro_fit_problem_n(const K26AstroFitProblem *p) { return p ? p->n : 0; }
int k26astro_fit_problem_m(const K26AstroFitProblem *p) { return p ? p->m : 0; }

static K26AstroFitProblem *k_alloc_(K_FitKind kind, int m, int n,
                                    const double *xs, const double *ys)
{
    if (m <= 0 || n <= 0 || m < n) return NULL;
    if (m > K26ASTRO_FIT_MAX_M_HANDLE) return NULL;
    if (n > K26ASTRO_FIT_MAX_N_HANDLE) return NULL;
    K26AstroFitProblem *p = (K26AstroFitProblem *)calloc(1, sizeof *p);
    if (!p) return NULL;
    p->kind = kind;
    p->m    = m;
    p->n    = n;
    p->xs   = (double *)calloc((size_t)m, sizeof(double));
    p->ys   = (double *)calloc((size_t)m, sizeof(double));
    if (!p->xs || !p->ys) {
        free(p->xs); free(p->ys); free(p);
        return NULL;
    }
    memcpy(p->xs, xs, (size_t)m * sizeof(double));
    memcpy(p->ys, ys, (size_t)m * sizeof(double));
    return p;
}

K26AstroFitProblem *k26astro_fit_problem_linear(int m, const double *xs, const double *ys)
{
    if (!xs || !ys) return NULL;
    return k_alloc_(K_FIT_KIND_LINEAR, m, 2, xs, ys);
}

K26AstroFitProblem *k26astro_fit_problem_quadratic(int m, const double *xs, const double *ys)
{
    if (!xs || !ys) return NULL;
    return k_alloc_(K_FIT_KIND_QUADRATIC, m, 3, xs, ys);
}

K26AstroFitProblem *k26astro_fit_problem_power(int m, const double *xs, const double *ys)
{
    if (!xs || !ys) return NULL;
    return k_alloc_(K_FIT_KIND_POWER, m, 2, xs, ys);
}

K26AstroFitProblem *k26astro_fit_problem_exp(int m, const double *xs, const double *ys)
{
    if (!xs || !ys) return NULL;
    return k_alloc_(K_FIT_KIND_EXP, m, 2, xs, ys);
}

K26AstroFitProblem *k26astro_fit_problem_gaussian(int m, const double *xs, const double *ys)
{
    if (!xs || !ys) return NULL;
    return k_alloc_(K_FIT_KIND_GAUSSIAN, m, 3, xs, ys);
}

void k26astro_fit_problem_destroy(K26AstroFitProblem *p)
{
    if (!p) return;
    free(p->xs);
    free(p->ys);
    free(p);
}

/* Residual dispatcher for canned families. Each writes (model - data)
 * into fvec[i] = predict(x; params) - ys[i]. */
static int k_fit_residual_(int m, int n, const double *params,
                           double *fvec, void *user)
{
    (void)n;
    const K26AstroFitProblem *p = (const K26AstroFitProblem *)user;
    if (!p || m != p->m) {
        for (int i = 0; i < m; ++i) fvec[i] = 0.0;
        return 0;
    }
    const double *xs = p->xs;
    const double *ys = p->ys;
    switch (p->kind) {
    case K_FIT_KIND_LINEAR: {
        /* params = [c0, c1]; model(x) = c0 + c1·x */
        const double c0 = params[0];
        const double c1 = params[1];
        for (int i = 0; i < m; ++i) {
            fvec[i] = (c0 + c1 * xs[i]) - ys[i];
        }
        return 0;
    }
    case K_FIT_KIND_QUADRATIC: {
        /* params = [c0, c1, c2]; model(x) = c0 + c1·x + c2·x² */
        const double c0 = params[0];
        const double c1 = params[1];
        const double c2 = params[2];
        for (int i = 0; i < m; ++i) {
            const double xi = xs[i];
            fvec[i] = (c0 + xi * (c1 + xi * c2)) - ys[i];
        }
        return 0;
    }
    case K_FIT_KIND_POWER: {
        /* params = [c, n]; model(x) = c · x^n. x must be positive. */
        const double c = params[0];
        const double n_pow = params[1];
        for (int i = 0; i < m; ++i) {
            const double xi = xs[i];
            const double model = (xi > 0.0) ? c * pow(xi, n_pow) : 0.0;
            fvec[i] = model - ys[i];
        }
        return 0;
    }
    case K_FIT_KIND_EXP: {
        /* params = [c, α]; model(x) = c · exp(-α·x) */
        const double c = params[0];
        const double a = params[1];
        for (int i = 0; i < m; ++i) {
            fvec[i] = c * exp(-a * xs[i]) - ys[i];
        }
        return 0;
    }
    case K_FIT_KIND_GAUSSIAN: {
        /* params = [μ, σ, amp]; model(x) = amp·exp(-(x-μ)²/(2σ²)) */
        const double mu    = params[0];
        const double sigma = params[1];
        const double amp   = params[2];
        const double s2    = sigma * sigma;
        if (s2 == 0.0) {
            for (int i = 0; i < m; ++i) fvec[i] = -ys[i];
            return 0;
        }
        for (int i = 0; i < m; ++i) {
            const double dx = xs[i] - mu;
            fvec[i] = amp * exp(-(dx * dx) / (2.0 * s2)) - ys[i];
        }
        return 0;
    }
    }
    for (int i = 0; i < m; ++i) fvec[i] = 0.0;
    return 0;
}

int k26astro_fit_lmdif_h(K26AstroFitProblem *p, double *x, double tol)
{
    if (!p || !x) return K26ASTRO_FIT_E_BAD_INPUT;
    return k26astro_fit_lmdif(k_fit_residual_, p, p->m, p->n,
                              x, NULL, tol);
}
