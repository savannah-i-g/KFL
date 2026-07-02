/* k26astro_quad.c — C-side wrapper for libk26astro_quad.
 *
 * Threads (fn, user) through thread-local storage to the Fortran
 * integrand shim in src/k26astro_quad_iface.f90, and implements the
 * KFL-callable opaque-handle integrand surface.
 *
 * The thread-local storage uses C11 `_Thread_local`. Alpine 3.21 musl
 * and Ubuntu 24.04 glibc both support it. Quadrature is sequential
 * within a single call; the TLS isolation only matters across calls
 * from different threads in the same process. */

#include "k26astro_quad/quad.h"
#include "k26astro_quad/quad_consts.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ----- Fortran bind(C) entries (defined in k26astro_quad_iface.f90) -- */

extern void k26astro_quad_dqags_call(double a, double b,
                                     double epsabs, double epsrel,
                                     double *out_result, double *out_abserr,
                                     int *out_neval, int *out_ier);

extern void k26astro_quad_dqagi_call(double bound, int inf,
                                     double epsabs, double epsrel,
                                     double *out_result, double *out_abserr,
                                     int *out_neval, int *out_ier);

/* ----- Thread-local (fn, user) for the Fortran integrand shim ------- */

static _Thread_local K26AstroQuadFn k_tls_fn   = NULL;
static _Thread_local void          *k_tls_user = NULL;

/* C trampoline invoked by the Fortran integrand_shim. Reads the TLS
 * pair set by the C-direct API entry points below. Symbol name is
 * referenced by an `interface ... bind(C, name=...)` block in
 * k26astro_quad_iface.f90. */
double k26astro_quad_tls_trampoline(double x)
{
    return k_tls_fn(x, k_tls_user);
}

/* ----- C-direct API ------------------------------------------------- */

static int k_dqags_inner_(K26AstroQuadFn f, void *user,
                          double a, double b,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr)
{
    K26AstroQuadFn save_fn   = k_tls_fn;
    void          *save_user = k_tls_user;
    k_tls_fn   = f;
    k_tls_user = user;

    double result = 0.0, abserr = 0.0;
    int    neval  = 0,   ier    = 0;
    k26astro_quad_dqags_call(a, b, epsabs, epsrel,
                             &result, &abserr, &neval, &ier);

    k_tls_fn   = save_fn;
    k_tls_user = save_user;

    if (out_result) *out_result = result;
    if (out_abserr) *out_abserr = abserr;
    return ier;
}

static int k_dqagi_inner_(K26AstroQuadFn f, void *user,
                          double bound, int infinite_dir,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr)
{
    K26AstroQuadFn save_fn   = k_tls_fn;
    void          *save_user = k_tls_user;
    k_tls_fn   = f;
    k_tls_user = user;

    double result = 0.0, abserr = 0.0;
    int    neval  = 0,   ier    = 0;
    k26astro_quad_dqagi_call(bound, infinite_dir, epsabs, epsrel,
                             &result, &abserr, &neval, &ier);

    k_tls_fn   = save_fn;
    k_tls_user = save_user;

    if (out_result) *out_result = result;
    if (out_abserr) *out_abserr = abserr;
    return ier;
}

int k26astro_quad_dqags(K26AstroQuadFn f, void *user,
                        double a, double b,
                        double epsabs, double epsrel,
                        double *out_result, double *out_abserr)
{
    if (!f) return K26ASTRO_QUAD_E_BAD_INPUT;
    return k_dqags_inner_(f, user, a, b, epsabs, epsrel,
                          out_result, out_abserr);
}

int k26astro_quad_dqagi(K26AstroQuadFn f, void *user,
                        double bound, int infinite_dir,
                        double epsabs, double epsrel,
                        double *out_result, double *out_abserr)
{
    if (!f) return K26ASTRO_QUAD_E_BAD_INPUT;
    if (infinite_dir != K26ASTRO_QUAD_INF_POS_INFTY &&
        infinite_dir != K26ASTRO_QUAD_INF_NEG_INFTY &&
        infinite_dir != K26ASTRO_QUAD_INF_BOTH) {
        return K26ASTRO_QUAD_E_BAD_INPUT;
    }
    return k_dqagi_inner_(f, user, bound, infinite_dir, epsabs, epsrel,
                          out_result, out_abserr);
}

/* ----- KFL-callable opaque-handle integrand surface ----------------- */

/* Canned-integrand discriminator. The handle's `kind` selects which
 * `params` slots are live; the matching `_at` helper below evaluates. */
typedef enum {
    K_KIND_CONST    = 1,
    K_KIND_POLY2    = 2,
    K_KIND_POLY3    = 3,
    K_KIND_POWER    = 4,
    K_KIND_EXP      = 5,
    K_KIND_GAUSSIAN = 6,
    K_KIND_PLANCK   = 7
} K_IntegrandKind;

struct K26AstroQuadIntegrand {
    K_IntegrandKind kind;
    double          params[8];
};

/* Planck constants - hex-literal IEEE-754. CODATA 2018 values:
 *   h = 6.62607015e-34 J·s (exact, SI definition since 2019)
 *   c = 2.99792458e8  m/s  (exact, SI definition since 1983)
 *   k_B = 1.380649e-23 J/K (exact, SI definition since 2019)
 * Decimal expansions converted via IEEE-754 round-to-nearest. */
#define K_PLANCK_H_BITS    0x390F00CB18C5398EULL  /* 6.62607015e-34 */
#define K_PLANCK_C_BITS    0x41B1DE784A000000ULL  /* 2.99792458e8 */
#define K_PLANCK_KB_BITS   0x3BFA28B4554D5D89ULL  /* 1.380649e-23 */

static double k_planck_h_(void)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = K_PLANCK_H_BITS;
    return cvt.d;
}

static double k_planck_c_(void)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = K_PLANCK_C_BITS;
    return cvt.d;
}

static double k_planck_kb_(void)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = K_PLANCK_KB_BITS;
    return cvt.d;
}

/* Single dispatch from K26AstroQuadFn into the canned integrand
 * family by inspecting the handle's `kind`. */
static double k_integrand_eval_(double x, void *user)
{
    const K26AstroQuadIntegrand *h = (const K26AstroQuadIntegrand *)user;
    if (!h) return 0.0;
    switch (h->kind) {
    case K_KIND_CONST:
        return h->params[0];
    case K_KIND_POLY2:
        return h->params[0] + x * (h->params[1] + x * h->params[2]);
    case K_KIND_POLY3:
        return h->params[0] +
               x * (h->params[1] +
                    x * (h->params[2] + x * h->params[3]));
    case K_KIND_POWER:
        /* c · x^n; handle x = 0 with the standard convention for
         * positive n (0) and negative n (returns inf, propagated). */
        if (x == 0.0) {
            if (h->params[1] > 0.0) return 0.0;
            if (h->params[1] == 0.0) return h->params[0];
            return INFINITY;
        }
        return h->params[0] * pow(x, h->params[1]);
    case K_KIND_EXP:
        return h->params[0] * exp(-h->params[1] * x);
    case K_KIND_GAUSSIAN: {
        const double dx = x - h->params[0];
        const double s2 = h->params[1] * h->params[1];
        return h->params[2] * exp(-(dx * dx) / (2.0 * s2));
    }
    case K_KIND_PLANCK: {
        /* B(λ; T) = (2hc² / λ⁵) · 1 / (exp(hc/(λkT)) - 1). */
        const double lambda = x;
        if (lambda <= 0.0) return 0.0;
        const double T  = h->params[0];
        const double h_ = k_planck_h_();
        const double c_ = k_planck_c_();
        const double kB = k_planck_kb_();
        const double l5 = lambda * lambda * lambda * lambda * lambda;
        const double a  = 2.0 * h_ * c_ * c_ / l5;
        const double e  = h_ * c_ / (lambda * kB * T);
        /* Guard against overflow in exp(); for huge e, return 0. */
        if (e > 700.0) return 0.0;
        return a / (exp(e) - 1.0);
    }
    }
    return 0.0;
}

static K26AstroQuadIntegrand *k_alloc_(K_IntegrandKind kind)
{
    K26AstroQuadIntegrand *h = (K26AstroQuadIntegrand *)
        calloc(1, sizeof *h);
    if (h) h->kind = kind;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_const(double c)
{
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_CONST);
    if (!h) return NULL;
    h->params[0] = c;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_poly2(double c0, double c1, double c2)
{
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_POLY2);
    if (!h) return NULL;
    h->params[0] = c0;
    h->params[1] = c1;
    h->params[2] = c2;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_poly3(double c0, double c1,
                                                     double c2, double c3)
{
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_POLY3);
    if (!h) return NULL;
    h->params[0] = c0;
    h->params[1] = c1;
    h->params[2] = c2;
    h->params[3] = c3;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_power(double c, double n)
{
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_POWER);
    if (!h) return NULL;
    h->params[0] = c;
    h->params[1] = n;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_exp(double c, double alpha)
{
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_EXP);
    if (!h) return NULL;
    h->params[0] = c;
    h->params[1] = alpha;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_gaussian(double mu, double sigma, double amp)
{
    if (sigma <= 0.0) return NULL;
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_GAUSSIAN);
    if (!h) return NULL;
    h->params[0] = mu;
    h->params[1] = sigma;
    h->params[2] = amp;
    return h;
}

K26AstroQuadIntegrand *k26astro_quad_integrand_planck(double T_K)
{
    if (T_K <= 0.0) return NULL;
    K26AstroQuadIntegrand *h = k_alloc_(K_KIND_PLANCK);
    if (!h) return NULL;
    h->params[0] = T_K;
    return h;
}

void k26astro_quad_integrand_destroy(K26AstroQuadIntegrand *h)
{
    free(h);
}

int k26astro_quad_dqags_h(K26AstroQuadIntegrand *h,
                          double a, double b,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr)
{
    if (!h) return K26ASTRO_QUAD_E_BAD_INPUT;
    return k_dqags_inner_(k_integrand_eval_, h, a, b, epsabs, epsrel,
                          out_result, out_abserr);
}

int k26astro_quad_dqagi_h(K26AstroQuadIntegrand *h,
                          double bound, int infinite_dir,
                          double epsabs, double epsrel,
                          double *out_result, double *out_abserr)
{
    if (!h) return K26ASTRO_QUAD_E_BAD_INPUT;
    if (infinite_dir != K26ASTRO_QUAD_INF_POS_INFTY &&
        infinite_dir != K26ASTRO_QUAD_INF_NEG_INFTY &&
        infinite_dir != K26ASTRO_QUAD_INF_BOTH) {
        return K26ASTRO_QUAD_E_BAD_INPUT;
    }
    return k_dqagi_inner_(k_integrand_eval_, h, bound, infinite_dir,
                          epsabs, epsrel, out_result, out_abserr);
}
