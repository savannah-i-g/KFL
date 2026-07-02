/* libk26compute — C numerics core for K26 research tooling.
 *
 * This library provides vector and matrix algebra (with LU-based
 * inverse, determinant, and linear-solve), descriptive and inferential
 * statistics (mean, variance, std, quantile, Pearson correlation, OLS
 * linear regression with uncertainty, chi-squared statistic), ODE
 * solvers (classical RK4 and adaptive Dormand-Prince RK45), 1-D and
 * N-D optimisation (Brent's method and Nelder-Mead simplex), and a
 * deterministic RNG with explicit seed (SplitMix64 plus Box-Muller
 * normals).
 *
 * Determinism contract:
 *   - No hidden global RNG state. Every primitive that touches
 *     random bits takes a uint64_t seed or a K26CRng pointer the
 *     caller owns.
 *   - All operations on the same inputs produce the same outputs
 *     bit-for-bit (modulo IEEE-754 rounding implementation).
 *
 * Not provided: GPU offload, sparse matrices, FFT, symbolic algebra,
 * dataset (CSV / Parquet / HDF5) loaders, eigendecomposition for
 * general matrices, nonlinear least squares with constraints. */
#ifndef K26COMPUTE_H
#define K26COMPUTE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define K26COMPUTE_LIB_VERSION "0.1.0"

/* ---- Status codes ---- */
typedef enum {
    K26C_OK = 0,
    K26C_ERR_OOM,
    K26C_ERR_INVAL,        /* invalid arg / dim mismatch / NaN input */
    K26C_ERR_SINGULAR,     /* matrix singular / pivot below tolerance */
    K26C_ERR_CONVERGE,     /* iterative method failed to converge */
    K26C_ERR_RANGE,        /* result outside representable range */
    K26C_ERR_INTERNAL,
} K26CStatus;

const char *k26compute_status_str(K26CStatus s);

/* ---- Vectors (heap-owned 1-D array of doubles) ---- */
typedef struct {
    double *data;
    size_t  n;
} K26CVector;

K26CStatus k26c_vec_alloc(K26CVector *v, size_t n);
void       k26c_vec_free (K26CVector *v);
K26CStatus k26c_vec_from (K26CVector *v, const double *src, size_t n);
K26CStatus k26c_vec_set  (K26CVector *v, size_t i, double x);
double     k26c_vec_get  (const K26CVector *v, size_t i);
K26CStatus k26c_vec_fill (K26CVector *v, double x);
K26CStatus k26c_vec_copy (K26CVector *dst, const K26CVector *src);

K26CStatus k26c_vec_add  (K26CVector *out, const K26CVector *a, const K26CVector *b);
K26CStatus k26c_vec_sub  (K26CVector *out, const K26CVector *a, const K26CVector *b);
K26CStatus k26c_vec_scale(K26CVector *v, double s);
double     k26c_vec_dot  (const K26CVector *a, const K26CVector *b);
double     k26c_vec_norm (const K26CVector *a);

/* ---- Matrices (row-major 2-D array of doubles) ---- */
typedef struct {
    double *data;          /* rows*cols entries, row-major */
    size_t  rows, cols;
} K26CMatrix;

K26CStatus k26c_mat_alloc(K26CMatrix *m, size_t rows, size_t cols);
void       k26c_mat_free (K26CMatrix *m);
K26CStatus k26c_mat_from (K26CMatrix *m, const double *src, size_t rows, size_t cols);
K26CStatus k26c_mat_set  (K26CMatrix *m, size_t r, size_t c, double x);
double     k26c_mat_get  (const K26CMatrix *m, size_t r, size_t c);
K26CStatus k26c_mat_identity (K26CMatrix *m);
K26CStatus k26c_mat_transpose(K26CMatrix *out, const K26CMatrix *in);
K26CStatus k26c_mat_mul      (K26CMatrix *out, const K26CMatrix *a, const K26CMatrix *b);
K26CStatus k26c_mat_vec      (K26CVector *out, const K26CMatrix *a, const K26CVector *x);

/* LU decomposition with partial pivoting. Singular matrices (any
 * pivot below 1e-300 in absolute value) return K26C_ERR_SINGULAR. */
K26CStatus k26c_mat_inverse(K26CMatrix *out, const K26CMatrix *in);
K26CStatus k26c_mat_det    (double *det_out, const K26CMatrix *in);
K26CStatus k26c_mat_solve  (K26CVector *x, const K26CMatrix *a, const K26CVector *b);

/* ---- Statistics ---- */
double k26c_stats_mean    (const double *xs, size_t n);
double k26c_stats_var     (const double *xs, size_t n, int sample);   /* sample=1 → /(n-1) */
double k26c_stats_std     (const double *xs, size_t n, int sample);
double k26c_stats_min     (const double *xs, size_t n);
double k26c_stats_max     (const double *xs, size_t n);
/* Quantile + median via in-place sort of `scratch[n]` (caller-owned). */
double k26c_stats_quantile(const double *xs, size_t n, double q, double *scratch);
double k26c_stats_median  (const double *xs, size_t n, double *scratch);
double k26c_stats_pearson (const double *xs, const double *ys, size_t n);

/* OLS linear regression y = a + b*x with uncertainty.
 *   r2          coefficient of determination
 *   a_se / b_se 1-σ standard errors of the intercept / slope
 *   n must be ≥ 2. */
typedef struct {
    double a, b;
    double a_se, b_se;
    double r2;
    size_t n;
} K26CLinReg;

K26CStatus k26c_stats_linreg(K26CLinReg *out, const double *xs,
                              const double *ys, size_t n);

/* χ² goodness-of-fit:  Σ (obs[i] - exp[i])² / exp[i]. */
double k26c_stats_chi2(const double *observed, const double *expected, size_t n);

/* ---- ODE solvers ---- */
typedef int (*K26CRhsFn)(double t, const K26CVector *y,
                         K26CVector *dydt, void *user);

K26CStatus k26c_ode_rk4 (K26CRhsFn rhs, void *user,
                         double t0, double t1, size_t n_steps,
                         K26CVector *y);
K26CStatus k26c_ode_rk45(K26CRhsFn rhs, void *user,
                         double t0, double t1,
                         double rtol, double atol,
                         K26CVector *y);

/* ---- Optimisation ---- */
typedef double (*K26CObj1Fn)(double x, void *user);
typedef double (*K26CObjNFn)(const K26CVector *x, void *user);

K26CStatus k26c_opt_brent       (K26CObj1Fn f, void *user,
                                  double lo, double hi, double tol,
                                  double *x_min, double *f_min);
K26CStatus k26c_opt_nelder_mead (K26CObjNFn f, void *user,
                                  const K26CVector *x0, double simplex_step,
                                  double tol, size_t max_iters,
                                  K26CVector *x_min, double *f_min);

/* ---- Deterministic RNG ----
 * SplitMix64 + Box-Muller normals. Tiny, fast, deterministic given
 * the seed. */
typedef struct {
    uint64_t s;
} K26CRng;

void   k26c_rng_init   (K26CRng *r, uint64_t seed);
double k26c_rng_uniform(K26CRng *r);                 /* [0, 1) */
double k26c_rng_normal (K26CRng *r, double mu, double sigma);

/* ---- Tagged result type ----
 * A K26CResult packages a computation's outcome (scalar, vector,
 * matrix, plot reference, or error) so a dispatcher can serialise
 * the value across a boundary without per-kind plumbing. */
typedef enum {
    K26C_RESULT_ERROR    = 0,
    K26C_RESULT_SCALAR   = 1,
    K26C_RESULT_VECTOR   = 2,
    K26C_RESULT_MATRIX   = 3,
    K26C_RESULT_PLOT_REF = 4,
} K26CResultKind;

typedef struct {
    K26CResultKind kind;
    union {
        char       *err_msg;          /* ERROR    — heap-owned */
        double      scalar;            /* SCALAR */
        K26CVector  vector;            /* VECTOR   — owns data */
        K26CMatrix  matrix;            /* MATRIX   — owns data */
        char       *plot_ref;          /* PLOT_REF — heap-owned (path / hash) */
    } u;
} K26CResult;

void k26c_result_init_error    (K26CResult *r, const char *msg);
void k26c_result_init_scalar   (K26CResult *r, double s);
/* The vector/matrix init steals ownership: caller MUST NOT free `v`
 * or `m` after handing it over (the result_free will). */
void k26c_result_init_vector   (K26CResult *r, K26CVector v);
void k26c_result_init_matrix   (K26CResult *r, K26CMatrix m);
void k26c_result_init_plot_ref (K26CResult *r, const char *ref);

void k26c_result_free(K26CResult *r);

#ifdef __cplusplus
}
#endif

#endif /* K26COMPUTE_H */
