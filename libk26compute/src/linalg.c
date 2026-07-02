/* libk26compute — LU-based linear algebra (inverse, determinant,
 * solve). Single in-place LU with partial pivoting drives all three
 * entry points; downstream callers benefit from a known-good
 * factorisation rather than re-implementing per operation. */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Pivot threshold — any pivot smaller than this in absolute value
 * means the matrix is singular (or near enough that inversion would
 * blow up). 1e-300 catches actual singular cases without rejecting
 * legitimately small pivots that arise from well-conditioned
 * matrices with extreme scale. */
#define K26C_LU_PIVOT_EPS 1e-300

/* In-place LU with partial pivoting.
 *   A[]    n×n matrix, overwritten with combined L (unit diagonal,
 *          stored below the diagonal) + U (on and above the diagonal).
 *   perm[] n-element row-permutation record.
 *   parity ±1; multiplied into the determinant by the caller. */
static K26CStatus lu_decompose(double *A, size_t n,
                               size_t *perm, int *parity)
{
    *parity = 1;
    for (size_t i = 0; i < n; i++) perm[i] = i;

    for (size_t k = 0; k < n; k++) {
        /* Find the largest-magnitude pivot in column k, rows k..n-1. */
        size_t piv  = k;
        double maxv = fabs(A[k * n + k]);
        for (size_t i = k + 1; i < n; i++) {
            double v = fabs(A[i * n + k]);
            if (v > maxv) { maxv = v; piv = i; }
        }
        if (maxv < K26C_LU_PIVOT_EPS) return K26C_ERR_SINGULAR;

        if (piv != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp        = A[k * n + j];
                A[k * n + j]      = A[piv * n + j];
                A[piv * n + j]    = tmp;
            }
            size_t pt = perm[k];
            perm[k]   = perm[piv];
            perm[piv] = pt;
            *parity   = -*parity;
        }

        double pivval = A[k * n + k];
        for (size_t i = k + 1; i < n; i++) {
            A[i * n + k] /= pivval;
            double f = A[i * n + k];
            for (size_t j = k + 1; j < n; j++) {
                A[i * n + j] -= f * A[k * n + j];
            }
        }
    }
    return K26C_OK;
}

/* Solve LU·x = P·b, in-place: b on entry holds the RHS, on return
 * holds x. Returns OOM if scratch alloc fails. */
static K26CStatus lu_solve(const double *LU, size_t n,
                           const size_t *perm, double *b)
{
    double *tmp = (double *)malloc(n * sizeof(double));
    if (!tmp) return K26C_ERR_OOM;

    /* Apply permutation: tmp = P·b. */
    for (size_t i = 0; i < n; i++) tmp[i] = b[perm[i]];

    /* Forward substitution: solve L·y = P·b (L has unit diagonal). */
    for (size_t i = 0; i < n; i++) {
        double s = tmp[i];
        for (size_t j = 0; j < i; j++) s -= LU[i * n + j] * tmp[j];
        tmp[i] = s;
    }

    /* Back substitution: solve U·x = y. */
    for (size_t ii = n; ii > 0; ii--) {
        size_t i = ii - 1;
        double s = tmp[i];
        for (size_t j = i + 1; j < n; j++) s -= LU[i * n + j] * tmp[j];
        tmp[i] = s / LU[i * n + i];
    }

    memcpy(b, tmp, n * sizeof(double));
    free(tmp);
    return K26C_OK;
}

K26CStatus k26c_mat_inverse(K26CMatrix *out, const K26CMatrix *in)
{
    if (!out || !in || in->rows != in->cols || in->rows == 0)
        return K26C_ERR_INVAL;
    size_t n = in->rows;

    double *LU   = (double *)malloc(n * n * sizeof(double));
    size_t *perm = (size_t *)malloc(n     * sizeof(size_t));
    if (!LU || !perm) { free(LU); free(perm); return K26C_ERR_OOM; }
    memcpy(LU, in->data, n * n * sizeof(double));

    int parity;
    K26CStatus rc = lu_decompose(LU, n, perm, &parity);
    if (rc != K26C_OK) { free(LU); free(perm); return rc; }

    K26CMatrix t;
    rc = k26c_mat_alloc(&t, n, n);
    if (rc != K26C_OK) { free(LU); free(perm); return rc; }

    /* Solve LU·x_c = e_c for each unit-vector e_c; columns of x form
     * the inverse. */
    double *col = (double *)malloc(n * sizeof(double));
    if (!col) {
        k26c_mat_free(&t); free(LU); free(perm);
        return K26C_ERR_OOM;
    }
    for (size_t c = 0; c < n; c++) {
        for (size_t i = 0; i < n; i++) col[i] = (i == c) ? 1.0 : 0.0;
        K26CStatus sr = lu_solve(LU, n, perm, col);
        if (sr != K26C_OK) {
            free(col); k26c_mat_free(&t); free(LU); free(perm);
            return sr;
        }
        for (size_t i = 0; i < n; i++) t.data[i * n + c] = col[i];
    }

    free(col); free(LU); free(perm);
    k26c_mat_free(out);
    *out = t;
    return K26C_OK;
}

K26CStatus k26c_mat_det(double *det_out, const K26CMatrix *in)
{
    if (!det_out || !in || in->rows != in->cols) return K26C_ERR_INVAL;
    size_t n = in->rows;
    if (n == 0) { *det_out = 1.0; return K26C_OK; }

    double *LU   = (double *)malloc(n * n * sizeof(double));
    size_t *perm = (size_t *)malloc(n     * sizeof(size_t));
    if (!LU || !perm) { free(LU); free(perm); return K26C_ERR_OOM; }
    memcpy(LU, in->data, n * n * sizeof(double));

    int parity;
    K26CStatus rc = lu_decompose(LU, n, perm, &parity);
    if (rc != K26C_OK) {
        /* A singular matrix has det == 0; surface that explicitly so
         * callers don't have to special-case the status. */
        free(LU); free(perm);
        if (rc == K26C_ERR_SINGULAR) { *det_out = 0.0; return K26C_OK; }
        return rc;
    }

    double det = (double)parity;
    for (size_t i = 0; i < n; i++) det *= LU[i * n + i];
    free(LU); free(perm);
    *det_out = det;
    return K26C_OK;
}

K26CStatus k26c_mat_solve(K26CVector *x, const K26CMatrix *a, const K26CVector *b)
{
    if (!x || !a || !b || a->rows != a->cols || a->rows != b->n)
        return K26C_ERR_INVAL;
    size_t n = a->rows;
    if (n == 0) {
        k26c_vec_free(x);
        return k26c_vec_alloc(x, 0);
    }

    double *LU   = (double *)malloc(n * n * sizeof(double));
    size_t *perm = (size_t *)malloc(n     * sizeof(size_t));
    if (!LU || !perm) { free(LU); free(perm); return K26C_ERR_OOM; }
    memcpy(LU, a->data, n * n * sizeof(double));

    int parity;
    K26CStatus rc = lu_decompose(LU, n, perm, &parity);
    if (rc != K26C_OK) { free(LU); free(perm); return rc; }

    K26CVector t;
    rc = k26c_vec_alloc(&t, n);
    if (rc != K26C_OK) { free(LU); free(perm); return rc; }
    memcpy(t.data, b->data, n * sizeof(double));

    rc = lu_solve(LU, n, perm, t.data);
    free(LU); free(perm);
    if (rc != K26C_OK) { k26c_vec_free(&t); return rc; }

    k26c_vec_free(x);
    *x = t;
    return K26C_OK;
}
