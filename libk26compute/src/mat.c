/* libk26compute — matrix primitives (row-major, heap-owned). */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

K26CStatus k26c_mat_alloc(K26CMatrix *m, size_t rows, size_t cols)
{
    if (!m) return K26C_ERR_INVAL;
    m->data = NULL;
    m->rows = m->cols = 0;
    if (rows == 0 || cols == 0) return K26C_OK;
    m->data = (double *)calloc(rows * cols, sizeof(double));
    if (!m->data) return K26C_ERR_OOM;
    m->rows = rows;
    m->cols = cols;
    return K26C_OK;
}

void k26c_mat_free(K26CMatrix *m)
{
    if (!m) return;
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

K26CStatus k26c_mat_from(K26CMatrix *m, const double *src, size_t rows, size_t cols)
{
    if (!m || (rows*cols > 0 && !src)) return K26C_ERR_INVAL;
    K26CStatus rc = k26c_mat_alloc(m, rows, cols);
    if (rc != K26C_OK) return rc;
    if (rows*cols > 0) memcpy(m->data, src, rows * cols * sizeof(double));
    return K26C_OK;
}

K26CStatus k26c_mat_set(K26CMatrix *m, size_t r, size_t c, double x)
{
    if (!m || !m->data || r >= m->rows || c >= m->cols) return K26C_ERR_INVAL;
    m->data[r * m->cols + c] = x;
    return K26C_OK;
}

double k26c_mat_get(const K26CMatrix *m, size_t r, size_t c)
{
    if (!m || !m->data || r >= m->rows || c >= m->cols) return 0.0;
    return m->data[r * m->cols + c];
}

K26CStatus k26c_mat_identity(K26CMatrix *m)
{
    if (!m || m->rows != m->cols) return K26C_ERR_INVAL;
    memset(m->data, 0, m->rows * m->cols * sizeof(double));
    for (size_t i = 0; i < m->rows; i++) m->data[i * m->cols + i] = 1.0;
    return K26C_OK;
}

K26CStatus k26c_mat_transpose(K26CMatrix *out, const K26CMatrix *in)
{
    if (!out || !in) return K26C_ERR_INVAL;
    K26CMatrix t;
    K26CStatus rc = k26c_mat_alloc(&t, in->cols, in->rows);
    if (rc != K26C_OK) return rc;
    for (size_t r = 0; r < in->rows; r++) {
        for (size_t c = 0; c < in->cols; c++) {
            t.data[c * t.cols + r] = in->data[r * in->cols + c];
        }
    }
    k26c_mat_free(out);
    *out = t;
    return K26C_OK;
}

K26CStatus k26c_mat_mul(K26CMatrix *out, const K26CMatrix *a, const K26CMatrix *b)
{
    if (!out || !a || !b || a->cols != b->rows) return K26C_ERR_INVAL;
    K26CMatrix t;
    K26CStatus rc = k26c_mat_alloc(&t, a->rows, b->cols);
    if (rc != K26C_OK) return rc;
    /* Naive triple-loop. Adequate for small / medium matrices; if
     * we ever wire BLAS in, this is the seam to swap. */
    for (size_t r = 0; r < a->rows; r++) {
        for (size_t c = 0; c < b->cols; c++) {
            double s = 0.0;
            for (size_t k = 0; k < a->cols; k++) {
                s += a->data[r * a->cols + k] * b->data[k * b->cols + c];
            }
            t.data[r * t.cols + c] = s;
        }
    }
    k26c_mat_free(out);
    *out = t;
    return K26C_OK;
}

K26CStatus k26c_mat_vec(K26CVector *out, const K26CMatrix *a, const K26CVector *x)
{
    if (!out || !a || !x || a->cols != x->n) return K26C_ERR_INVAL;
    K26CVector t;
    K26CStatus rc = k26c_vec_alloc(&t, a->rows);
    if (rc != K26C_OK) return rc;
    for (size_t r = 0; r < a->rows; r++) {
        double s = 0.0;
        for (size_t c = 0; c < a->cols; c++) {
            s += a->data[r * a->cols + c] * x->data[c];
        }
        t.data[r] = s;
    }
    k26c_vec_free(out);
    *out = t;
    return K26C_OK;
}
