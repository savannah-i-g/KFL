/* libk26compute — vector primitives (heap-owned 1-D doubles). */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

K26CStatus k26c_vec_alloc(K26CVector *v, size_t n)
{
    if (!v) return K26C_ERR_INVAL;
    v->data = NULL;
    v->n    = 0;
    if (n == 0) return K26C_OK;
    v->data = (double *)calloc(n, sizeof(double));
    if (!v->data) return K26C_ERR_OOM;
    v->n = n;
    return K26C_OK;
}

void k26c_vec_free(K26CVector *v)
{
    if (!v) return;
    free(v->data);
    v->data = NULL;
    v->n    = 0;
}

K26CStatus k26c_vec_from(K26CVector *v, const double *src, size_t n)
{
    if (!v || (n > 0 && !src)) return K26C_ERR_INVAL;
    K26CStatus rc = k26c_vec_alloc(v, n);
    if (rc != K26C_OK) return rc;
    if (n > 0) memcpy(v->data, src, n * sizeof(double));
    return K26C_OK;
}

K26CStatus k26c_vec_set(K26CVector *v, size_t i, double x)
{
    if (!v || !v->data || i >= v->n) return K26C_ERR_INVAL;
    v->data[i] = x;
    return K26C_OK;
}

double k26c_vec_get(const K26CVector *v, size_t i)
{
    if (!v || !v->data || i >= v->n) return 0.0;
    return v->data[i];
}

K26CStatus k26c_vec_fill(K26CVector *v, double x)
{
    if (!v) return K26C_ERR_INVAL;
    for (size_t i = 0; i < v->n; i++) v->data[i] = x;
    return K26C_OK;
}

K26CStatus k26c_vec_copy(K26CVector *dst, const K26CVector *src)
{
    if (!dst || !src) return K26C_ERR_INVAL;
    if (dst == src) return K26C_OK;
    K26CVector t;
    K26CStatus rc = k26c_vec_alloc(&t, src->n);
    if (rc != K26C_OK) return rc;
    if (src->n > 0) memcpy(t.data, src->data, src->n * sizeof(double));
    k26c_vec_free(dst);
    *dst = t;
    return K26C_OK;
}

/* Helper: ensure `out` has length n; resize if needed. */
static K26CStatus vec_ensure(K26CVector *out, size_t n)
{
    if (out->n == n && (n == 0 || out->data)) return K26C_OK;
    k26c_vec_free(out);
    return k26c_vec_alloc(out, n);
}

K26CStatus k26c_vec_add(K26CVector *out, const K26CVector *a, const K26CVector *b)
{
    if (!out || !a || !b || a->n != b->n) return K26C_ERR_INVAL;
    K26CStatus rc = vec_ensure(out, a->n);
    if (rc != K26C_OK) return rc;
    for (size_t i = 0; i < a->n; i++) out->data[i] = a->data[i] + b->data[i];
    return K26C_OK;
}

K26CStatus k26c_vec_sub(K26CVector *out, const K26CVector *a, const K26CVector *b)
{
    if (!out || !a || !b || a->n != b->n) return K26C_ERR_INVAL;
    K26CStatus rc = vec_ensure(out, a->n);
    if (rc != K26C_OK) return rc;
    for (size_t i = 0; i < a->n; i++) out->data[i] = a->data[i] - b->data[i];
    return K26C_OK;
}

K26CStatus k26c_vec_scale(K26CVector *v, double s)
{
    if (!v) return K26C_ERR_INVAL;
    for (size_t i = 0; i < v->n; i++) v->data[i] *= s;
    return K26C_OK;
}

double k26c_vec_dot(const K26CVector *a, const K26CVector *b)
{
    if (!a || !b || a->n != b->n) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i < a->n; i++) s += a->data[i] * b->data[i];
    return s;
}

double k26c_vec_norm(const K26CVector *a)
{
    if (!a) return 0.0;
    /* Two-pass scaling to avoid overflow on huge components — same
     * trick BLAS dnrm2 uses. */
    double scale = 0.0;
    for (size_t i = 0; i < a->n; i++) {
        double v = fabs(a->data[i]);
        if (v > scale) scale = v;
    }
    if (scale == 0.0) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i < a->n; i++) {
        double r = a->data[i] / scale;
        s += r * r;
    }
    return scale * sqrt(s);
}
