/* libk26compute — K26CResult tagged-union helpers.
 *
 * Packages a computation's outcome (scalar, vector, matrix, plot
 * reference, or error) so callers can hand a single value back across
 * a serialisation boundary without per-kind plumbing. */

#include "k26compute.h"

#include <stdlib.h>
#include <string.h>

static char *xstrdup_or_null(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *r = (char *)malloc(n);
    if (r) memcpy(r, s, n);
    return r;
}

void k26c_result_init_error(K26CResult *r, const char *msg)
{
    if (!r) return;
    r->kind = K26C_RESULT_ERROR;
    r->u.err_msg = xstrdup_or_null(msg ? msg : "unknown error");
}

void k26c_result_init_scalar(K26CResult *r, double s)
{
    if (!r) return;
    r->kind = K26C_RESULT_SCALAR;
    r->u.scalar = s;
}

void k26c_result_init_vector(K26CResult *r, K26CVector v)
{
    if (!r) return;
    r->kind = K26C_RESULT_VECTOR;
    r->u.vector = v;        /* steals ownership */
}

void k26c_result_init_matrix(K26CResult *r, K26CMatrix m)
{
    if (!r) return;
    r->kind = K26C_RESULT_MATRIX;
    r->u.matrix = m;        /* steals ownership */
}

void k26c_result_init_plot_ref(K26CResult *r, const char *ref)
{
    if (!r) return;
    r->kind = K26C_RESULT_PLOT_REF;
    r->u.plot_ref = xstrdup_or_null(ref ? ref : "");
}

void k26c_result_free(K26CResult *r)
{
    if (!r) return;
    switch (r->kind) {
    case K26C_RESULT_ERROR:
        free(r->u.err_msg);
        r->u.err_msg = NULL;
        break;
    case K26C_RESULT_VECTOR:
        k26c_vec_free(&r->u.vector);
        break;
    case K26C_RESULT_MATRIX:
        k26c_mat_free(&r->u.matrix);
        break;
    case K26C_RESULT_PLOT_REF:
        free(r->u.plot_ref);
        r->u.plot_ref = NULL;
        break;
    case K26C_RESULT_SCALAR:
        /* nothing to free */
        break;
    }
    r->kind = K26C_RESULT_ERROR;       /* leave in a defined-empty state */
    r->u.err_msg = NULL;
}
