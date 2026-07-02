/* mat4 — column-major 4x4 matrix algebra.
 *
 * Storage: m.m[col][row]. The four columns are the four basis vectors
 * the matrix maps the standard basis to. This matches GL/Vulkan
 * memory layout so a future GPU backend can memcpy through.
 *
 * The inverse routine is a textbook cofactor expansion (16 cofactors
 * plus a 4-term determinant); for 4x4 this is faster than LU + back
 * substitution and avoids any allocation. Determinism: no FMA, no
 * fused multiply-add reordering. */
#include "k26m3d.h"

#include <math.h>
#include <string.h>

void k26m3d_mat4_zero(K26M4 *out)
{
    if (!out) return;
    memset(out->m, 0, sizeof out->m);
}

void k26m3d_mat4_identity(K26M4 *out)
{
    if (!out) return;
    k26m3d_mat4_zero(out);
    out->m[0][0] = 1.0;
    out->m[1][1] = 1.0;
    out->m[2][2] = 1.0;
    out->m[3][3] = 1.0;
}

void k26m3d_mat4_copy(K26M4 *out, const K26M4 *in)
{
    if (!out || !in || out == in) return;
    memcpy(out->m, in->m, sizeof out->m);
}

void k26m3d_mat4_mul(K26M4 *out, const K26M4 *a, const K26M4 *b)
{
    if (!out || !a || !b) return;
    /* (out)_rc = sum_k (a)_rk * (b)_kc; columns of b dotted with rows
     * of a. Use a scratch result so out can alias a or b. */
    K26M4 t;
    for (int c = 0; c < 4; c++) {
        for (int r = 0; r < 4; r++) {
            double s = 0.0;
            for (int k = 0; k < 4; k++) {
                s += a->m[k][r] * b->m[c][k];
            }
            t.m[c][r] = s;
        }
    }
    *out = t;
}

void k26m3d_mat4_transpose(K26M4 *out, const K26M4 *a)
{
    if (!out || !a) return;
    K26M4 t;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            t.m[r][c] = a->m[c][r];
    *out = t;
}

K26V4 k26m3d_mat4_mul_v4(const K26M4 *a, K26V4 v)
{
    if (!a) { K26V4 z = { 0, 0, 0, 0 }; return z; }
    return k26m3d_v4(
        a->m[0][0]*v.x + a->m[1][0]*v.y + a->m[2][0]*v.z + a->m[3][0]*v.w,
        a->m[0][1]*v.x + a->m[1][1]*v.y + a->m[2][1]*v.z + a->m[3][1]*v.w,
        a->m[0][2]*v.x + a->m[1][2]*v.y + a->m[2][2]*v.z + a->m[3][2]*v.w,
        a->m[0][3]*v.x + a->m[1][3]*v.y + a->m[2][3]*v.z + a->m[3][3]*v.w
    );
}

K26V3 k26m3d_mat4_mul_point(const K26M4 *a, K26V3 p)
{
    K26V4 v = k26m3d_mat4_mul_v4(a, k26m3d_v4(p.x, p.y, p.z, 1.0));
    if (v.w != 0.0 && v.w != 1.0) {
        double iw = 1.0 / v.w;
        return k26m3d_v3(v.x * iw, v.y * iw, v.z * iw);
    }
    return k26m3d_v3(v.x, v.y, v.z);
}

K26V3 k26m3d_mat4_mul_dir(const K26M4 *a, K26V3 d)
{
    K26V4 v = k26m3d_mat4_mul_v4(a, k26m3d_v4(d.x, d.y, d.z, 0.0));
    return k26m3d_v3(v.x, v.y, v.z);
}

/* 4x4 cofactor expansion. Returns 0 on success; non-zero (-1) when
 * |det| < 1e-300. */
int k26m3d_mat4_inverse(K26M4 *out, const K26M4 *a)
{
    if (!out || !a) return -1;

    /* Flatten to row-major doubles for the cofactor formula; cleaner
     * than 16 separate locals. m[i] = a->m[col=i/4][row=i%4]? — no,
     * we want m[row*4 + col]. Use indexing that matches the textbook
     * cofactor table. */
    double m[16];
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            m[r * 4 + c] = a->m[c][r];

    double inv[16];

    inv[0]  =  m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
             + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
    inv[4]  = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
             - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
    inv[8]  =  m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
             + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
    inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
             - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];

    inv[1]  = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
             - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
    inv[5]  =  m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
             + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
    inv[9]  = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
             - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
    inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
             + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];

    inv[2]  =  m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
             + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
    inv[6]  = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
             - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
    inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
             + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
    inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
             - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];

    inv[3]  = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
             - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
    inv[7]  =  m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
             + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
    inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
             - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
    inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
             + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

    double det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
    if (fabs(det) < 1e-300) return -1;

    double inv_det = 1.0 / det;
    K26M4 t;
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            t.m[c][r] = inv[r * 4 + c] * inv_det;
    *out = t;
    return 0;
}
