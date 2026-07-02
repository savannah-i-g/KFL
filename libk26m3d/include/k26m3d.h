/* libk26m3d — fixed-size 3-D linear algebra for K26 rendering.
 *
 * This library provides vec2 / vec3 / vec4 / mat4 (column-major) and
 * quaternion types for double-precision 3-D graphics math. All types
 * are POD doubles, passed by value where they fit in two SSE registers
 * (V2/V3/V4) and by const-pointer for mat4 (128 bytes). Hot vec / vec
 * operations are `static inline` here; heavier mat / quat routines
 * live in libk26m3d.a.
 *
 * Convention:
 *   - Column-major mat4: m->m[col][row]. Matches GL/Vulkan memory
 *     layout so a GPU backend can memcpy through.
 *   - Right-handed coordinate system. +X right, +Y up, -Z forward.
 *   - Angles in radians at the C ABI (the KFL surface accepts degrees
 *     and converts at codegen).
 *
 * Determinism:
 *   - IEEE-754 strict; never depends on fast-math or FMA contraction.
 *   - The build flags in libk26m3d/Makefile pin this (-ffp-contract=off
 *     -fexcess-precision=standard). Downstream callers may keep their
 *     own flags; the library remains correct either way, only
 *     bit-exactness across compilers depends on it.
 *
 * Not provided: float32 vec / mat operations, SIMD intrinsics,
 * sparse / sparse-symmetric storage, eigendecomposition, dual
 * quaternions. */
#ifndef K26M3D_H
#define K26M3D_H

#include <math.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define K26M3D_LIB_VERSION "0.1.0"
#define K26M3D_PI          3.14159265358979323846

/* ----------------------------- Types ----------------------------- */

typedef struct { double x, y;             } K26V2;
typedef struct { double x, y, z;          } K26V3;
typedef struct { double x, y, z, w;       } K26V4;

/* Column-major 4x4. Element at row r, column c is m.m[c][r]. */
typedef struct { double m[4][4];          } K26M4;

/* Row-major 3x3. Element at row r, column c is m.m[r][c]. Row-major
 * matches the convention used by physics references (Markley &
 * Crassidis, Wertz, Hughes) for inertia tensors and rotation
 * matrices in body / inertial frames. Helpers for K26M3 (multiply,
 * invert) live with their primary consumers until enough call sites
 * accrue to warrant a libk26m3d-wide K26M3 algebra block. */
typedef struct { double m[3][3];          } K26M3;

/* Quaternion. w is the real part (scalar). */
typedef struct { double x, y, z, w;       } K26Quat;

/* ------------------------- Constructors -------------------------- */

static inline K26V2 k26m3d_v2(double x, double y)
{ K26V2 v = { x, y }; return v; }

static inline K26V3 k26m3d_v3(double x, double y, double z)
{ K26V3 v = { x, y, z }; return v; }

static inline K26V4 k26m3d_v4(double x, double y, double z, double w)
{ K26V4 v = { x, y, z, w }; return v; }

static inline K26V3 k26m3d_v3_from_v4(K26V4 v)
{ return k26m3d_v3(v.x, v.y, v.z); }

static inline K26V4 k26m3d_v4_from_v3(K26V3 v, double w)
{ return k26m3d_v4(v.x, v.y, v.z, w); }

/* --------------------------- V2 algebra -------------------------- */

static inline K26V2 k26m3d_v2_add(K26V2 a, K26V2 b)
{ return k26m3d_v2(a.x + b.x, a.y + b.y); }

static inline K26V2 k26m3d_v2_sub(K26V2 a, K26V2 b)
{ return k26m3d_v2(a.x - b.x, a.y - b.y); }

static inline K26V2 k26m3d_v2_scale(K26V2 a, double s)
{ return k26m3d_v2(a.x * s, a.y * s); }

static inline double k26m3d_v2_dot(K26V2 a, K26V2 b)
{ return a.x * b.x + a.y * b.y; }

static inline double k26m3d_v2_len2(K26V2 a) { return k26m3d_v2_dot(a, a); }
static inline double k26m3d_v2_len (K26V2 a) { return sqrt(k26m3d_v2_len2(a)); }

static inline K26V2 k26m3d_v2_norm(K26V2 a)
{
    double n = k26m3d_v2_len(a);
    return n > 0.0 ? k26m3d_v2_scale(a, 1.0 / n) : k26m3d_v2(0, 0);
}

/* --------------------------- V3 algebra -------------------------- */

static inline K26V3 k26m3d_v3_add(K26V3 a, K26V3 b)
{ return k26m3d_v3(a.x + b.x, a.y + b.y, a.z + b.z); }

static inline K26V3 k26m3d_v3_sub(K26V3 a, K26V3 b)
{ return k26m3d_v3(a.x - b.x, a.y - b.y, a.z - b.z); }

static inline K26V3 k26m3d_v3_neg(K26V3 a)
{ return k26m3d_v3(-a.x, -a.y, -a.z); }

static inline K26V3 k26m3d_v3_scale(K26V3 a, double s)
{ return k26m3d_v3(a.x * s, a.y * s, a.z * s); }

static inline K26V3 k26m3d_v3_mul(K26V3 a, K26V3 b)
{ return k26m3d_v3(a.x * b.x, a.y * b.y, a.z * b.z); }

static inline double k26m3d_v3_dot(K26V3 a, K26V3 b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

static inline K26V3 k26m3d_v3_cross(K26V3 a, K26V3 b)
{
    return k26m3d_v3(a.y * b.z - a.z * b.y,
                     a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x);
}

static inline double k26m3d_v3_len2(K26V3 a) { return k26m3d_v3_dot(a, a); }
static inline double k26m3d_v3_len (K26V3 a) { return sqrt(k26m3d_v3_len2(a)); }

static inline K26V3 k26m3d_v3_norm(K26V3 a)
{
    double n = k26m3d_v3_len(a);
    return n > 0.0 ? k26m3d_v3_scale(a, 1.0 / n) : k26m3d_v3(0, 0, 0);
}

static inline K26V3 k26m3d_v3_lerp(K26V3 a, K26V3 b, double t)
{ return k26m3d_v3_add(a, k26m3d_v3_scale(k26m3d_v3_sub(b, a), t)); }

/* --------------------------- V4 algebra -------------------------- */

K26V4  k26m3d_v4_add  (K26V4 a, K26V4 b);
K26V4  k26m3d_v4_sub  (K26V4 a, K26V4 b);
K26V4  k26m3d_v4_scale(K26V4 a, double s);
double k26m3d_v4_dot  (K26V4 a, K26V4 b);

/* --------------------------- M4 algebra -------------------------- */

void k26m3d_mat4_zero      (K26M4 *out);
void k26m3d_mat4_identity  (K26M4 *out);
void k26m3d_mat4_copy      (K26M4 *out, const K26M4 *in);
void k26m3d_mat4_mul       (K26M4 *out, const K26M4 *a, const K26M4 *b);
void k26m3d_mat4_transpose (K26M4 *out, const K26M4 *a);

/* Affine inverse — accepts general 4x4. Returns 0 on success, non-zero
 * if the matrix is singular (determinant below 1e-300 in absolute value).
 * `out` is left untouched on failure. */
int  k26m3d_mat4_inverse   (K26M4 *out, const K26M4 *a);

/* `out` = `a` * `v`. Treats v as a column. */
K26V4 k26m3d_mat4_mul_v4   (const K26M4 *a, K26V4 v);

/* Apply transform to a point (w=1) or direction (w=0). */
K26V3 k26m3d_mat4_mul_point  (const K26M4 *a, K26V3 p);
K26V3 k26m3d_mat4_mul_dir    (const K26M4 *a, K26V3 d);

/* ------------------------ Transform builders -------------------- */

void k26m3d_mat4_translate    (K26M4 *out, K26V3 t);
void k26m3d_mat4_scale_v      (K26M4 *out, K26V3 s);
void k26m3d_mat4_rotate_x     (K26M4 *out, double angle_rad);
void k26m3d_mat4_rotate_y     (K26M4 *out, double angle_rad);
void k26m3d_mat4_rotate_z     (K26M4 *out, double angle_rad);
void k26m3d_mat4_rotate_axis  (K26M4 *out, K26V3 axis_unit, double angle_rad);

/* Right-handed view matrix. Camera looks from eye toward center, with
 * +Y oriented along `up`. Output transforms world-space into
 * view-space where the camera is at the origin looking down -Z. */
void k26m3d_mat4_look_at(K26M4 *out, K26V3 eye, K26V3 center, K26V3 up);

/* Right-handed perspective. NDC z in [-1, 1] (GL convention). `fov_y`
 * in radians, `aspect` = width/height, both planes strictly positive
 * with `near < far`. */
void k26m3d_mat4_perspective(K26M4 *out,
                             double fov_y_rad, double aspect,
                             double near_plane, double far_plane);

/* Right-handed orthographic. NDC z in [-1, 1]. */
void k26m3d_mat4_orthographic(K26M4 *out,
                              double left, double right,
                              double bottom, double top,
                              double near_plane, double far_plane);

/* ----------------------------- Quat ----------------------------- */

static inline K26Quat k26m3d_quat(double x, double y, double z, double w)
{ K26Quat q = { x, y, z, w }; return q; }

static inline K26Quat k26m3d_quat_identity(void)
{ return k26m3d_quat(0, 0, 0, 1); }

K26Quat k26m3d_quat_from_axis_angle(K26V3 axis_unit, double angle_rad);
K26Quat k26m3d_quat_mul            (K26Quat a, K26Quat b);
K26Quat k26m3d_quat_conj           (K26Quat q);
K26Quat k26m3d_quat_norm           (K26Quat q);
K26Quat k26m3d_quat_slerp          (K26Quat a, K26Quat b, double t);
void    k26m3d_quat_to_mat4        (K26M4 *out, K26Quat q);
K26V3   k26m3d_quat_rotate_v3      (K26Quat q, K26V3 v);

/* ---------------------------- Helpers --------------------------- */

static inline double k26m3d_deg2rad(double d) { return d * (K26M3D_PI / 180.0); }
static inline double k26m3d_rad2deg(double r) { return r * (180.0 / K26M3D_PI); }
static inline double k26m3d_clamp  (double x, double lo, double hi)
{ return x < lo ? lo : (x > hi ? hi : x); }
static inline double k26m3d_lerp   (double a, double b, double t)
{ return a + (b - a) * t; }

#ifdef __cplusplus
}
#endif

#endif /* K26M3D_H */
