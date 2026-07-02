/* libk26m3d smoke test. Verifies the load-bearing identities the
 * renderer relies on: matrix inverse round-trip, look_at maps the eye
 * to the origin, perspective sends the near plane to z=-1 and far to
 * z=+1, quat slerp endpoints, axis-angle round-trip via to_mat4.
 *
 * Prints PASS/FAIL per test; exits 0 if all pass, 1 otherwise. */
#include "k26m3d.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int g_failures = 0;

static void check(const char *name, int ok, const char *detail)
{
    if (ok) {
        printf("PASS  %s  (%s)\n", name, detail);
    } else {
        printf("FAIL  %s  (%s)\n", name, detail);
        g_failures++;
    }
}

static int approx_eq(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}

static int v3_eq(K26V3 a, K26V3 b, double tol)
{
    return approx_eq(a.x, b.x, tol) && approx_eq(a.y, b.y, tol)
        && approx_eq(a.z, b.z, tol);
}

/* ---- 1. v3 cross product / dot identities. */
static void test_v3_basics(void)
{
    K26V3 x = k26m3d_v3(1, 0, 0);
    K26V3 y = k26m3d_v3(0, 1, 0);
    K26V3 z = k26m3d_v3(0, 0, 1);

    K26V3 xy = k26m3d_v3_cross(x, y);
    K26V3 yz = k26m3d_v3_cross(y, z);
    K26V3 zx = k26m3d_v3_cross(z, x);

    int ok = v3_eq(xy, z, 1e-15)
          && v3_eq(yz, x, 1e-15)
          && v3_eq(zx, y, 1e-15)
          && approx_eq(k26m3d_v3_dot(x, x), 1.0, 1e-15)
          && approx_eq(k26m3d_v3_len(k26m3d_v3(3, 4, 0)), 5.0, 1e-15);
    check("v3 cross+dot+len", ok, "xyz, dot, len");
}

/* ---- 2. mat4 identity and identity * M == M. */
static void test_mat4_identity(void)
{
    K26M4 I, M, P;
    k26m3d_mat4_identity(&I);
    /* M is a non-trivial mat4. */
    k26m3d_mat4_zero(&M);
    M.m[0][0] = 2;  M.m[1][0] = 3;  M.m[2][0] = 5;  M.m[3][0] = 7;
    M.m[0][1] = 11; M.m[1][1] = 13; M.m[2][1] = 17; M.m[3][1] = 19;
    M.m[0][2] = 23; M.m[1][2] = 29; M.m[2][2] = 31; M.m[3][2] = 37;
    M.m[0][3] = 41; M.m[1][3] = 43; M.m[2][3] = 47; M.m[3][3] = 53;

    k26m3d_mat4_mul(&P, &I, &M);
    double err = 0.0;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            err += fabs(P.m[c][r] - M.m[c][r]);
    check("mat4 identity*M == M", err < 1e-12, "elementwise abs sum");
}

/* ---- 3. mat4 inverse round-trip: M * M^-1 == I (within tol). */
static void test_mat4_inverse(void)
{
    K26M4 M, Mi, P, I;
    k26m3d_mat4_identity(&I);
    k26m3d_mat4_zero(&M);
    /* A well-conditioned matrix: rotate around (1,1,1) by 1.3 rad
     * then translate by (5, -3, 2). */
    K26M4 R, T;
    k26m3d_mat4_rotate_axis(&R, k26m3d_v3(1, 1, 1), 1.3);
    k26m3d_mat4_translate(&T, k26m3d_v3(5, -3, 2));
    k26m3d_mat4_mul(&M, &T, &R);

    int rc = k26m3d_mat4_inverse(&Mi, &M);
    if (rc != 0) {
        check("mat4 inverse", 0, "k26m3d_mat4_inverse failed");
        return;
    }
    k26m3d_mat4_mul(&P, &M, &Mi);
    double err = 0.0;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            err += fabs(P.m[c][r] - I.m[c][r]);

    char detail[64];
    snprintf(detail, sizeof detail, "||M.M^-1 - I||_1 = %.3e", err);
    check("mat4 inverse round-trip", err < 1e-10, detail);
}

/* ---- 4. look_at maps eye to origin, center to negative-z axis. */
static void test_look_at(void)
{
    K26V3 eye    = k26m3d_v3(3, 4, 5);
    K26V3 center = k26m3d_v3(0, 0, 0);
    K26V3 up     = k26m3d_v3(0, 1, 0);
    K26M4 V;
    k26m3d_mat4_look_at(&V, eye, center, up);

    /* The eye in view-space should be the origin. */
    K26V3 eye_v = k26m3d_mat4_mul_point(&V, eye);
    int ok_eye = approx_eq(eye_v.x, 0.0, 1e-12)
              && approx_eq(eye_v.y, 0.0, 1e-12)
              && approx_eq(eye_v.z, 0.0, 1e-12);
    char detail[96];
    snprintf(detail, sizeof detail, "eye_v=(%.3e,%.3e,%.3e)",
             eye_v.x, eye_v.y, eye_v.z);
    check("look_at: eye -> origin", ok_eye, detail);

    /* The center in view-space should lie along -Z. */
    K26V3 c_v = k26m3d_mat4_mul_point(&V, center);
    int ok_c = approx_eq(c_v.x, 0.0, 1e-12)
            && approx_eq(c_v.y, 0.0, 1e-12)
            && c_v.z < 0.0;
    snprintf(detail, sizeof detail, "center_v=(%.3e,%.3e,%.3e)",
             c_v.x, c_v.y, c_v.z);
    check("look_at: center along -Z", ok_c, detail);
}

/* ---- 5. perspective: near plane → ndc z = -1, far → +1 (GL). */
static void test_perspective(void)
{
    K26M4 P;
    double near_plane = 0.5, far_plane = 50.0;
    k26m3d_mat4_perspective(&P, k26m3d_deg2rad(60.0), 16.0/9.0,
                            near_plane, far_plane);

    /* A point at view-space (0,0,-near). After perspective, clip-space
     * w should be near. NDC z = clip.z / clip.w should be -1.0. */
    K26V4 near_pt = k26m3d_mat4_mul_v4(&P, k26m3d_v4(0, 0, -near_plane, 1));
    K26V4 far_pt  = k26m3d_mat4_mul_v4(&P, k26m3d_v4(0, 0, -far_plane,  1));
    double ndc_z_near = near_pt.z / near_pt.w;
    double ndc_z_far  = far_pt.z  / far_pt.w;

    char detail[96];
    snprintf(detail, sizeof detail,
             "ndc_z(near)=%.4f ndc_z(far)=%.4f", ndc_z_near, ndc_z_far);
    int ok = approx_eq(ndc_z_near, -1.0, 1e-12)
          && approx_eq(ndc_z_far,   1.0, 1e-12);
    check("perspective: near→-1, far→+1", ok, detail);
}

/* ---- 6. Quat slerp endpoints. */
static void test_quat_slerp(void)
{
    K26Quat a = k26m3d_quat_from_axis_angle(k26m3d_v3(0,1,0), 0.0);
    K26Quat b = k26m3d_quat_from_axis_angle(k26m3d_v3(0,1,0), K26M3D_PI / 2);

    K26Quat s0 = k26m3d_quat_slerp(a, b, 0.0);
    K26Quat s1 = k26m3d_quat_slerp(a, b, 1.0);

    int ok0 = approx_eq(s0.x, a.x, 1e-12) && approx_eq(s0.y, a.y, 1e-12)
           && approx_eq(s0.z, a.z, 1e-12) && approx_eq(s0.w, a.w, 1e-12);
    int ok1 = approx_eq(s1.x, b.x, 1e-12) && approx_eq(s1.y, b.y, 1e-12)
           && approx_eq(s1.z, b.z, 1e-12) && approx_eq(s1.w, b.w, 1e-12);
    check("quat slerp endpoints", ok0 && ok1, "t=0 and t=1");
}

/* ---- 7. Axis-angle → quat → mat4 → vector ≡ rotate_axis(mat4) applied. */
static void test_quat_mat4(void)
{
    K26V3 axis = k26m3d_v3_norm(k26m3d_v3(1, 2, 3));
    double a = 0.7;
    K26Quat q = k26m3d_quat_from_axis_angle(axis, a);
    K26M4 R_quat, R_axis;
    k26m3d_quat_to_mat4(&R_quat, q);
    k26m3d_mat4_rotate_axis(&R_axis, axis, a);

    K26V3 v = k26m3d_v3(0.5, -0.25, 1.0);
    K26V3 rq = k26m3d_mat4_mul_dir(&R_quat, v);
    K26V3 ra = k26m3d_mat4_mul_dir(&R_axis, v);
    int ok = v3_eq(rq, ra, 1e-12);
    char detail[160];
    snprintf(detail, sizeof detail,
             "rq=(%.4f,%.4f,%.4f) ra=(%.4f,%.4f,%.4f)",
             rq.x, rq.y, rq.z, ra.x, ra.y, ra.z);
    check("quat→mat4 == rotate_axis", ok, detail);
}

int main(void)
{
    test_v3_basics();
    test_mat4_identity();
    test_mat4_inverse();
    test_look_at();
    test_perspective();
    test_quat_slerp();
    test_quat_mat4();

    printf("\n%s — %d failure(s)\n",
           g_failures == 0 ? "OK" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
