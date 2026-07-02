/* quat — unit quaternion algebra for K26 rotations.
 *
 * Convention: q = (x, y, z, w) with w the real (scalar) part.
 * Quaternion multiplication is non-commutative; (q*r)(v) applies r
 * first then q. */
#include "k26m3d.h"

#include <math.h>

K26Quat k26m3d_quat_from_axis_angle(K26V3 axis, double angle)
{
    K26V3 u = k26m3d_v3_norm(axis);
    double s = sin(angle * 0.5);
    return k26m3d_quat(u.x * s, u.y * s, u.z * s, cos(angle * 0.5));
}

K26Quat k26m3d_quat_mul(K26Quat a, K26Quat b)
{
    return k26m3d_quat(
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    );
}

K26Quat k26m3d_quat_conj(K26Quat q)
{
    return k26m3d_quat(-q.x, -q.y, -q.z, q.w);
}

K26Quat k26m3d_quat_norm(K26Quat q)
{
    double n2 = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
    if (n2 <= 0.0) return k26m3d_quat_identity();
    double inv = 1.0 / sqrt(n2);
    return k26m3d_quat(q.x * inv, q.y * inv, q.z * inv, q.w * inv);
}

K26Quat k26m3d_quat_slerp(K26Quat a, K26Quat b, double t)
{
    /* Cosine of angle between unit quaternions. Negate b if it
     * crossed the hemisphere — interpolating the short arc is the
     * convention. */
    double cos_h = a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
    if (cos_h < 0.0) {
        b.x = -b.x; b.y = -b.y; b.z = -b.z; b.w = -b.w;
        cos_h = -cos_h;
    }
    if (cos_h > 0.9995) {
        /* Quaternions very close — linear interpolate + normalize to
         * avoid sin(theta) → 0. */
        K26Quat r = k26m3d_quat(
            a.x + (b.x - a.x) * t,
            a.y + (b.y - a.y) * t,
            a.z + (b.z - a.z) * t,
            a.w + (b.w - a.w) * t
        );
        return k26m3d_quat_norm(r);
    }
    double theta_0  = acos(cos_h);
    double theta    = theta_0 * t;
    double sin_t    = sin(theta);
    double sin_t0   = sin(theta_0);
    double s0 = cos(theta) - cos_h * sin_t / sin_t0;
    double s1 = sin_t / sin_t0;
    return k26m3d_quat(
        s0 * a.x + s1 * b.x,
        s0 * a.y + s1 * b.y,
        s0 * a.z + s1 * b.z,
        s0 * a.w + s1 * b.w
    );
}

void k26m3d_quat_to_mat4(K26M4 *out, K26Quat q)
{
    K26Quat n = k26m3d_quat_norm(q);
    double x = n.x, y = n.y, z = n.z, w = n.w;
    double xx = x*x, yy = y*y, zz = z*z;
    double xy = x*y, xz = x*z, yz = y*z;
    double wx = w*x, wy = w*y, wz = w*z;

    k26m3d_mat4_identity(out);
    out->m[0][0] = 1.0 - 2.0*(yy + zz);
    out->m[0][1] =        2.0*(xy + wz);
    out->m[0][2] =        2.0*(xz - wy);
    out->m[1][0] =        2.0*(xy - wz);
    out->m[1][1] = 1.0 - 2.0*(xx + zz);
    out->m[1][2] =        2.0*(yz + wx);
    out->m[2][0] =        2.0*(xz + wy);
    out->m[2][1] =        2.0*(yz - wx);
    out->m[2][2] = 1.0 - 2.0*(xx + yy);
}

K26V3 k26m3d_quat_rotate_v3(K26Quat q, K26V3 v)
{
    /* v' = q * (v, 0) * q^-1, but with the conjugate trick. */
    K26V3 u = k26m3d_v3(q.x, q.y, q.z);
    double s = q.w;
    K26V3 t  = k26m3d_v3_scale(k26m3d_v3_cross(u, v), 2.0);
    K26V3 a  = k26m3d_v3_add(v, k26m3d_v3_scale(t, s));
    return k26m3d_v3_add(a, k26m3d_v3_cross(u, t));
}
