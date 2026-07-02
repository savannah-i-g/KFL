/* transform — affine + camera matrix builders.
 *
 * Right-handed, column-major. Conventions:
 *   - look_at: camera at `eye` looking at `center`, +Y aligned with
 *     `up`. View matrix maps world to camera space (camera at origin
 *     looking down -Z, +X right, +Y up).
 *   - perspective: NDC z in [-1, 1] (GL convention). fov_y in radians.
 *   - orthographic: NDC z in [-1, 1]. left/right/bottom/top in scene
 *     space. */
#include "k26m3d.h"

#include <math.h>

void k26m3d_mat4_translate(K26M4 *out, K26V3 t)
{
    k26m3d_mat4_identity(out);
    out->m[3][0] = t.x;
    out->m[3][1] = t.y;
    out->m[3][2] = t.z;
}

void k26m3d_mat4_scale_v(K26M4 *out, K26V3 s)
{
    k26m3d_mat4_zero(out);
    out->m[0][0] = s.x;
    out->m[1][1] = s.y;
    out->m[2][2] = s.z;
    out->m[3][3] = 1.0;
}

void k26m3d_mat4_rotate_x(K26M4 *out, double a)
{
    double c = cos(a), s = sin(a);
    k26m3d_mat4_identity(out);
    out->m[1][1] =  c;   out->m[2][1] = -s;
    out->m[1][2] =  s;   out->m[2][2] =  c;
}

void k26m3d_mat4_rotate_y(K26M4 *out, double a)
{
    double c = cos(a), s = sin(a);
    k26m3d_mat4_identity(out);
    out->m[0][0] =  c;   out->m[2][0] =  s;
    out->m[0][2] = -s;   out->m[2][2] =  c;
}

void k26m3d_mat4_rotate_z(K26M4 *out, double a)
{
    double c = cos(a), s = sin(a);
    k26m3d_mat4_identity(out);
    out->m[0][0] =  c;   out->m[1][0] = -s;
    out->m[0][1] =  s;   out->m[1][1] =  c;
}

void k26m3d_mat4_rotate_axis(K26M4 *out, K26V3 axis, double a)
{
    K26V3 u = k26m3d_v3_norm(axis);
    double c = cos(a);
    double s = sin(a);
    double t = 1.0 - c;
    double x = u.x, y = u.y, z = u.z;

    /* Rodrigues' rotation formula, column-major layout. */
    k26m3d_mat4_identity(out);
    out->m[0][0] = t*x*x + c;     out->m[1][0] = t*x*y - s*z;   out->m[2][0] = t*x*z + s*y;
    out->m[0][1] = t*x*y + s*z;   out->m[1][1] = t*y*y + c;     out->m[2][1] = t*y*z - s*x;
    out->m[0][2] = t*x*z - s*y;   out->m[1][2] = t*y*z + s*x;   out->m[2][2] = t*z*z + c;
}

void k26m3d_mat4_look_at(K26M4 *out, K26V3 eye, K26V3 center, K26V3 up)
{
    /* Camera looks toward center from eye. f = forward (eye -> center),
     * s = side = f x up_world, u = up = s x f. View matrix has rows
     * (s, u, -f) and translation -(s.eye, u.eye, -f.eye). */
    K26V3 f = k26m3d_v3_norm(k26m3d_v3_sub(center, eye));
    K26V3 s = k26m3d_v3_norm(k26m3d_v3_cross(f, up));
    K26V3 u = k26m3d_v3_cross(s, f);

    k26m3d_mat4_identity(out);
    out->m[0][0] =  s.x;   out->m[1][0] =  s.y;   out->m[2][0] =  s.z;
    out->m[0][1] =  u.x;   out->m[1][1] =  u.y;   out->m[2][1] =  u.z;
    out->m[0][2] = -f.x;   out->m[1][2] = -f.y;   out->m[2][2] = -f.z;
    out->m[3][0] = -k26m3d_v3_dot(s, eye);
    out->m[3][1] = -k26m3d_v3_dot(u, eye);
    out->m[3][2] =  k26m3d_v3_dot(f, eye);
}

void k26m3d_mat4_perspective(K26M4 *out,
                             double fov_y_rad, double aspect,
                             double near_plane, double far_plane)
{
    double f = 1.0 / tan(fov_y_rad * 0.5);
    double nf = 1.0 / (near_plane - far_plane);

    k26m3d_mat4_zero(out);
    out->m[0][0] = f / aspect;
    out->m[1][1] = f;
    out->m[2][2] = (far_plane + near_plane) * nf;
    out->m[2][3] = -1.0;
    out->m[3][2] = 2.0 * far_plane * near_plane * nf;
}

void k26m3d_mat4_orthographic(K26M4 *out,
                              double left, double right,
                              double bottom, double top,
                              double near_plane, double far_plane)
{
    double rl = 1.0 / (right - left);
    double tb = 1.0 / (top - bottom);
    double fn = 1.0 / (far_plane - near_plane);

    k26m3d_mat4_zero(out);
    out->m[0][0] =  2.0 * rl;
    out->m[1][1] =  2.0 * tb;
    out->m[2][2] = -2.0 * fn;
    out->m[3][0] = -(right + left)       * rl;
    out->m[3][1] = -(top + bottom)       * tb;
    out->m[3][2] = -(far_plane + near_plane) * fn;
    out->m[3][3] =  1.0;
}
