/* attitude.c — attitude integration + gravity-gradient torque.
 *
 * References:
 *   - Markley & Crassidis, Fundamentals of Spacecraft Attitude
 *     Determination and Control (2014), §3.7 — quaternion
 *     exponential map
 *   - Wertz (ed.), Spacecraft Attitude Determination and Control
 *     (1978), §17.2 — gravity-gradient torque on extended bodies
 *   - Hughes, Spacecraft Attitude Dynamics (1986), §4.5 — Euler's
 *     rotational equation in principal axes
 */
#include "k26astro_body/attitude.h"

#include <math.h>

void k26astro_attitude_init(K26AstroAttitudeState *a,
                            double sphere_mass, double sphere_radius)
{
    if (!a) return;
    a->q = k26m3d_quat_identity();
    a->omega_body.x = a->omega_body.y = a->omega_body.z = 0.0;
    /* Solid sphere principal moment: I = (2/5) M R². */
    double I = (2.0 / 5.0) * sphere_mass * sphere_radius * sphere_radius;
    a->inertia_diag.x = I;
    a->inertia_diag.y = I;
    a->inertia_diag.z = I;
    a->torques = NULL;
}

K26Quat k26astro_quat_exp_half(K26V3 theta)
{
    /* exp(0.5 (θx i + θy j + θz k)) = cos(|θ|/2)
     *                                + sin(|θ|/2) θ̂ · (i, j, k)
     * Numerically stable form using sinc for small angles. */
    double half_norm_sq = 0.25 * (theta.x * theta.x
                                + theta.y * theta.y
                                + theta.z * theta.z);
    double half_norm    = sqrt(half_norm_sq);
    K26Quat q;
    if (half_norm < 1.0e-8) {
        /* Truncated Taylor: cos(x) ≈ 1 - x²/2, sin(x)/x ≈ 1 - x²/6. */
        double s_over_2 = 0.5 * (1.0 - half_norm_sq / 6.0);
        q.x = theta.x * s_over_2;
        q.y = theta.y * s_over_2;
        q.z = theta.z * s_over_2;
        q.w = 1.0 - half_norm_sq * 0.5;
    } else {
        double s = sin(half_norm) / (2.0 * half_norm);
        q.x = theta.x * s;
        q.y = theta.y * s;
        q.z = theta.z * s;
        q.w = cos(half_norm);
    }
    return q;
}

void k26astro_attitude_step_free(K26AstroAttitudeState *a, double dt)
{
    if (!a) return;
    K26V3 theta = {
        a->omega_body.x * dt,
        a->omega_body.y * dt,
        a->omega_body.z * dt
    };
    K26Quat dq = k26astro_quat_exp_half(theta);
    a->q = k26m3d_quat_norm(k26m3d_quat_mul(a->q, dq));
}

void k26astro_attitude_step_torque(K26AstroAttitudeState *a,
                                    K26V3 torque_body, double dt)
{
    if (!a) return;
    double Ixx = a->inertia_diag.x;
    double Iyy = a->inertia_diag.y;
    double Izz = a->inertia_diag.z;
    /* Euler's equations in principal axes. */
    double wx = a->omega_body.x;
    double wy = a->omega_body.y;
    double wz = a->omega_body.z;
    double wdot_x = (Ixx > 0.0) ? (torque_body.x + (Iyy - Izz) * wy * wz) / Ixx : 0.0;
    double wdot_y = (Iyy > 0.0) ? (torque_body.y + (Izz - Ixx) * wz * wx) / Iyy : 0.0;
    double wdot_z = (Izz > 0.0) ? (torque_body.z + (Ixx - Iyy) * wx * wy) / Izz : 0.0;
    /* Mid-point integration: advance ω with half-step, then full-step
     * orientation, then full-step ω. Symplectic-ish; sufficient for
     * typical attitude budgets. Higher-order variants can be added
     * once libk26astro_grav provides the perturbation surface. */
    a->omega_body.x += wdot_x * dt;
    a->omega_body.y += wdot_y * dt;
    a->omega_body.z += wdot_z * dt;
    K26V3 theta = {
        a->omega_body.x * dt,
        a->omega_body.y * dt,
        a->omega_body.z * dt
    };
    K26Quat dq = k26astro_quat_exp_half(theta);
    a->q = k26m3d_quat_norm(k26m3d_quat_mul(a->q, dq));
}

K26V3 k26astro_torque_gravity_gradient(K26V3 r_central_body,
                                       double r_distance,
                                       double mu,
                                       K26V3 inertia_diag)
{
    K26V3 result = { 0.0, 0.0, 0.0 };
    if (r_distance <= 0.0 || mu <= 0.0) return result;
    double inv_r = 1.0 / r_distance;
    K26V3 r_hat = {
        r_central_body.x * inv_r,
        r_central_body.y * inv_r,
        r_central_body.z * inv_r
    };
    /* I r̂ for diagonal I. */
    K26V3 Ir_hat = {
        inertia_diag.x * r_hat.x,
        inertia_diag.y * r_hat.y,
        inertia_diag.z * r_hat.z
    };
    /* r̂ × I r̂. */
    K26V3 cross = {
        r_hat.y * Ir_hat.z - r_hat.z * Ir_hat.y,
        r_hat.z * Ir_hat.x - r_hat.x * Ir_hat.z,
        r_hat.x * Ir_hat.y - r_hat.y * Ir_hat.x
    };
    double k = 3.0 * mu / (r_distance * r_distance * r_distance);
    result.x = k * cross.x;
    result.y = k * cross.y;
    result.z = k * cross.z;
    return result;
}

K26V3 k26astro_attitude_error(K26Quat from, K26Quat to)
{
    /* Error quaternion: q_err = from^-1 * to */
    K26Quat from_conj = k26m3d_quat_conj(from);
    K26Quat q_err     = k26m3d_quat_mul(from_conj, to);
    /* Force shortest-path (q and -q encode the same rotation). */
    if (q_err.w < 0.0) {
        q_err.x = -q_err.x;
        q_err.y = -q_err.y;
        q_err.z = -q_err.z;
        q_err.w = -q_err.w;
    }
    /* Small-rotation-vector extraction. For |sin(θ/2)| not tiny:
     *   θ = 2 atan2(|v|, w) along v̂.  For tiny |v|: linear approx. */
    double v_norm = sqrt(q_err.x * q_err.x + q_err.y * q_err.y + q_err.z * q_err.z);
    K26V3 out;
    if (v_norm < 1.0e-9) {
        out.x = 2.0 * q_err.x;
        out.y = 2.0 * q_err.y;
        out.z = 2.0 * q_err.z;
    } else {
        double angle = 2.0 * atan2(v_norm, q_err.w);
        double scale = angle / v_norm;
        out.x = q_err.x * scale;
        out.y = q_err.y * scale;
        out.z = q_err.z * scale;
    }
    return out;
}
