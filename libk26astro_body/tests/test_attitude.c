/* test_attitude.c — quaternion + ω integration + gravity-gradient
 * torque tests. */
#include "k26astro_body/attitude.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

static double quat_norm_sq_(K26Quat q)
{
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

int main(void)
{
    /* ---- Identity init ------------------------------------- */
    K26AstroAttitudeState a;
    k26astro_attitude_init(&a, 1.0, 1.0);
    assert(near_(quat_norm_sq_(a.q), 1.0, 1e-12));
    assert(a.omega_body.x == 0.0);
    assert(a.omega_body.y == 0.0);
    assert(a.omega_body.z == 0.0);
    /* Solid-sphere inertia: I = 2/5 * 1 * 1² = 0.4. */
    assert(near_(a.inertia_diag.x, 0.4, 1e-12));

    /* ---- Free spin: rotate 90° around z over 1 second ------ */
    a.omega_body.x = 0.0;
    a.omega_body.y = 0.0;
    a.omega_body.z = K26A_HALF_PI;
    k26astro_attitude_step_free(&a, 1.0);
    /* q ≈ (0, 0, sin(45°), cos(45°)). */
    assert(fabs(a.q.z - sin(K26A_PI / 4.0)) < 1e-10);
    assert(fabs(a.q.w - cos(K26A_PI / 4.0)) < 1e-10);
    assert(near_(quat_norm_sq_(a.q), 1.0, 1e-12));

    /* ---- Gravity-gradient torque on an elongated body ------ *
     * Body inertia: long axis along z. Place along +x of a
     * gravitating body at 1 AU. Torque should rotate the body
     * to align long axis radially. */
    K26AstroAttitudeState eb;
    k26astro_attitude_init(&eb, 1.0, 1.0);
    eb.inertia_diag.x = 1.0;
    eb.inertia_diag.y = 1.0;
    eb.inertia_diag.z = 2.0;   /* elongated along z */

    K26V3 r_central = { K26A_AU_M, 0.0, 0.0 };   /* central body at +x */
    double r_dist = K26A_AU_M;
    K26V3 tau = k26astro_torque_gravity_gradient(r_central, r_dist,
                                                  K26A_GM_SUN,
                                                  eb.inertia_diag);
    /* For long axis along z and r along x:
     *   r̂ × I r̂ = x̂ × (Ixx * x̂)
     * Since x̂ × x̂ = 0, the dominant term cancels;
     * but tiny torque due to off-diagonal terms = 0 (diagonal I).
     * Actually r̂ × (I r̂) = (1,0,0) × (1,0,0) = 0. The torque is
     * zero when r̂ is aligned with a principal axis. */
    assert(fabs(tau.x) < 1e-9);
    assert(fabs(tau.y) < 1e-9);
    assert(fabs(tau.z) < 1e-9);

    /* Tilt: r̂ = (cos45°, sin45°, 0). Now r̂ × I r̂ should be
     * nonzero in z. */
    double s45 = sqrt(2.0) / 2.0;
    K26V3 r_tilt = { K26A_AU_M * s45, K26A_AU_M * s45, 0.0 };
    tau = k26astro_torque_gravity_gradient(r_tilt, K26A_AU_M,
                                            K26A_GM_SUN,
                                            eb.inertia_diag);
    /* (Ixx - Iyy) = 0, so torque has only z component zero too.
     * But (Iyy - Izz) ≠ 0; so x-component should be nonzero when
     * r has y and z components. With z=0, x-component is also 0.
     * Re-tilt out of plane. */
    K26V3 r_3d = { K26A_AU_M * 0.5, K26A_AU_M * 0.5, K26A_AU_M * 0.5 };
    double r_3d_mag = sqrt(0.25 + 0.25 + 0.25) * K26A_AU_M;
    tau = k26astro_torque_gravity_gradient(r_3d, r_3d_mag,
                                            K26A_GM_SUN,
                                            eb.inertia_diag);
    /* Torque magnitude > 0 because I is asymmetric and r̂ is generic. */
    double tau_mag = sqrt(tau.x * tau.x + tau.y * tau.y + tau.z * tau.z);
    assert(tau_mag > 0.0);

    /* ---- Attitude error round-trip ------------------------ */
    K26Quat q1 = k26m3d_quat_identity();
    K26Quat q2 = k26m3d_quat_from_axis_angle(
        k26m3d_v3(0, 0, 1), K26A_HALF_PI);
    K26V3 err = k26astro_attitude_error(q1, q2);
    assert(fabs(err.x) < 1e-10);
    assert(fabs(err.y) < 1e-10);
    assert(fabs(err.z - K26A_HALF_PI) < 1e-10);

    printf("test_attitude: OK\n");
    return 0;
}
