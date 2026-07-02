/* k26astro_body/attitude.h — body attitude state + integration.
 *
 * Attitude is a unit quaternion q (body → inertial) plus an angular
 * velocity ω in the body frame. The integrator advances both each
 * substep:
 *
 *   q_{t+dt} = q_t * exp(0.5 * ω_body * dt)        (quaternion exponential)
 *   ω_{t+dt} = ω_t + I⁻¹ (τ - ω × I ω) dt           (Euler equation)
 *
 * The body's moment-of-inertia tensor I is needed for the second
 * equation; I is represented as a diagonal tensor (principal-axis
 * frame), expressed as a K26V3 of the three principal moments.
 * The diagonal form is exact for symmetric bodies and a reasonable
 * approximation for the rest in the absence of high-fidelity inertia
 * measurements.
 *
 * Torques accumulate per substep. Built-in:
 *   - Gravity-gradient torque (libk26astro_body)
 *
 * Caller-supplied torques (SRP, magnetic, thrust, aero moment, etc.)
 * register through k26astro_attitude_register_torque. The registry
 * is owned by the K26AstroAttitudeState; k26astro_attitude_destroy
 * frees it. */
#ifndef K26ASTRO_BODY_ATTITUDE_H
#define K26ASTRO_BODY_ATTITUDE_H

#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Attitude state shape ------------------------------------- *
 *
 * Mirror of the (K26Quat attitude, K26V3 omega) pair stored on
 * K26AstroBody, plus the body's principal-axis inertia tensor.
 *
 * The bundle is exposed as a separate struct so the integrator can
 * operate on it independent of the full K26AstroBody — useful for
 * stand-alone attitude propagation tests + detached attitude state
 * for spacecraft dynamics carried by libk26astro_vehicle. */

/* Opaque torque-source registry. NULL until the first
 * k26astro_attitude_register_torque call; freed by
 * k26astro_attitude_destroy or k26astro_attitude_clear_torques. */
typedef struct K26AstroTorqueList K26AstroTorqueList;

typedef struct {
    K26Quat              q;            /* body → inertial */
    K26V3                omega_body;   /* angular velocity, body frame, rad/s */
    K26V3                inertia_diag; /* principal moments (Ixx, Iyy, Izz), kg·m² */
    K26AstroTorqueList  *torques;      /* nullable; lazy-allocated */
} K26AstroAttitudeState;

/* Identity initialisation: q = identity, omega = zero, inertia
 * uniform = mass * radius² * 2/5 (sphere). */
void k26astro_attitude_init(K26AstroAttitudeState *a,
                            double sphere_mass, double sphere_radius);

/* ---- Quaternion utilities ------------------------------------- */
/* Quaternion exponential — for a small rotation vector θ (rad),
 * returns exp(0.5 θ_quat) where θ_quat = (θ.x, θ.y, θ.z, 0).
 * Numerically stable for both |θ| → 0 and large rotations. */
K26Quat k26astro_quat_exp_half(K26V3 theta);

/* ---- Integration ---------------------------------------------- *
 *
 * Free-body integration step (no external torques): rotate the
 * orientation by ω * dt, leave ω unchanged.
 * Applies the small-angle exponential map then renormalises. */
void k26astro_attitude_step_free(K26AstroAttitudeState *a, double dt);

/* Integration step with an external body-frame torque. Updates both
 * orientation and angular velocity per the Euler rotational
 * equations:
 *   Ixx ω̇x = τx + (Iyy - Izz) ωy ωz
 *   Iyy ω̇y = τy + (Izz - Ixx) ωz ωx
 *   Izz ω̇z = τz + (Ixx - Iyy) ωx ωy
 * Then orientation advances via the quaternion exponential map. */
void k26astro_attitude_step_torque(K26AstroAttitudeState *a,
                                    K26V3 torque_body, double dt);

/* ---- Gravity-gradient torque -------------------------------- *
 *
 * For a body orbiting in a gravitational potential V, the gravity
 * gradient on an extended (non-point-mass) body produces a torque
 * proportional to the inertia tensor and the radial direction:
 *
 *   τ_gg = (3μ / r³) (r̂ × I r̂)
 *
 * where r̂ is the body-frame unit vector pointing from the body
 * centre to the central body, μ is the central body's GM, and I is
 * the body's inertia tensor (in body frame). Returns the torque in
 * body frame. */
K26V3 k26astro_torque_gravity_gradient(K26V3 r_central_body,
                                       double r_distance,
                                       double mu,
                                       K26V3 inertia_diag);

/* ---- Torque-source registry --------------------------------- *
 *
 * Each registered function adds its body-frame contribution to
 * `torque_body_out` (additive — do not overwrite). The summed
 * registry torque is the input to a single Euler step at substep
 * close. Registration order is preserved; the registry is dispatched
 * sequentially in that order. */
typedef void (*K26AstroTorqueFn)(const K26AstroAttitudeState *state,
                                 double t,
                                 K26V3 *torque_body_out,
                                 void  *ctx);

/* Register a torque source. Returns 0 on success, 1 on null input,
 * 2 on allocation failure. Lifetime of `ctx` is the caller's
 * responsibility. */
int  k26astro_attitude_register_torque(K26AstroAttitudeState *state,
                                       K26AstroTorqueFn fn,
                                       void *ctx);

/* Drop all registered torque sources and free the registry buffer.
 * No-op if no sources were ever registered. */
void k26astro_attitude_clear_torques  (K26AstroAttitudeState *state);

/* Free any allocated registry storage. Safe to call on stack-
 * allocated POD states that never registered a torque (no-op when
 * `state->torques == NULL`). New code that registers torques must
 * pair the lifecycle with this call. */
void k26astro_attitude_destroy        (K26AstroAttitudeState *state);

/* Sum the registered torques (calling each with `t` and the
 * accumulator) and apply one Euler step of duration `dt`. Equivalent
 * to k26astro_attitude_step_torque with the registry sum as input. */
void k26astro_attitude_step_torque_registry(K26AstroAttitudeState *state,
                                            double t, double dt);

/* ---- Full 3×3-inertia attitude state (spacecraft) ----------- *
 *
 * Sibling of K26AstroAttitudeState carrying a full 3×3 inertia
 * tensor (body frame) and its precomputed inverse. For vehicles
 * whose mass distribution is not aligned with principal axes — the
 * inertia tensor has off-diagonal coupling — and for vehicles whose
 * inertia changes mid-simulation due to stage events or propellant
 * depletion.
 *
 * Celestial bodies keep the diagonal-inertia K26AstroAttitudeState
 * path; this struct is the spacecraft path. The Ext-path torque
 * registry uses its own K26AstroTorqueFnExt signature (below) so
 * torque sources see the full Ext state and can read the inertia
 * tensor directly. */
typedef struct {
    K26Quat              q;               /* body → inertial */
    K26V3                omega_body;      /* angular velocity, body frame, rad/s */
    K26M3                inertia;         /* full 3×3, body frame, kg·m² */
    K26M3                inertia_inverse; /* precomputed; updated on inertia change */
    K26AstroTorqueList  *torques;         /* nullable; lazy-allocated */
} K26AstroAttitudeStateExt;

/* Initialise q = identity, ω = zero, inertia = `inertia` (caller-
 * supplied; assumed symmetric positive-definite), inverse computed.
 * If `inertia` is singular (|det| < 1e-30), inertia_inverse is set
 * to the zero matrix and the state is unusable for torque steps. */
void k26astro_attitude_init_ext   (K26AstroAttitudeStateExt *a,
                                    K26M3 inertia);

/* Free the torque-source registry; no-op if none was registered.
 * Required for any state that called k26astro_attitude_register_
 * torque_ext. */
void k26astro_attitude_destroy_ext(K26AstroAttitudeStateExt *a);

/* Free-body integration step (no external torques): rotate the
 * orientation by ω · dt; ω evolves under the symmetric-inertia
 * cross-coupling I⁻¹(ω × Iω). For non-principal-axis ω on an
 * asymmetric inertia the body-frame ω traces a polhode on the
 * inertia ellipsoid (Euler nutation). */
void k26astro_attitude_step_free_ext  (K26AstroAttitudeStateExt *a,
                                        double dt);

/* Integration step with an external body-frame torque. Implements
 * the full Euler rotational equation in body frame:
 *
 *   ω̇ = I⁻¹ (τ_body − ω × (I ω))
 *
 * then advances orientation via the small-angle quaternion
 * exponential map. Uses the precomputed inertia_inverse. */
void k26astro_attitude_step_torque_ext(K26AstroAttitudeStateExt *a,
                                        K26V3 torque_body, double dt);

/* Replace the inertia tensor (mid-simulation event — stage drop,
 * propellant depletion, deployable extension) and recompute the
 * inverse. q + ω are preserved. The discontinuity in inertia is
 * the caller's responsibility to model physically (angular momentum
 * is in general not conserved across an inertia change unless the
 * caller explicitly enforces it by rescaling ω). */
void k26astro_attitude_update_inertia (K26AstroAttitudeStateExt *a,
                                        K26M3 new_inertia);

/* Ext-state torque-source signature. Mirror of K26AstroTorqueFn for
 * the full-inertia path; each registered fn adds its body-frame
 * contribution to `torque_body_out`. */
typedef void (*K26AstroTorqueFnExt)(const K26AstroAttitudeStateExt *state,
                                    double t,
                                    K26V3 *torque_body_out,
                                    void  *ctx);

/* Register / step a per-state torque registry. Independent of the
 * K26AstroAttitudeState registry. Returns 0/1/2 per the diagonal
 * path's register call. */
int  k26astro_attitude_register_torque_ext(K26AstroAttitudeStateExt *a,
                                            K26AstroTorqueFnExt fn,
                                            void *ctx);

void k26astro_attitude_clear_torques_ext  (K26AstroAttitudeStateExt *a);

void k26astro_attitude_step_torque_registry_ext(
    K26AstroAttitudeStateExt *a, double t, double dt);

/* ---- Attitude error ----------------------------------------- */
/* Returns the small-rotation-vector error from `from` to `to`:
 * a vector whose magnitude is the rotation angle and whose
 * direction is the rotation axis. Useful for control-loop residuals
 * and integrator divergence checks. */
K26V3 k26astro_attitude_error(K26Quat from, K26Quat to);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_BODY_ATTITUDE_H */
