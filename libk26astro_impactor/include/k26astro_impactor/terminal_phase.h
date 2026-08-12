/* terminal_phase.h — closing-velocity + impact-geometry helpers.
 *
 * At terminal phase the impactor / swarm unit is treated as a ballistic
 * body; the world propagator advances its state. This lib
 * provides the small kinematic helpers consumers use to derive
 * impact velocity, time-to-close, and cos(impact angle) from the
 * relative state vector at intercept. */
#ifndef K26ASTRO_IMPACTOR_TERMINAL_PHASE_H
#define K26ASTRO_IMPACTOR_TERMINAL_PHASE_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Closing speed (m/s) = magnitude of the relative velocity. */
double k26astro_impactor_closing_speed(K26V3 v_projectile, K26V3 v_target);

/* Time-to-closest-approach (s) for two ballistic bodies. Returns
 * the time at which |r_proj(t) + v_proj·t - (r_tgt + v_tgt·t)|
 * is minimal. Negative return → the projectile is already past
 * its closest approach (overshoot). */
double k26astro_impactor_time_to_closest_approach(
    K26V3 r_projectile, K26V3 v_projectile,
    K26V3 r_target,     K26V3 v_target);

/* Closest-approach distance (m) over a straight-line ballistic
 * encounter. */
double k26astro_impactor_closest_approach_distance(
    K26V3 r_projectile, K26V3 v_projectile,
    K26V3 r_target,     K26V3 v_target);

/* Cosine of the impact angle between the closing velocity vector
 * and the target's surface normal. cos θ in (0, 1]. */
double k26astro_impactor_impact_cos_angle(K26V3 closing_velocity_unit,
                                      K26V3 surface_normal_unit);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_IMPACTOR_TERMINAL_PHASE_H */
