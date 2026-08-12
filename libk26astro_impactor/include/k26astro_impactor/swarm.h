/* swarm.h — multi-unit swarm dispatch geometry.
 *
 * A swarm-pattern impactor releases N small projectiles in a
 * cone of half-angle θ at terminal phase. Each unit follows
 * an independent ballistic trajectory under the same N-body
 * propagator; the lib provides the geometric dispatch helpers
 * the driver uses to spawn the world bodies.
 *
 * Dispatch sampling: the dispersion cone is uniformly sampled in
 * solid angle (each unit's direction is independent and
 * isotropic within the cone). Per-unit mass is the swarm's
 * total released mass divided by N.
 *
 * Determinism: the per-unit direction draws are seeded
 * from the caller-supplied K26CRng (world RNG in REFERENCED
 * mode). */
#ifndef K26ASTRO_IMPACTOR_SWARM_H
#define K26ASTRO_IMPACTOR_SWARM_H

#include "k26m3d.h"
#include "k26compute.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Sample a single swarm-unit direction inside a cone of
 * half-angle θ aligned with `cone_axis_unit`. The sampling is
 * uniform on the cap's solid angle.
 *
 * Returns a unit vector in the same coordinate frame as
 * cone_axis_unit. NULL rng returns the cone axis itself
 * (deterministic on-axis launch). */
K26V3 k26astro_swarm_sample_direction(K26V3 cone_axis_unit,
                                           double half_angle_rad,
                                           K26CRng *rng);

/* Per-unit mass (kg) for a swarm of `count` units totalling
 * `total_mass_kg`. */
double k26astro_swarm_per_unit_mass(int count,
                                         double total_mass_kg);

/* Effective coverage area of the swarm at range R from the
 * launcher (m^2). The cone's projected footprint at range R is
 * π (R sin θ)^2; with N units, the per-unit projected area is
 * footprint / N. */
double k26astro_swarm_footprint_m2(double range_m,
                                              double half_angle_rad);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_IMPACTOR_SWARM_H */
