/* k26astro_conics/soi.h — sphere of influence: radius + hierarchy
 * + crossing detection.
 *
 * Sphere-of-influence (SOI) is the locus around a secondary body where
 * its gravitational pull dominates over its primary's. Inside the SOI,
 * patched-conic propagation treats the secondary as the central body;
 * outside, the primary. The on-rails / patched-conic propagation
 * mode is a stack of SOI containments: each body has a parent (the
 * body whose SOI it's currently inside) and may itself contain child
 * bodies in its own SOI.
 *
 * Two radius formulae:
 *
 *   Hill radius:    r_H = a · (m₂ / (3·m₁))^(1/3)
 *     The classical three-body Jacobi-integral-derived radius. Most
 *     common in dynamical astronomy literature. Slightly larger than
 *     the Laplace radius for the same body pair.
 *
 *   Laplace radius: r_L = a · (m₂ / m₁)^(2/5)
 *     Derived from the perturbation-equals-central-force criterion.
 *     The traditional spacecraft-trajectory choice; matches Lagrange
 *     and Laplace's original formulation.
 *
 * Both are first-order approximations to a fundamentally fuzzy
 * concept. K26 defaults to Laplace for backward compatibility with
 * mission-planning tooling; Hill is provided for cases where the
 * literature being matched uses it (e.g. solar-system dynamics
 * papers).
 *
 * Hierarchy: each K26AstroBody carries a `parent_body_idx` field
 * (-1 = no parent, body is on its own around the universe origin).
 * The hierarchy walker resolves a body's current parent (by examining
 * its position relative to each candidate's SOI) and returns the
 * deepest containing parent.
 *
 * Crossing detection: given a body's state at t and at t+dt, plus
 * the parent's SOI radius, return whether the body's distance from
 * parent crossed the SOI boundary, and if so a refined crossing
 * time. Linear interpolation in r(t); sufficient for the patched-
 * conic handoff threshold (typical SOI radius is a few % of orbital
 * radius; sub-step refinement is the MERCURIUS orchestrator's
 * concern). */
#ifndef K26ASTRO_CONICS_SOI_H
#define K26ASTRO_CONICS_SOI_H

#include "k26astro_body/body.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- SOI radii --------------------------------------------------- */

/* Hill radius: r_H = a · ∛(m₂ / (3·m₁)) = a · ∛(GM₂ / (3·GM₁)). */
double k26astro_soi_radius_hill(double parent_gm, double child_gm,
                                double semi_major_axis);

/* Laplace radius: r_L = a · (m₂ / m₁)^(2/5) = a · (GM₂/GM₁)^(2/5). */
double k26astro_soi_radius_laplace(double parent_gm, double child_gm,
                                   double semi_major_axis);

/* ---- Hierarchy walk --------------------------------------------- */

/* Resolve `child_idx`'s parent at time `t`. Walks the parent_body_idx
 * chain and confirms each candidate's SOI containment via Laplace
 * radius (using the parent-of-parent for the radius's "a" argument).
 *
 * Returns the parent body's index, or -1 if `child_idx` has no
 * parent. The `t` argument is reserved for time-varying SOI hierarchies
 * (moons whose parent depends on the phase of motion); it is
 * currently accepted but unused; the static parent_body_idx field
 * is taken as authoritative. */
int k26astro_soi_parent(const K26AstroBody *bodies, int n_bodies,
                        int child_idx, double t);

/* ---- Crossing detection ----------------------------------------- *
 *
 * Given a body's pre-step and post-step state, plus the parent's
 * GM and SOI radius, return whether the body crossed the SOI
 * boundary during dt. If it did, *out_t_cross is a fraction in
 * [0, 1] representing the linear-interpolated crossing time
 * (multiply by dt to get seconds since pre-step).
 *
 * Returns:
 *    1 - crossing detected (out_t_cross populated)
 *    0 - no crossing (body stayed inside or stayed outside)
 *   -1 - invalid input */
int k26astro_soi_crossing_detect(K26V3 pos_before_rel_parent,
                                  K26V3 pos_after_rel_parent,
                                  double soi_radius,
                                  double *out_t_cross);

#ifdef __cplusplus
}
#endif

#endif
