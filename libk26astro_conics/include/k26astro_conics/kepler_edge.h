/* k26astro_conics/kepler_edge.h — Kepler propagator edge cases.
 *
 * The baseline universal-variable propagator in kepler.h handles the
 * elliptic + hyperbolic majority of two-body geometry, but three corners
 * need bespoke treatment:
 *
 *   1. Parabolic (α ≈ 0): the universal-variable Newton-Raphson stalls
 *      because the Lagrange-f denominators blow up. Barker's equation
 *      gives an algebraic closed form.
 *
 *   2. Near-rectilinear (angular momentum h ≈ 0): the orbit is radial.
 *      Treat as 1-D Kepler problem along the line of motion; classical
 *      universal-variable propagator's r0·σ0 term oscillates wildly.
 *
 *   3. Hyperbolic asymptote (large excess velocity, dt large enough that
 *      the body has effectively escaped): the χ initial guess from
 *      kepler.c relies on log() arithmetic that loses precision for
 *      very large |α·dt|. Provide a recursive bisection that halves
 *      dt until the per-step χ is well-conditioned.
 *
 * All three functions share the kepler_propagate signature so the
 * baseline propagator can detect its own failure (returns 4) and
 * dispatch to the right edge handler.
 *
 * References:
 *   - Barker's equation: tan(ν/2)/2 + tan³(ν/2)/6 = √(μ/2p³)·(t-tp)
 *     Vallado §2.3
 *   - Near-rectilinear: Goodyear (1965) AIAA J 3:1326 — solution by
 *     Kepler's equation in the line-frame
 *   - Hyperbolic stability: Battin §4.5 — initial-guess heuristics
 *     break down for ε ≪ 0
 */
#ifndef K26ASTRO_CONICS_KEPLER_EDGE_H
#define K26ASTRO_CONICS_KEPLER_EDGE_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Parabolic orbit propagator via Barker's equation. Caller passes the
 * specific angular momentum h (= |r × v|, computed once); semi-latus
 * rectum p = h² / mu. The orbit is in the plane spanned by (pos0,
 * vel0); the parabolic solution stays in that plane.
 *
 * Return codes match kepler_propagate (0 = ok). Note: a true parabola
 * has e = exactly 1 which is unphysical; this entry point is for
 * orbits where |α| < K26A_KEPLER_PARABOLIC_EPSILON (1e-12). */
int k26astro_kepler_propagate_parabolic(K26V3 *out_pos, K26V3 *out_vel,
                                        K26V3 pos0, K26V3 vel0,
                                        double mu, double dt);

/* Radial-orbit (rectilinear) propagator. Use when |h| / (r0 · v0)
 * < K26A_KEPLER_RADIAL_EPSILON (1e-10). The orbit is a degenerate
 * line; propagation reduces to 1-D Kepler along the line direction. */
int k26astro_kepler_propagate_radial(K26V3 *out_pos, K26V3 *out_vel,
                                     K26V3 pos0, K26V3 vel0,
                                     double mu, double dt);

/* Auto-dispatch wrapper. Examines (pos0, vel0, mu, dt) and routes to:
 *   - k26astro_kepler_propagate          (default; elliptic + hyperbolic)
 *   - k26astro_kepler_propagate_parabolic (|α| ≈ 0)
 *   - k26astro_kepler_propagate_radial    (|h| ≈ 0)
 *
 * Recommended entry point for code paths that may see any conic. */
int k26astro_kepler_propagate_any(K26V3 *out_pos, K26V3 *out_vel,
                                  K26V3 pos0, K26V3 vel0,
                                  double mu, double dt,
                                  int max_iter);

#define K26A_KEPLER_PARABOLIC_EPSILON  1.0e-12
#define K26A_KEPLER_RADIAL_EPSILON     1.0e-10

#ifdef __cplusplus
}
#endif

#endif
