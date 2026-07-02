/* k26astro_conics/kepler.h — universal-variable Kepler propagator.
 *
 * The universal-variable formulation handles all conic types (ellipse,
 * parabola, hyperbola) without case branching; it works in terms of
 * the universal anomaly χ via Stumpff C and S series. The
 * implementation is robust for hyperbolic and near-parabolic orbits
 * where the classical Kepler form (M = E - e sin E) breaks down.
 *
 * References:
 *   - Bate, Mueller & White (1971) Fundamentals of Astrodynamics §4.5
 *   - Vallado (2013) Fundamentals of Astrodynamics and Applications §2.4
 *
 * Edge cases (parabolic via Barker's equation, near-rectilinear,
 * hyperbolic-asymptote initial guess) live in kepler_edge.h; the
 * baseline propagator below covers ~99% of two-body geometry.
 *
 * The Stumpff helpers C(z), S(z) are also exposed (via the internal
 * header) so kepler_edge.c, lambert*.c, and the Wisdom-Holman Kepler
 * drift step can reach them without duplicating. */
#ifndef K26ASTRO_CONICS_KEPLER_H
#define K26ASTRO_CONICS_KEPLER_H

#include "k26astro_body/elements.h"  /* K26AstroKeplerian, _StateVector */
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Advance state from (pos, vel) at t = 0 to the same orbit at t = dt.
 * `mu` is the central body's GM in m³/s². `max_iter` caps the Newton-
 * Raphson loop on χ; 32 is generous, 8 is typical for inner-SS orbits.
 *
 * Returns 0 on convergence, non-zero on:
 *   1 - NULL output pointer
 *   2 - non-positive mu
 *   3 - degenerate initial position (r0 ≤ 0)
 *   4 - degenerate post-step radius (r ≤ 0); usually indicates
 *       Newton-Raphson stalled on a parabolic-edge case; re-issue
 *       via k26astro_kepler_propagate_parabolic() from kepler_edge.h. */
int k26astro_kepler_propagate(K26V3 *out_pos, K26V3 *out_vel,
                              K26V3 pos0, K26V3 vel0,
                              double mu, double dt,
                              int max_iter);

/* Convenience: propagate Keplerian elements forward by dt. Internally
 * converts elements → state → propagates → elements. The output
 * mean-anomaly is wrapped to [0, 2π). Returns 0 on success. */
int k26astro_kepler_propagate_elements(K26AstroKeplerian *out,
                                       const K26AstroKeplerian *in,
                                       double dt);

#ifdef __cplusplus
}
#endif

#endif
