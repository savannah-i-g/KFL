/* k26astro_conics/lambert.h — Lambert's problem.
 *
 * Lambert's problem: given two position vectors r₁, r₂ and a time-of-
 * flight tof, find the connecting transfer orbit. The classical
 * single-revolution form has a unique solution per direction (short-way
 * vs long-way). The multi-revolution form (kepler_multi.h) allows
 * n ≥ 1 full orbits between the two endpoints.
 *
 * Implementation: universal-variable form of Vallado §7.6 / Algorithm 58,
 * which is itself a clean-room re-derivation of Battin's original method
 * (Battin 1987, "An Introduction to the Mathematics and Methods of
 * Astrodynamics" §7.3). The universal variable z = α χ² unifies elliptic
 * + hyperbolic transfers without case branching; convergence is
 * Newton-Raphson on z.
 *
 * Edge cases handled by the multi-rev solver (lambert_multi.h):
 *   - Multi-revolution branches (low Δv vs high Δv per n_rev)
 *   - Near-rectilinear geometry (parabolic limit on a single rev)
 *   - 180° transfers (z indeterminate; falls back to a direction prompt)
 *
 * Reference: Vallado (2013) §7.6 + Algorithm 58. */
#ifndef K26ASTRO_CONICS_LAMBERT_H
#define K26ASTRO_CONICS_LAMBERT_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Direction modes. */
#define K26A_LAMBERT_SHORT_WAY  0   /* prograde / θ < π */
#define K26A_LAMBERT_LONG_WAY   1   /* retrograde / θ > π */

/* Solve Lambert's problem (single revolution).
 *
 * Inputs:
 *   r1, r2       — position vectors (m), same inertial frame
 *   mu           — GM of central body (m³/s²)
 *   tof          — time of flight (s), > 0
 *   direction    — K26A_LAMBERT_SHORT_WAY or _LONG_WAY
 *
 * Outputs:
 *   *out_v1      — velocity at r1 on the transfer orbit (m/s)
 *   *out_v2      — velocity at r2 on the transfer orbit (m/s)
 *
 * Returns 0 on convergence; non-zero on:
 *   1 — NULL output
 *   2 — non-positive mu or tof
 *   3 — degenerate geometry (|r1| ≤ 0 or |r2| ≤ 0)
 *   4 — Newton-Raphson failed to converge in max_iter=64
 *   5 — 180° transfer with no direction hint (caller must pick a plane) */
int k26astro_lambert(K26V3 *out_v1, K26V3 *out_v2,
                     K26V3 r1, K26V3 r2,
                     double mu, double tof,
                     int direction);

#ifdef __cplusplus
}
#endif

#endif
