/* k26astro_conics/lambert_multi.h — multi-revolution Lambert solver.
 *
 * Izzo's algorithm (Izzo 2014/2015 CeMec 121:1, "Revisiting Lambert's
 * problem") solves Lambert's problem in a numerically robust form
 * across single + multi-revolution cases. The variable x parameterises
 * the orbital geometry: x ∈ (-1, 1) for elliptic, x > 1 for hyperbolic;
 * λ encodes the chord-vs-semiperimeter ratio with sign for
 * short-way/long-way.
 *
 * For n_rev ≥ 1, two distinct elliptic transfers exist per revolution
 * count: low Δv (low TOF on the right branch) vs high Δv (high TOF
 * on the left branch). The `branch` parameter selects.
 *
 * Caveats:
 *   - For n_rev = 0, the result matches the single-rev solver from
 *     lambert.h (within iteration tolerance) for all non-degenerate
 *     geometries.
 *   - n_rev ≥ 1 requires elliptic geometry — hyperbolic orbits don't
 *     close, so they have no multi-rev branches. The solver returns
 *     K26ASTRO_LAMBERT_NO_SOLUTION if the supplied tof is too short
 *     for the requested n_rev.
 *   - "Degenerate" return code is reserved for genuinely
 *     near-rectilinear geometry where λ → 1; falls back to caller's
 *     single-rev path. */
#ifndef K26ASTRO_CONICS_LAMBERT_MULTI_H
#define K26ASTRO_CONICS_LAMBERT_MULTI_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Branch selector (n_rev ≥ 1). */
#define K26A_LAMBERT_LOW_DV    0   /* right branch (low Δv) */
#define K26A_LAMBERT_HIGH_DV   1   /* left branch  (high Δv) */

/* Direction (shared with lambert.h). */
#define K26A_LAMBERT_PROGRADE   0
#define K26A_LAMBERT_RETROGRADE 1

/* Error codes specific to multi-rev. */
#define K26A_LAMBERT_OK              0
#define K26A_LAMBERT_NULL_OUT        1
#define K26A_LAMBERT_BAD_INPUT       2
#define K26A_LAMBERT_DEGENERATE      3
#define K26A_LAMBERT_NO_SOLUTION     4
#define K26A_LAMBERT_NO_CONVERGE     5

/* Solve Lambert's problem allowing n_rev ≥ 0 revolutions.
 *
 * Inputs:
 *   r1, r2     — position vectors (m), same inertial frame
 *   mu         — GM of central body (m³/s²)
 *   tof        — time of flight (s), > 0
 *   n_rev      — number of full revolutions, ≥ 0
 *   direction  — K26A_LAMBERT_PROGRADE or _RETROGRADE
 *   branch     — K26A_LAMBERT_LOW_DV or _HIGH_DV (only used if n_rev ≥ 1)
 *
 * Outputs:
 *   *out_v1, *out_v2 - transfer-orbit velocities at r1, r2 (m/s)
 *
 * Returns K26A_LAMBERT_OK (0) on convergence; non-zero per the
 * K26A_LAMBERT_* codes above. */
int k26astro_lambert_multi_rev(K26V3 *out_v1, K26V3 *out_v2,
                               K26V3 r1, K26V3 r2,
                               double mu, double tof,
                               int n_rev,
                               int direction,
                               int branch);

/* Diagnostic test hook exposing the internal T(x, λ, n_rev)
 * evaluator. Validates against the reference CSV (truth from the
 * single-branch Lagrange form at
 * libk26astro_conics/tests/data/lambert_tof_reference.csv). Any
 * divergence from the truth at points where the Lagrange branch
 * itself is accurate (|x| < 0.7) localises the Battin/Lagrange
 * boundary scaling discrepancy.
 *
 * Returns T(x, λ, n_rev). x must be in (-1, 1), |λ| < 1. Not part
 * of the public Lambert API surface; exposed only for test
 * validation. The function is implemented inline at the bottom of
 * lambert_multi.c. */
double k26astro_lambert_tof_for_test(double x, double lambda, int n_rev);

#ifdef __cplusplus
}
#endif

#endif
