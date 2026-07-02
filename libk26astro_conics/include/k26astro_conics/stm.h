/* k26astro_conics/stm.h — Cartesian state-transition matrix for the
 * universal-variable Kepler propagator.
 *
 * Returns the 6 x 6 Jacobian M of the final Cartesian state with
 * respect to the initial Cartesian state under pure two-body
 * dynamics:
 *
 *     state(t)   = [ r(t), v(t) ]
 *     state(0)   = [ r0,   v0   ]
 *     M[i][j]    = d state(t)[i] / d state(0)[j]
 *
 * The (d r / d v0) sub-block (rows 0-2, columns 3-5) is the
 * impulsive-delta-v -> arrival-position partial that consumers of
 * the Newton-iterate trajectory targeting (libk26astro_traj) use
 * for analytic Jacobian columns. The full 6 x 6 is exposed so
 * orbit-determination partial derivatives can reuse the same
 * surface.
 *
 * The implementation drives centered finite differences through the
 * existing universal-variable propagator (kepler.h). Battin (1999)
 * paragraph 10.5 and Vallado (2013) paragraph 10.7.2 derive a
 * closed-form expression that is mathematically equivalent and
 * avoids the FD step-size choice; the closed form is significantly
 * more code and offers no asymptotic speed advantage over the
 * centered FD path because the bulk cost is the Kepler propagation
 * itself, which the FD path amortises across all twelve perturbed
 * evaluations. Centered O(h^2) FD at h = 1e-6 times the input
 * scale; the regression suite empirically holds the dr/dv0
 * sub-block agreement to ~1e-6 relative across LEO / GTO /
 * hyperbolic test fixtures. */
#ifndef K26ASTRO_CONICS_STM_H
#define K26ASTRO_CONICS_STM_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Compute the 6 x 6 STM from (r0, v0) advanced by dt under pure
 * two-body dynamics with central GM mu. The output matrix is
 * row-major: M_out[i][j] = d state_final[i] / d state_initial[j],
 * where state = (r.x, r.y, r.z, v.x, v.y, v.z).
 *
 * Returns 0 on success, non-zero on:
 *   1 — NULL output, non-positive mu, or |r0| == 0
 *   2 — internal Kepler propagation failed on a perturbed column
 *       (degenerate edge case; the caller should treat this as a
 *       request to fall back to a robust integrator). */
int k26astro_kepler_stm(K26V3 r0, K26V3 v0, double dt, double mu,
                        double M_out[6][6]);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CONICS_STM_H */
