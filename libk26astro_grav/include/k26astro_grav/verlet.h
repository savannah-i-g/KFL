/* k26astro_grav/verlet.h — DKD (Drift-Kick-Drift) Velocity Verlet.
 *
 * The simplest symplectic integrator: second-order accurate, energy
 * bounded over arbitrary integration spans, and dirt cheap. Step:
 *
 *     x → x + (dt/2)·v
 *     v → v + dt·a(x_new)
 *     x → x + (dt/2)·v_new
 *
 * The right choice when the system has no dominant central mass
 * (equal-mass clusters, dust around minor bodies) or when WH's
 * Kepler-drift assumption breaks down. Much cheaper than IAS15 per
 * step; the trade-off is that absolute precision is order-2 versus
 * IAS15's machine-precision.
 */
#ifndef K26ASTRO_GRAV_VERLET_H
#define K26ASTRO_GRAV_VERLET_H

#include "k26astro_grav/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

int k26astro_grav_step_verlet(K26AstroGravState *state, double dt);

#ifdef __cplusplus
}
#endif

#endif
