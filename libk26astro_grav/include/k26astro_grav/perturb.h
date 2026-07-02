/* k26astro_grav/perturb.h — perturbation registry + built-ins.
 *
 * A perturbation is an additive acceleration on top of point-mass
 * Newtonian gravity. Built-ins:
 *
 *   J2          - oblateness for bodies with body.j2 != 0.
 *                 Effective within ~10 body radii.
 *
 *   SRP         - solar radiation pressure on bodies with non-zero
 *                 srp_area_over_mass. Conical shadow geometry
 *                 detects occlusion by user-marked occluder bodies
 *                 (Earth, Moon by default).
 *
 *   GR PPN-1    - Schwarzschild post-Newtonian correction in the
 *                 Einstein-Infeld-Hoffmann form (Will 1993 Eq 8.62).
 *                 Produces the ~43"/century Mercury perihelion
 *                 advance with the right sign.
 *
 *   drag (stub) - placeholder hook for libk26astro_atmos.
 *                 Returns zero acceleration; the hook exists so the
 *                 force pipeline can adopt atmospheric drag without
 *                 ABI churn.
 *
 * User-registered perturbations are called in registration order.
 * The registry is owned by K26AstroGravState; freed on destroy. */
#ifndef K26ASTRO_GRAV_PERTURB_H
#define K26ASTRO_GRAV_PERTURB_H

#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Perturbation function: take state + view, ADD acceleration into
 * accel_out (do not overwrite — the registry is additive). */
typedef void (*K26AstroPerturbFn)(const K26AstroGravState *state,
                                   const K26AstroGravView  *view,
                                   K26V3 *accel_out,
                                   void  *ctx);

struct K26AstroPerturbList {
    K26AstroPerturbFn *fns;
    void             **ctxs;
    int                count;
    int                capacity;
};

/* Register a user perturbation. Returns 0 on success. Lifetime
 * of `ctx` is the caller's responsibility. */
int k26astro_grav_register_perturb(K26AstroGravState *state,
                                    K26AstroPerturbFn fn,
                                    void *ctx);

/* Built-in perturbation entry points (called by accel_total when the
 * corresponding use_* flag is on). Exposed in case callers want to
 * invoke them directly outside the integrator pipeline. */
void k26astro_perturb_j2     (const K26AstroGravState *state,
                              const K26AstroGravView *view,
                              K26V3 *accel_out, void *ctx);

void k26astro_perturb_srp    (const K26AstroGravState *state,
                              const K26AstroGravView *view,
                              K26V3 *accel_out, void *ctx);

void k26astro_perturb_gr_ppn1(const K26AstroGravState *state,
                              const K26AstroGravView *view,
                              K26V3 *accel_out, void *ctx);

/* SRP shadow detector: returns 1 if body `i` is in the shadow of any
 * occluder (currently bodies whose kind == K26ASTRO_BODY_PLANET or
 * K26ASTRO_BODY_MOON within ~1 AU of the Sun). */
int  k26astro_srp_shadow_test(const K26AstroGravView *view,
                              int sun_idx, int body_idx);

#ifdef __cplusplus
}
#endif

#endif
