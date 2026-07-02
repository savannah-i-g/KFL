/* k26astro_grav/perturb_outer_planets.h — outer-planet on-rails
 * perturbation interface.
 *
 * Adds GM·r̂/r² acceleration from a list of ephemeris-driven outer
 * planets (Jupiter, Saturn, optionally Uranus + Neptune) to every
 * body in the integrated view. The outer planets are NOT integrated
 * by the K26 N-body engine; they are queried from a supplied
 * K26AstroEphem (typically the DE441 inner-SS kernel extended with
 * outer planets, or a separate outer-planet SPK) at the integrator's
 * current sim epoch.
 *
 * Use case: dropping the inner-planet 100-year integration residual
 * from ~1e10-1e11 m (no outer planets) to ~100-1000 km (with
 * Jupiter + Saturn on rails). The dominant outer-planet secular
 * forcing comes from Jupiter (~318 Earth masses) and Saturn
 * (~95 Earth masses) at their 5.2 / 9.5 AU heliocentric distances. */
#ifndef K26ASTRO_GRAV_PERTURB_OUTER_PLANETS_H
#define K26ASTRO_GRAV_PERTURB_OUTER_PLANETS_H

#include "k26astro_grav/grav.h"
#include "k26astro_grav/perturb.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of outer planets supported in a single context. */
#define K26ASTRO_MAX_OUTER_PLANETS 8

/* Forward declaration; full type defined in libk26astro_ephem. */
typedef struct K26AstroEphem K26AstroEphem;

typedef struct K26AstroOuterPlanetsCtx {
    K26AstroEphem *ephem;
    int            n_outer;
    int            naif_ids[K26ASTRO_MAX_OUTER_PLANETS];
    double         gms[K26ASTRO_MAX_OUTER_PLANETS];   /* m³/s² */
} K26AstroOuterPlanetsCtx;

/* Perturbation function compatible with K26AstroPerturbFn. Reads
 * outer-planet positions from `ctx->ephem` at the integrator's
 * `state->t` and applies GM·r̂/r² to every body in the view's
 * accel_out array. Bodies whose ephem query returns the zero-state
 * sentinel are skipped silently. */
void k26astro_perturb_outer_planets(const K26AstroGravState *state,
                                    const K26AstroGravView  *view,
                                    K26V3 *accel_out, void *ctx);

/* Convenience registrar: register the outer-planets perturbation
 * with `state`'s perturbation registry. The caller retains
 * ownership of `ctx`; `state` does not copy or free it.
 *
 * Returns 0 on success; -1 on invalid args (NULL state/ctx, or
 * n_outer out of [1, K26ASTRO_MAX_OUTER_PLANETS]).
 *
 * Typical use:
 *   K26AstroOuterPlanetsCtx opc = {
 *       .ephem    = my_ephem,
 *       .n_outer  = 2,
 *       .naif_ids = { 5, 6 },       // Jupiter barycentre, Saturn barycentre
 *       .gms      = { 1.26686534e17, 3.7931187e16 },
 *   };
 *   k26astro_grav_enable_outer_planets(&grav_state, &opc);
 */
int k26astro_grav_enable_outer_planets(K26AstroGravState *state,
                                       K26AstroOuterPlanetsCtx *ctx);

/* KFL-friendly convenience wrapper: registers a Jupiter + Saturn
 * outer-planets perturbation with default NAIF ids (5, 6) and
 * IAU 2015 GM values. The internal ctx is a process-scoped
 * singleton; multi-world programs needing distinct configs
 * (different ephem files, more outer planets) must use
 * k26astro_grav_enable_outer_planets() with their own ctx
 * allocation. */
int k26astro_grav_enable_outer_planets_default(K26AstroGravState *state,
                                                K26AstroEphem    *ephem);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_GRAV_PERTURB_OUTER_PLANETS_H */
