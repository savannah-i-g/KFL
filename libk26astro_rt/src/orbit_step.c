/* orbit_step.c — orbit-channel callback (paper-faithful MERCURIUS).
 *
 * Runs every wallclock advance with the actual elapsed dt. Drives
 * the world's integrator (WH, IAS15, Verlet, RK4, RK45, or MERCURIUS).
 *
 * Paper-faithful MERCURIUS, Rein-Tamayo 2019 eq. 12-14:
 *   1. Detect close-encounter pairs (k26astro_mercurius_detect).
 *      Each pair carries a K(y) weight computed at detect time.
 *   2. If no encounters or active integrator is not WH/Verlet:
 *      single full-force step (the integrator handles dynamics
 *      uniformly).
 *   3. Otherwise (paper-faithful split):
 *        a. Set state->mercurius = { FAR, weights, n }.
 *        b. Run the chosen outer integrator (WH or Verlet) over
 *           dt; it sees only the (1-K)-weighted portion of each
 *           encounter pair plus full-force on non-encounter pairs.
 *           This handles the smooth bulk dynamics.
 *        c. Set state->mercurius = { NEAR, weights, n }.
 *        d. Run IAS15 over the same dt — it sees only the
 *           K-weighted portion of each encounter pair (zero for
 *           non-encounter pairs). This handles the encounter
 *           dynamics with adaptive precision.
 *        e. Clear state->mercurius.
 *
 *      The two integrators contribute additively to position +
 *      velocity because (a) they integrate disjoint pieces of the
 *      total acceleration (Σ((1-K)+K) = Σ identity) and (b) the
 *      MERCURIUS context applies its K filter inside accel_total,
 *      so each integrator's internal kick/drift sees the right
 *      force field. Position and velocity updates accumulate in
 *      place — no separate delta-summation buffer needed. */
#include "encounter_internal.h"

#include "k26astro_grav/grav.h"
#include "k26astro_grav/close_encounter.h"
#include "k26astro_grav/forces.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_vehicle/vehicle.h"

#include <stdlib.h>

/* Walk the world's vehicle registry and consume each per-substep
 * mass accumulator. Propulsion thrust callbacks add into the
 * accumulator during accel_total; this is the closing step that
 * propagates accumulated dot_m into vehicle.basic_mass_kg and the
 * bound body's mass/GM via k26astro_body_set_mass. */
static void commit_vehicle_mass_(K26AstroWorld *world, double dt_s)
{
    for (int i = 0; i < world->n_vehicles; i++) {
        K26AstroVehicle *v = world->vehicles[i];
        if (v) k26astro_vehicle_commit_mass_step(v, dt_s);
    }
}

/* Build a contiguous K26AstroPairWeight array from world->encounters.
 * Returns a heap-allocated buffer the caller must free; *out_n holds
 * the count. Returns NULL on n=0 or allocation failure. */
static K26AstroPairWeight *build_pair_weights_(const K26AstroWorld *world,
                                                int *out_n)
{
    int n = world->n_encounters;
    *out_n = 0;
    if (n <= 0) return NULL;
    K26AstroPairWeight *w = (K26AstroPairWeight *)
        malloc((size_t)n * sizeof(K26AstroPairWeight));
    if (!w) return NULL;
    for (int k = 0; k < n; k++) {
        w[k].i = world->encounters[k].i;
        w[k].j = world->encounters[k].j;
        w[k].k_weight = world->encounters[k].k_weight;
    }
    *out_n = n;
    return w;
}

void k26astro_rt_orbit_step_cb(double dt_s, void *user)
{
    K26AstroWorld *world = (K26AstroWorld *)user;
    if (!world) return;
    if (!(dt_s > 0.0)) return;

    /* Encounter detection — populates world->encounters with K
     * weights baked in. */
    int n_enc = k26astro_mercurius_detect(world);

    K26AstroIntegrator base = world->grav.integrator;
    int do_split = (n_enc > 0)
        && (base == K26ASTRO_INTEGRATOR_WH ||
            base == K26ASTRO_INTEGRATOR_VERLET);

    if (!do_split) {
        /* Standard single-integrator step. */
        (void)k26astro_grav_step(&world->grav, dt_s);
        commit_vehicle_mass_(world, dt_s);
        return;
    }

    /* Paper-faithful MERCURIUS split (Rein-Tamayo 2019 eq. 12-14). */
    int n_w = 0;
    K26AstroPairWeight *w = build_pair_weights_(world, &n_w);
    if (!w) {
        /* Fallback: full single-integrator step on alloc failure. */
        (void)k26astro_grav_step(&world->grav, dt_s);
        commit_vehicle_mass_(world, dt_s);
        return;
    }

    K26AstroMercuriusContext far_ctx  = {
        .mode = K26ASTRO_MERCURIUS_FAR, .pair_weights = w, .n_pair_weights = n_w };
    K26AstroMercuriusContext near_ctx = {
        .mode = K26ASTRO_MERCURIUS_NEAR, .pair_weights = w, .n_pair_weights = n_w };

    /* Step 1: outer (WH or Verlet) on FAR. */
    world->grav.mercurius = &far_ctx;
    (void)k26astro_grav_step(&world->grav, dt_s);

    /* Step 2: IAS15 sub-step on NEAR. Switch the integrator for
     * the inner pass, then restore. */
    world->grav.mercurius = &near_ctx;
    (void)k26astro_grav_set_integrator(&world->grav,
                                         K26ASTRO_INTEGRATOR_IAS15);
    (void)k26astro_grav_step(&world->grav, dt_s);
    (void)k26astro_grav_set_integrator(&world->grav, base);

    world->grav.mercurius = NULL;
    free(w);

    /* MERCURIUS split: commit mass once per outer substep (FAR pass).
     * The NEAR sub-step is internal to IAS15 and shouldn't double-count
     * the dot_m accumulator. */
    commit_vehicle_mass_(world, dt_s);
}
