/* world.c — K26AstroWorld lifecycle + body collection.
 *
 * Composes K26AstroGravState under a single opaque handle.
 * Body storage is owned by the world (the world's reallocating buffer
 * is what grav_state's `bodies` pointer points at). */
#include "k26astro_rt/world.h"
#include "k26astro_rt/scheduler.h"
#include "k26astro_rt/referenced.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/pos.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_grav/perturb.h"
#include "k26astro_grav/forces.h"
#include "world_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* ---- IMU non-grav-accel capture ------------------------------- *
 *
 * Context owned by the world (heap; freed at destroy). Borrowed by
 * the gravity-state perturb list via k26astro_grav_register_perturb. */
typedef struct {
    struct K26AstroVehicle *vehicle;
    int                     body_idx;
    uint64_t                vehicle_gen_at_register;
} K26AstroWorldImuCacheCtx_;

static void rt_imu_accel_capture_(const K26AstroGravState *state,
                                  const K26AstroGravView *view,
                                  K26V3 *accel_out, void *ctx_v)
{
    (void)state;
    K26AstroWorldImuCacheCtx_ *ctx = (K26AstroWorldImuCacheCtx_ *)ctx_v;
    if (!ctx || !ctx->vehicle || !accel_out || !view) return;
    if (ctx->body_idx < 0 || ctx->body_idx >= view->n) return;
    /* Stale-owner gate. */
    if (ctx->vehicle_gen_at_register !=
        k26astro_vehicle_generation(ctx->vehicle)) return;
    k26astro_vehicle_set_last_non_grav_accel_inertial(
        ctx->vehicle, accel_out[ctx->body_idx]);
}

K26AstroWorld *k26astro_world_create(K26AstroWorldMode  mode,
                                      K26AstroCoordsMode coord_mode)
{
    if (mode != K26ASTRO_MODE_FAST &&
        mode != K26ASTRO_MODE_PORTABLE &&
        mode != K26ASTRO_MODE_REFERENCED) {
        return NULL;
    }
    if (coord_mode != K26ASTRO_COORDS_SECTOR_GRID &&
        coord_mode != K26ASTRO_COORDS_Q64_64) {
        return NULL;
    }

    K26AstroWorld *world = (K26AstroWorld *)calloc(1, sizeof(*world));
    if (!world) return NULL;

    world->mode       = mode;
    world->coord_mode = coord_mode;

    k26astro_fpu_save_and_pin(&world->fpu);

    /* Empty grav state. We open-code the init (grav_state_init rejects
     * NULL bodies / zero count) — we'll point `bodies` at our growing
     * buffer as bodies are added. */
    world->grav.bodies      = NULL;
    world->grav.n_bodies    = 0;
    world->grav.integrator  = K26ASTRO_INTEGRATOR_WH;
    world->grav.softening   = 0.0;
    world->grav.t           = k26astro_epoch_j2000_tt();
    world->grav.dt_last     = 0.0;
    world->grav.use_j2      = 0;
    world->grav.use_srp     = 0;
    world->grav.use_gr_ppn1 = 0;
    world->grav.perturbs    = NULL;
    world->grav.ias15_carry = NULL;
    world->grav.wh_carry    = NULL;

    /* Frame registry: leave next_user_id at the user-base; per-world
     * user frames live in `world->frames` (the libk26astro_core
     * registry is process-global, so we also push registrations
     * there when k26astro_world_register_frame is called). */
    world->n_frames     = 0;
    world->next_user_id = K26A_FRAME_USER_BASE;

    /* MERCURIUS Rein-Tamayo 2019 §3.1 defaults. */
    world->mercurius_hill_factor  = 3.0;
    world->mercurius_outer_factor = 5.0;

    world->observer_mode = K26ASTRO_OBS_ASTROMETRIC;
    world->atmos         = NULL;
    world->snapshot_id   = 0;
    world->user          = NULL;
    world->ref_ctx       = NULL;

    if (k26astro_rt_scheduler_init(world) != 0) {
        k26astro_grav_state_destroy(&world->grav);
        k26astro_fpu_restore(&world->fpu);
        free(world);
        return NULL;
    }

    return world;
}

void k26astro_world_destroy(K26AstroWorld *world)
{
    if (!world) return;
    /* Flush + close the referenced-mode log first so any pending
     * step records reach disk before scheduler teardown. */
    if (world->ref_ctx) (void)k26astro_world_ref_log_close(world);
    k26astro_rt_scheduler_destroy(world);

    free(world->encounters);
    world->encounters     = NULL;
    world->n_encounters   = 0;
    world->cap_encounters = 0;

    /* Vehicles are non-owning; just drop the pointer array. */
    free(world->vehicles);
    world->vehicles     = NULL;
    world->n_vehicles   = 0;
    world->cap_vehicles = 0;

    /* Inertial-measurement-unit capture contexts. The perturb list
     * (freed below as part of grav_state_destroy) holds borrowed
     * pointers; we free the underlying storage here. */
    for (int i = 0; i < world->n_imu_cache; i++) {
        free(world->imu_cache_ctxs[i]);
    }
    free(world->imu_cache_ctxs);
    world->imu_cache_ctxs = NULL;
    world->n_imu_cache    = 0;
    world->cap_imu_cache  = 0;

    /* Grav holds a pointer into our body buffer; tear down its
     * carry buffers first. */
    K26AstroBody *bodies = world->grav.bodies;
    k26astro_grav_state_destroy(&world->grav);
    free(bodies);

    k26astro_fpu_restore(&world->fpu);
    free(world);
}

/* Body collection ------------------------------------------------- */

int k26astro_world_add_body(K26AstroWorld *world, K26AstroBody b)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    int n = world->grav.n_bodies;
    K26AstroBody *grown = (K26AstroBody *)realloc(
        world->grav.bodies, (size_t)(n + 1) * sizeof(K26AstroBody));
    if (!grown) return -K26ASTRO_RT_E_OOM;
    grown[n] = b;
    world->grav.bodies   = grown;
    world->grav.n_bodies = n + 1;
    return n;
}

int k26astro_world_remove_body(K26AstroWorld *world, int idx)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    int n = world->grav.n_bodies;
    if (idx < 0 || idx >= n) return -K26ASTRO_RT_E_BAD_ARG;
    /* Tail-swap. */
    if (idx != n - 1) {
        world->grav.bodies[idx] = world->grav.bodies[n - 1];
    }
    world->grav.n_bodies = n - 1;
    /* Don't shrink the allocation — capacity sticks around. */
    return K26ASTRO_RT_OK;
}

int k26astro_world_body_count(const K26AstroWorld *world)
{
    if (!world) return 0;
    return world->grav.n_bodies;
}

int k26astro_world_find_body(const K26AstroWorld *world, const char *name)
{
    if (!world || !name) return -1;
    for (int i = 0; i < world->grav.n_bodies; i++) {
        if (strcasecmp(world->grav.bodies[i].name, name) == 0) return i;
    }
    return -1;
}

void k26astro_world_for_each_body(const K26AstroWorld *world,
                                   K26AstroBodyVisitorFn visit,
                                   void *user)
{
    if (!world || !visit) return;
    for (int i = 0; i < world->grav.n_bodies; i++) {
        visit(&world->grav.bodies[i], i, user);
    }
}

K26AstroBody *k26astro_world_body_at(K26AstroWorld *world, int idx)
{
    if (!world) return NULL;
    if (idx < 0 || idx >= world->grav.n_bodies) return NULL;
    return &world->grav.bodies[idx];
}

const char *k26astro_body_name(const K26AstroBody *b)
{
    return b ? b->name : NULL;
}

/* Stepping --------------------------------------------------------- */

int k26astro_world_step(K26AstroWorld *world, double wallclock_dt_s)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    if (!(wallclock_dt_s >= 0.0)) return -K26ASTRO_RT_E_BAD_ARG;
    k26astro_rt_ref_emit_step_begin(world, wallclock_dt_s);
    (void)k26tick_advance(world->tick, wallclock_dt_s);
    /* World time after the step — seconds-past-J2000 (TDB). */
    double t_s = (double)world->grav.t.days_since_J2000 * 86400.0
               + world->grav.t.seconds_of_day;
    k26astro_rt_ref_emit_step_end(world, t_s);
    return K26ASTRO_RT_OK;
}

int k26astro_world_body_step(K26AstroWorld *world, int body_idx,
                              double dt)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    if (body_idx < 0 || body_idx >= world->grav.n_bodies)
        return -K26ASTRO_RT_E_BAD_ARG;
    if (!isfinite(dt)) return -K26ASTRO_RT_E_BAD_ARG;

    K26AstroBody *body = &world->grav.bodies[body_idx];
    int parent_idx = body->parent_body_idx;
    /* Heliocentric convention: when no SOI parent set, fall back to
     * body 0 (typically the Sun). */
    if (parent_idx < 0) parent_idx = 0;
    if (parent_idx >= world->grav.n_bodies)
        return -K26ASTRO_RT_E_BAD_ARG;
    if (parent_idx == body_idx) return -K26ASTRO_RT_E_BAD_ARG;

    const K26AstroBody *parent = &world->grav.bodies[parent_idx];
    if (!(parent->gm > 0.0)) return -K26ASTRO_RT_E_BAD_ARG;

    /* Body state in the parent's instantaneous rest frame. Position
     * via k26astro_pos_sub preserves sector-grid precision (the
     * subtraction lives in fixed-point sector arithmetic; result is
     * a K26V3 in meters). */
    K26V3 rel_pos = k26astro_pos_sub(&body->pos, &parent->pos);
    K26V3 rel_vel = { body->vel.x - parent->vel.x,
                      body->vel.y - parent->vel.y,
                      body->vel.z - parent->vel.z };

    K26V3 new_rel_pos, new_rel_vel;
    int rc = k26astro_kepler_propagate(&new_rel_pos, &new_rel_vel,
                                        rel_pos, rel_vel,
                                        parent->gm, dt, 32);
    if (rc != 0) return -K26ASTRO_RT_E_INTEGRATOR;

    /* Rebuild world-frame pos from parent + delta to keep precision. */
    body->pos = parent->pos;
    k26astro_pos_add(&body->pos, new_rel_pos);
    body->vel.x = parent->vel.x + new_rel_vel.x;
    body->vel.y = parent->vel.y + new_rel_vel.y;
    body->vel.z = parent->vel.z + new_rel_vel.z;
    return K26ASTRO_RT_OK;
}

int k26astro_world_now(const K26AstroWorld *world, K26AstroEpoch *out)
{
    if (!world || !out) return -K26ASTRO_RT_E_NULL;
    *out = world->grav.t;
    return K26ASTRO_RT_OK;
}

/* Configuration --------------------------------------------------- */

int k26astro_world_set_ephem(K26AstroWorld *world, K26AstroEphem *ephem)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    world->ephem = ephem;
    return K26ASTRO_RT_OK;
}

int k26astro_world_set_observer_mode(K26AstroWorld *world,
                                      K26AstroObserverMode mode)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    if (mode > K26ASTRO_OBS_TOPOCENTRIC) return -K26ASTRO_RT_E_BAD_ARG;
    world->observer_mode = mode;
    return K26ASTRO_RT_OK;
}

int k26astro_world_set_atmos(K26AstroWorld *world,
                              struct K26AstroAtmos *atmos)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    world->atmos = atmos;
    return K26ASTRO_RT_OK;
}

int k26astro_world_set_mercurius_factors(K26AstroWorld *world,
                                          double y_inner, double y_outer)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    if (!(y_inner > 0.0) || !(y_outer > y_inner))
        return -K26ASTRO_RT_E_BAD_ARG;
    world->mercurius_hill_factor  = y_inner;
    world->mercurius_outer_factor = y_outer;
    return K26ASTRO_RT_OK;
}

K26AstroGravState *k26astro_world_grav(K26AstroWorld *world)
{
    return world ? &world->grav : NULL;
}

/* Vehicle registry ------------------------------------------------ */

int k26astro_world_register_vehicle(K26AstroWorld *world,
                                    struct K26AstroVehicle *v)
{
    if (!world || !v) return -K26ASTRO_RT_E_NULL;

    /* Idempotent on duplicates. */
    for (int i = 0; i < world->n_vehicles; i++) {
        if (world->vehicles[i] == v) return K26ASTRO_RT_OK;
    }

    if (world->n_vehicles >= world->cap_vehicles) {
        int new_cap = world->cap_vehicles ? (world->cap_vehicles * 2) : 4;
        struct K26AstroVehicle **grown = (struct K26AstroVehicle **)realloc(
            world->vehicles,
            (size_t)new_cap * sizeof(struct K26AstroVehicle *));
        if (!grown) return -K26ASTRO_RT_E_OOM;
        world->vehicles     = grown;
        world->cap_vehicles = new_cap;
    }
    world->vehicles[world->n_vehicles++] = v;
    return K26ASTRO_RT_OK;
}

void k26astro_world_unregister_vehicle(K26AstroWorld *world,
                                       struct K26AstroVehicle *v)
{
    if (!world || !v) return;
    for (int i = 0; i < world->n_vehicles; i++) {
        if (world->vehicles[i] == v) {
            /* Tail-swap; registry order isn't load-bearing for
             * commit_mass_step (per-vehicle accumulators are
             * independent scalars). */
            int last = world->n_vehicles - 1;
            if (i != last) world->vehicles[i] = world->vehicles[last];
            world->vehicles[last] = NULL;
            world->n_vehicles--;
            return;
        }
    }
}

int k26astro_world_enable_imu_accel_cache(K26AstroWorld *world,
                                          struct K26AstroVehicle *v,
                                          int body_idx)
{
    if (!world || !v) return -K26ASTRO_RT_E_NULL;
    if (body_idx < 0 || body_idx >= world->grav.n_bodies)
        return -K26ASTRO_RT_E_BAD_ARG;

    K26AstroWorldImuCacheCtx_ *ctx = (K26AstroWorldImuCacheCtx_ *)
        calloc(1, sizeof(*ctx));
    if (!ctx) return -K26ASTRO_RT_E_OOM;
    ctx->vehicle                 = v;
    ctx->body_idx                = body_idx;
    ctx->vehicle_gen_at_register = k26astro_vehicle_generation(v);

    if (world->n_imu_cache >= world->cap_imu_cache) {
        int new_cap = world->cap_imu_cache ? (world->cap_imu_cache * 2) : 4;
        void **grown = (void **)realloc(world->imu_cache_ctxs,
                                        (size_t)new_cap * sizeof(void *));
        if (!grown) {
            free(ctx);
            return -K26ASTRO_RT_E_OOM;
        }
        world->imu_cache_ctxs = grown;
        world->cap_imu_cache  = new_cap;
    }
    world->imu_cache_ctxs[world->n_imu_cache++] = ctx;

    int rc = k26astro_grav_register_perturb(&world->grav,
                                            rt_imu_accel_capture_, ctx);
    if (rc != 0) {
        /* The ctx storage stays owned by the world (will be freed at
         * destroy) even though the perturb registration failed —
         * keeping the array invariant simple. */
        return -K26ASTRO_RT_E_OOM;
    }
    return K26ASTRO_RT_OK;
}
