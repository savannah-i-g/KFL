/* test_vehicle_perturb_capture.c — inertial-measurement-unit cache
 * wiring: k26astro_world_enable_imu_accel_cache installs a perturb
 * callback that snapshots accel_out[body_idx] into the vehicle's
 * last_non_grav_accel_inertial field.
 *
 * Scenarios:
 *   (1) Register a vehicle on body 1 and enable the cache. Install
 *       a synthetic perturb that adds {1.0, 2.0, 3.0} to
 *       accel_out[1]. Step the world. Verify the vehicle's cache
 *       reads {1.0, 2.0, 3.0} exactly (after registration order
 *       guarantees the capture runs last).
 *   (2) Bad inputs: NULL world, NULL vehicle, out-of-range body_idx
 *       return non-zero without crash.
 *   (3) Reuse: enabling twice for the same vehicle is allowed; the
 *       last-registered capture wins. */

#include "k26astro_rt/world.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/perturb.h"
#include "k26astro_grav/forces.h"
#include "k26m3d.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

/* Synthetic perturb: ADDS a constant accel to accel_out[body_idx]. */
typedef struct {
    int    body_idx;
    K26V3  accel;
} ConstAccelCtx;

static void const_accel_perturb_(const K26AstroGravState *state,
                                  const K26AstroGravView  *view,
                                  K26V3 *accel_out,
                                  void  *ctx_v)
{
    (void)state;
    ConstAccelCtx *ctx = (ConstAccelCtx *)ctx_v;
    if (!ctx || !accel_out || !view) return;
    if (ctx->body_idx < 0 || ctx->body_idx >= view->n) return;
    accel_out[ctx->body_idx].x += ctx->accel.x;
    accel_out[ctx->body_idx].y += ctx->accel.y;
    accel_out[ctx->body_idx].z += ctx->accel.z;
}

int main(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_PORTABLE,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    /* Body 0: massive primary. */
    K26AstroBody primary;
    k26astro_body_init(&primary);
    strncpy(primary.name, "primary", sizeof primary.name - 1);
    primary.mass = 1.989e30;
    primary.gm   = 1.32712440018e20;
    int idx_primary = k26astro_world_add_body(w, primary);
    assert(idx_primary == 0);

    /* Body 1: spacecraft body. Position left at default zero — the
     * test exercises the perturb-capture flow, not the integrator's
     * Newtonian gravity reduction. */
    K26AstroBody craft;
    k26astro_body_init(&craft);
    strncpy(craft.name, "spacecraft", sizeof craft.name - 1);
    craft.mass = 1000.0;
    int idx_craft = k26astro_world_add_body(w, craft);
    assert(idx_craft == 1);

    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v);
    k26astro_vehicle_bind_body(v, k26astro_world_body_at(w, idx_craft));
    k26astro_vehicle_set_dry_mass(v, 1000.0);

    /* ---- Scenario 2: bad inputs ----------------------------- */
    assert(k26astro_world_enable_imu_accel_cache(NULL, v, 1)   < 0);
    assert(k26astro_world_enable_imu_accel_cache(w, NULL, 1)   < 0);
    assert(k26astro_world_enable_imu_accel_cache(w, v, -1)     < 0);
    assert(k26astro_world_enable_imu_accel_cache(w, v, 999)    < 0);

    /* ---- Scenario 1: const-accel + capture ------------------ *
     *
     * Register the synthetic perturb FIRST so it runs first; then
     * enable the capture which registers a perturb LAST. The
     * capture observes the running accel_out sum, which by then
     * includes the synthetic contribution. */
    K26AstroGravState *grav = k26astro_world_grav(w);
    assert(grav);

    ConstAccelCtx pctx = { .body_idx = 1,
                           .accel = { 1.0, 2.0, 3.0 } };
    int prc = k26astro_grav_register_perturb(grav,
                                              const_accel_perturb_,
                                              &pctx);
    assert(prc == 0);

    int rc = k26astro_world_register_vehicle(w, v);
    assert(rc == 0);
    rc = k26astro_world_enable_imu_accel_cache(w, v, 1);
    assert(rc == 0);

    /* Pre-step: cache reads zero (vehicle_new calloc'd it). */
    K26V3 pre = k26astro_vehicle_last_non_grav_accel_inertial(v);
    assert(pre.x == 0.0 && pre.y == 0.0 && pre.z == 0.0);

    /* Step the world. The orbit channel drives k26astro_grav_step,
     * which calls k26astro_grav_accel_total, which invokes the
     * perturb list (synthetic first, capture second). */
    rc = k26astro_world_step(w, 1.0);
    assert(rc == 0);

    K26V3 cached = k26astro_vehicle_last_non_grav_accel_inertial(v);
    /* The capture writes the running accel_out sum. With only the
     * synthetic perturb contributing to body 1 (no J2/SRP/GR by
     * default), the captured value must equal the synthetic accel. */
    assert(cached.x == 1.0);
    assert(cached.y == 2.0);
    assert(cached.z == 3.0);

    /* ---- Scenario 3: reuse / double-enable ----------------- */
    rc = k26astro_world_enable_imu_accel_cache(w, v, 1);
    assert(rc == 0);
    rc = k26astro_world_step(w, 1.0);
    assert(rc == 0);
    cached = k26astro_vehicle_last_non_grav_accel_inertial(v);
    /* Same accel; both captures observe the same sum. */
    assert(cached.x == 1.0);
    assert(cached.y == 2.0);
    assert(cached.z == 3.0);

    /* ---- Cleanup --------------------------------------------- */
    k26astro_vehicle_destroy(v);
    k26astro_world_destroy(w);

    printf("test_vehicle_perturb_capture: OK\n");
    return 0;
}
