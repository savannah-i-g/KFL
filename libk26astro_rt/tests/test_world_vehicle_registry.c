/* test_world_vehicle_registry.c — world-side vehicle registry +
 * commit_mass_step wiring.
 *
 * Verifies the registry that orbit_step.c walks at substep close:
 *   (1) Register a vehicle on a world; world_step's orbit channel
 *       drives the registry walk so the vehicle's mass accumulator
 *       commits each substep.
 *   (2) A synthetic perturb that adds a constant dot_m = -1 kg/s
 *       must shrink vehicle.basic_mass_kg as expected over a 10 s
 *       integration window.
 *   (3) Unregister then re-step: mass stops decreasing.
 *   (4) Duplicate register is idempotent.
 *   (5) NULL inputs are safe. */

#include "k26astro_rt/world.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/perturb.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Synthetic perturb: pulls a vehicle pointer out of ctx and adds
 * a constant dot_m to its accumulator. No acceleration written. */
typedef struct {
    K26AstroVehicle *v;
    double           dot_m;     /* kg/s, negative for mass leaving */
} TestPerturbCtx;

static void test_dot_m_perturb_(const K26AstroGravState *state,
                                const K26AstroGravView  *view,
                                K26V3 *accel_out,
                                void  *ctx_v)
{
    (void)state; (void)view; (void)accel_out;
    TestPerturbCtx *ctx = ctx_v;
    if (ctx && ctx->v) {
        k26astro_vehicle_mass_accum_add(ctx->v, ctx->dot_m);
    }
}

int main(void)
{
    /* ---- Build a minimal world + body + vehicle. ------------- */
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_PORTABLE,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    /* Body 0: massive primary (Sun-like, ~1e30 kg). */
    K26AstroBody primary;
    k26astro_body_init(&primary);
    strncpy(primary.name, "primary", sizeof primary.name - 1);
    primary.mass = 1.989e30;
    primary.gm   = 1.32712440018e20;
    int idx_primary = k26astro_world_add_body(w, primary);
    assert(idx_primary == 0);

    /* Body 1: spacecraft body. Mass set via vehicle later; bound to
     * vehicle for mass-step propagation. */
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

    /* ---- Scenario 4: NULL safety. ---------------------------- */
    assert(k26astro_world_register_vehicle(NULL, v) < 0);
    assert(k26astro_world_register_vehicle(w, NULL) < 0);
    k26astro_world_unregister_vehicle(NULL, v);  /* no crash */
    k26astro_world_unregister_vehicle(w, NULL);  /* no crash */

    /* ---- Scenario 1+2: register + step + verify mass shrinks. */
    int rc = k26astro_world_register_vehicle(w, v);
    assert(rc == 0);

    /* Register synthetic perturb that drives mass_accum_add. */
    TestPerturbCtx ctx = { .v = v, .dot_m = -1.0 };  /* -1 kg/s */
    K26AstroGravState *grav = k26astro_world_grav(w);
    assert(grav);
    int prc = k26astro_grav_register_perturb(grav,
                                              test_dot_m_perturb_,
                                              &ctx);
    (void)prc;
    /* Note: the perturb registry fires inside accel_total, called
     * by k26astro_grav_step. The synthetic perturb writes no accel
     * (single body still drifts as Keplerian); it only mutates the
     * accumulator, which orbit_step.c commits after the grav step. */

    double mass_initial = k26astro_vehicle_mass_now(v);
    assert(mass_initial == 1000.0);

    /* Advance 10 s of simulated time at 1 Hz default orbit rate. */
    for (int i = 0; i < 10; i++) {
        rc = k26astro_world_step(w, 1.0);
        assert(rc == 0);
    }

    double mass_after_10s = k26astro_vehicle_mass_now(v);
    /* Each orbit-channel tick should fire one perturb call adding
     * dot_m=-1 kg/s to the accumulator; commit_mass_step then
     * applies dm = -1 * dt where dt is the orbit tick interval.
     * Total expected: roughly -10 kg over 10 s of wallclock at any
     * sane orbit rate. The exact value depends on the tick schedule,
     * but the direction + bounded magnitude check is the gate. */
    assert(mass_after_10s < mass_initial);
    assert(mass_after_10s > mass_initial - 100.0);

    /* ---- Scenario 4 cont.: duplicate register is idempotent. - */
    rc = k26astro_world_register_vehicle(w, v);
    assert(rc == 0);  /* still ok, no leak */

    /* ---- Scenario 3: unregister + step + mass stays. --------- */
    double mass_at_unregister = k26astro_vehicle_mass_now(v);
    k26astro_world_unregister_vehicle(w, v);

    for (int i = 0; i < 10; i++) {
        rc = k26astro_world_step(w, 1.0);
        assert(rc == 0);
    }

    /* The perturb still fires, still adds dot_m, but
     * commit_mass_step no longer runs against this vehicle. So
     * basic_mass_kg stays at mass_at_unregister; mass_accum
     * accumulates but isn't consumed. */
    double mass_after_unregister = k26astro_vehicle_mass_now(v);
    assert(mass_after_unregister == mass_at_unregister);

    /* ---- Cleanup. -------------------------------------------- */
    k26astro_vehicle_destroy(v);
    k26astro_world_destroy(w);

    printf("test_world_vehicle_registry: OK\n");
    return 0;
}
