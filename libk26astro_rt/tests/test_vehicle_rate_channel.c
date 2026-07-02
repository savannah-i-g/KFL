/* test_vehicle_rate_channel.c — per-vehicle sub-orbit-channel
 * cadence gate.
 *
 * Validates the wrapper around libk26tick that gates RCS-pulse /
 * engine-ignition transients on a per-vehicle predicate. The real
 * K26AstroVehicle struct is defined in libk26astro_vehicle; this
 * test only requires that the rate-channel layer carry the vehicle
 * handle without inspecting it, so it defines a stub `struct
 * K26AstroVehicle` inline. The matching forward declaration in
 * scheduler.h pins the API surface.
 *
 * Coverage:
 *   (1) Step count after world advance matches floor(elapsed · hz)
 *       within the libk26tick accumulator's ±1-step boundary.
 *   (2) Predicate-gated inactive vehicle: step count is flat while
 *       the predicate returns 0.
 *   (3) Live Hz change via _set_hz takes effect mid-run.
 *   (4) Destroy disables the slot — subsequent world advances do
 *       not call the step. */
#include "k26astro_rt/world.h"
#include "k26astro_rt/scheduler.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

/* Stub — the real struct lives in libk26astro_vehicle. The
 * rate-channel layer takes K26AstroVehicle by forward-declared
 * pointer and only forwards it to the caller's callbacks. */
struct K26AstroVehicle {
    int active;       /* 0 = rails, 1 = powered */
    int step_count;
    double last_dt;
};

static int vehicle_is_active_(const struct K26AstroVehicle *v, void *ctx)
{
    (void)ctx;
    return v ? v->active : 0;
}

static void vehicle_step_(struct K26AstroVehicle *v, double dt, void *ctx)
{
    (void)ctx;
    if (!v) return;
    v->step_count++;
    v->last_dt = dt;
}

int main(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    struct K26AstroVehicle veh = { .active = 1, .step_count = 0, .last_dt = 0.0 };

    /* ---- (1) Step count matches floor(elapsed · hz) ----------- */
    K26AstroVehicleRateChannel *ch =
        k26astro_rt_vehicle_rate_channel_new(w, &veh, 100.0,
                                              vehicle_is_active_,
                                              vehicle_step_, NULL);
    assert(ch != NULL);
    assert(fabs(k26astro_rt_vehicle_rate_channel_hz(ch) - 100.0) < 1e-9);
    assert(k26astro_rt_vehicle_rate_channel_active(ch) == 1);

    /* Advance world by 1.0 s in 10 wallclock chunks of 0.1 s. */
    for (int i = 0; i < 10; i++) {
        int rc = k26astro_world_step(w, 0.1);
        assert(rc == 0);
    }
    /* 100 Hz · 1 s = 100 steps; accumulator may carry ±1. */
    assert(veh.step_count >= 99 && veh.step_count <= 101);
    /* Per-step dt is 1/hz = 0.01 s. */
    assert(fabs(veh.last_dt - 0.01) < 1e-9);

    /* ---- (2) Inactive predicate gates the step ---------------- */
    int before = veh.step_count;
    veh.active = 0;
    for (int i = 0; i < 10; i++) {
        int rc = k26astro_world_step(w, 0.1);
        assert(rc == 0);
    }
    /* No steps should have fired while the predicate said inactive. */
    assert(veh.step_count == before);
    assert(k26astro_rt_vehicle_rate_channel_active(ch) == 0);

    /* Reactivate; steps resume. */
    veh.active = 1;
    int before2 = veh.step_count;
    for (int i = 0; i < 10; i++) {
        int rc = k26astro_world_step(w, 0.1);
        assert(rc == 0);
    }
    int fired = veh.step_count - before2;
    assert(fired >= 99 && fired <= 101);

    /* ---- (3) Live rate change ---------------------------------- */
    int rc = k26astro_rt_vehicle_rate_channel_set_hz(ch, 50.0);
    assert(rc == 0);
    assert(fabs(k26astro_rt_vehicle_rate_channel_hz(ch) - 50.0) < 1e-9);
    int before3 = veh.step_count;
    for (int i = 0; i < 10; i++) {
        int rc2 = k26astro_world_step(w, 0.1);
        assert(rc2 == 0);
    }
    int fired3 = veh.step_count - before3;
    /* 50 Hz · 1 s = 50 steps ± accumulator boundary. */
    assert(fired3 >= 49 && fired3 <= 51);
    assert(fabs(veh.last_dt - 0.02) < 1e-9);

    /* ---- (4) Destroy disables further dispatch ---------------- */
    k26astro_rt_vehicle_rate_channel_destroy(ch);
    int before4 = veh.step_count;
    for (int i = 0; i < 10; i++) {
        int rc2 = k26astro_world_step(w, 0.1);
        assert(rc2 == 0);
    }
    assert(veh.step_count == before4);

    /* ---- Null guards ------------------------------------------ */
    assert(k26astro_rt_vehicle_rate_channel_new(NULL, &veh, 100.0,
                                                  vehicle_is_active_,
                                                  vehicle_step_, NULL) == NULL);
    assert(k26astro_rt_vehicle_rate_channel_new(w, NULL, 100.0,
                                                  vehicle_is_active_,
                                                  vehicle_step_, NULL) == NULL);
    assert(k26astro_rt_vehicle_rate_channel_new(w, &veh, 100.0,
                                                  vehicle_is_active_,
                                                  NULL, NULL) == NULL);
    k26astro_rt_vehicle_rate_channel_destroy(NULL);
    assert(k26astro_rt_vehicle_rate_channel_set_hz(NULL, 100.0) != 0);
    assert(k26astro_rt_vehicle_rate_channel_hz(NULL) == 0.0);
    assert(k26astro_rt_vehicle_rate_channel_active(NULL) == 0);

    /* ---- NULL-predicate fallback: always-active ---------------- */
    struct K26AstroVehicle veh2 = { 0, 0, 0.0 };
    K26AstroVehicleRateChannel *ch2 =
        k26astro_rt_vehicle_rate_channel_new(w, &veh2, 10.0, NULL,
                                              vehicle_step_, NULL);
    assert(ch2 != NULL);
    assert(k26astro_rt_vehicle_rate_channel_active(ch2) == 1);
    for (int i = 0; i < 10; i++) k26astro_world_step(w, 0.1);
    /* 10 Hz · 1 s = 10 steps. */
    assert(veh2.step_count >= 9 && veh2.step_count <= 11);
    k26astro_rt_vehicle_rate_channel_destroy(ch2);

    k26astro_world_destroy(w);
    printf("test_vehicle_rate_channel: OK\n");
    return 0;
}
