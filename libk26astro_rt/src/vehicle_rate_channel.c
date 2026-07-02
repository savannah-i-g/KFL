/* vehicle_rate_channel.c — per-vehicle sub-orbit-channel cadence.
 *
 * Wraps libk26tick. One libk26tick channel per vehicle rate channel;
 * the channel's user-pointer carries a wrapper that holds the
 * vehicle handle + the predicate + step + ctx. The internal callback
 * dispatches predicate-then-step on each fixed-rate firing.
 *
 * Slot conservation: libk26tick caps channels per world at 32 and
 * has no remove primitive. Destroying a vehicle rate channel
 * disables the slot via k26tick_set_enabled(0) but does not reclaim
 * it. Practical limit per world lifetime: ~29 distinct vehicle rate
 * channels (orbit + spin + render claim three). */
#include "k26astro_rt/scheduler.h"
#include "world_internal.h"

#include <stdlib.h>

struct K26AstroVehicleRateChannel {
    K26AstroWorld                  *world;
    struct K26AstroVehicle         *vehicle;
    K26TickChannel                  tick_ch;
    K26AstroVehicleRatePredicateFn  is_active;
    K26AstroVehicleRateStepFn       step;
    void                           *ctx;
    double                          hz;
};

static int channel_handle_valid_(K26TickChannel ch)
{
    return ch.id != K26TICK_CHANNEL_NULL.id;
}

static void vehicle_rate_dispatch_(double step_dt_s, void *user)
{
    K26AstroVehicleRateChannel *c = (K26AstroVehicleRateChannel *)user;
    if (!c || !c->step) return;
    if (c->is_active && !c->is_active(c->vehicle, c->ctx)) return;
    c->step(c->vehicle, step_dt_s, c->ctx);
}

K26AstroVehicleRateChannel *
k26astro_rt_vehicle_rate_channel_new(K26AstroWorld                  *world,
                                     struct K26AstroVehicle         *vehicle,
                                     double                          hz,
                                     K26AstroVehicleRatePredicateFn  is_active,
                                     K26AstroVehicleRateStepFn       step,
                                     void                           *ctx)
{
    if (!world || !world->tick || !vehicle || !step) return NULL;

    K26AstroVehicleRateChannel *c =
        (K26AstroVehicleRateChannel *)calloc(1, sizeof(*c));
    if (!c) return NULL;

    c->world     = world;
    c->vehicle   = vehicle;
    c->is_active = is_active;
    c->step      = step;
    c->ctx       = ctx;
    c->hz        = hz;

    c->tick_ch = k26tick_add_channel(world->tick, "vehicle_rate",
                                      hz, vehicle_rate_dispatch_, c);
    if (!channel_handle_valid_(c->tick_ch)) {
        free(c);
        return NULL;
    }
    /* Read the channel's clamped Hz back so the wrapper reports the
     * effective rate, not the requested one. */
    c->hz = k26tick_channel_get_hz(world->tick, c->tick_ch);
    return c;
}

void k26astro_rt_vehicle_rate_channel_destroy(K26AstroVehicleRateChannel *ch)
{
    if (!ch) return;
    if (ch->world && ch->world->tick && channel_handle_valid_(ch->tick_ch)) {
        k26tick_set_enabled(ch->world->tick, ch->tick_ch, 0);
    }
    free(ch);
}

int k26astro_rt_vehicle_rate_channel_set_hz(K26AstroVehicleRateChannel *ch,
                                            double hz)
{
    if (!ch || !ch->world || !ch->world->tick) return -1;
    int rc = k26tick_channel_set_hz(ch->world->tick, ch->tick_ch, hz);
    if (rc != 0) return rc;
    ch->hz = k26tick_channel_get_hz(ch->world->tick, ch->tick_ch);
    return 0;
}

double k26astro_rt_vehicle_rate_channel_hz(const K26AstroVehicleRateChannel *ch)
{
    return ch ? ch->hz : 0.0;
}

int k26astro_rt_vehicle_rate_channel_active(const K26AstroVehicleRateChannel *ch)
{
    if (!ch) return 0;
    if (!ch->is_active) return 1;
    return ch->is_active(ch->vehicle, ch->ctx) ? 1 : 0;
}
