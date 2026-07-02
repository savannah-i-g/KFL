/* scheduler.c — libk26tick adapter.
 *
 * Registers three named channels on the world's tick context:
 *   orbit   (hz=0)   — fires once per advance with the wallclock dt
 *   spin    (10 Hz)  — attitude propagation (libk26astro_body)
 *   render  (60 Hz)  — scene update hook (libk26astro_render reads it)
 *
 * Rates are clamped to [0.01, 10000] Hz to avoid spiral-of-death even
 * if a caller passes nonsense. */
#include "k26astro_rt/scheduler.h"
#include "world_internal.h"

static double clamp_hz_(double hz)
{
    if (!(hz >= 0.01))    return 0.01;
    if (!(hz <= 10000.0)) return 10000.0;
    return hz;
}

int k26astro_rt_scheduler_init(K26AstroWorld *world)
{
    world->tick = k26tick_open();
    if (!world->tick) return -1;

    world->spin_hz   = K26ASTRO_DEFAULT_SPIN_HZ;
    world->render_hz = K26ASTRO_DEFAULT_RENDER_HZ;

    world->chan_orbit  = k26tick_add_channel(world->tick, "orbit", 0.0,
                                              k26astro_rt_orbit_step_cb, world);
    world->chan_spin   = k26tick_add_channel(world->tick, "spin",
                                              world->spin_hz,
                                              k26astro_rt_spin_step_cb, world);
    world->chan_render = k26tick_add_channel(world->tick, "render",
                                              world->render_hz,
                                              k26astro_rt_render_step_cb, world);
    return 0;
}

void k26astro_rt_scheduler_destroy(K26AstroWorld *world)
{
    if (!world || !world->tick) return;
    k26tick_close(world->tick);
    world->tick = NULL;
}

/* Live rate change via libk26tick's set_hz. The world's
 * channel handles are stashed at scheduler_init time, so we can
 * push the new rate into the channel in-place without re-
 * registering. libk26tick preserves the accumulator across the
 * change, so partial-tick work isn't lost and the channel doesn't
 * burst-catchup. */
int k26astro_scheduler_set_spin_hz(K26AstroWorld *world, double hz)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    double clamped = clamp_hz_(hz);
    world->spin_hz = clamped;
    if (world->tick) {
        if (k26tick_channel_set_hz(world->tick, world->chan_spin,
                                     clamped) != 0) {
            return -K26ASTRO_RT_E_BAD_ARG;
        }
    }
    return K26ASTRO_RT_OK;
}

int k26astro_scheduler_set_render_hz(K26AstroWorld *world, double hz)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    double clamped = clamp_hz_(hz);
    world->render_hz = clamped;
    if (world->tick) {
        if (k26tick_channel_set_hz(world->tick, world->chan_render,
                                     clamped) != 0) {
            return -K26ASTRO_RT_E_BAD_ARG;
        }
    }
    return K26ASTRO_RT_OK;
}

double k26astro_scheduler_spin_hz(const K26AstroWorld *world)
{
    return world ? world->spin_hz : 0.0;
}

double k26astro_scheduler_render_hz(const K26AstroWorld *world)
{
    return world ? world->render_hz : 0.0;
}

double k26astro_scheduler_alpha(const K26AstroWorld *world)
{
    /* libk26tick exposes alpha via the return value of k26tick_advance,
     * not as a query. In v0.1 we cache the last advance return; for now
     * return 0 as a stable placeholder. v0.2 housekeeping: cache it. */
    (void)world;
    return 0.0;
}

/* Spin / render callbacks. v0.1 are no-ops — attitude propagation
 * lives in libk26astro_body (per-step within orbit_step), render
 * hook lives in libk26astro_render (via the plot3d extra-nodes cb,
 * not the scheduler). These exist so the libk26tick channels have
 * something to dispatch. */
void k26astro_rt_spin_step_cb(double dt_s, void *user)
{
    (void)dt_s; (void)user;
}

void k26astro_rt_render_step_cb(double dt_s, void *user)
{
    (void)dt_s; (void)user;
}
