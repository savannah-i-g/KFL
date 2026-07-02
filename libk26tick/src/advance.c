/* advance.c — accumulator dispatch + spiral-of-death clamp.
 *
 * Per channel:
 *   - hz > 0: classic accumulator. accum += dt; while accum >= step:
 *     fire(step); accum -= step;
 *   - hz == 0: fire(actual_dt) once. No accumulator.
 *
 * The returned alpha is the interpolation fraction left in the
 * SLOWEST fixed channel's accumulator. A consumer doing visual blend
 * between two physics states should interpolate by (1 - alpha). */
#include "k26tick.h"
#include "tick_internal.h"

double k26tick_advance(K26TickWorld *w, double wallclock_dt_s)
{
    if (!w) return 0.0;
    if (wallclock_dt_s < 0.0) wallclock_dt_s = 0.0;

    /* Spiral-of-death clamp. 5× the longest registered step caps how
     * many catch-up iterations any single advance can queue. With a
     * 1/120s longest step that's 41.7 ms — beyond that we drop the
     * excess wall time (visible as a one-off slow-mo on hitches, never
     * a runaway). */
    double max_dt = w->longest_step_dt_s * 5.0;
    if (max_dt > 0.0 && wallclock_dt_s > max_dt) {
        wallclock_dt_s = max_dt;
    }

    /* Hard outer cap for the render-rate (hz==0) case where no
     * fixed channel exists to bound longest_step_dt_s. */
    if (wallclock_dt_s > 0.5) wallclock_dt_s = 0.5;

    double slowest_step = 0.0;
    double slowest_accum_after = 0.0;

    for (size_t i = 0; i < w->n_channels; i++) {
        K26TickChan *c = &w->ch[i];
        if (!c->enabled) continue;

        if (c->hz <= 0.0) {
            /* Render-rate channel — one call per advance. */
            c->fn(wallclock_dt_s, c->user);
            continue;
        }

        c->accum_s += wallclock_dt_s;
        while (c->accum_s >= c->step_dt_s) {
            c->fn(c->step_dt_s, c->user);
            c->accum_s -= c->step_dt_s;
        }
        if (c->step_dt_s > slowest_step) {
            slowest_step = c->step_dt_s;
            slowest_accum_after = c->accum_s;
        }
    }

    if (slowest_step > 0.0) {
        double alpha = slowest_accum_after / slowest_step;
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;
        return alpha;
    }
    return 0.0;
}
