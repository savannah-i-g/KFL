/* world.c — K26TickWorld lifecycle and channel registry.
 *
 * Channels live in a fixed-size array (cap 32) so registration is
 * branch-free and we never allocate after k26tick_open. The cap is a
 * compile-time choice — bump K26TICK_CHANNEL_CAP if a consumer wants
 * more. k26atlas uses ~6 channels. */
#include "k26tick.h"
#include "tick_internal.h"

#include <stdlib.h>
#include <string.h>

K26TickWorld *k26tick_open(void)
{
    K26TickWorld *w = (K26TickWorld *)calloc(1, sizeof(*w));
    if (!w) return NULL;
    w->n_channels = 0;
    w->longest_step_dt_s = 0.0;
    return w;
}

void k26tick_close(K26TickWorld *w)
{
    if (!w) return;
    free(w);
}

K26TickChannel k26tick_add_channel(K26TickWorld *w, const char *name,
                                   double hz, K26TickStepFn fn, void *user)
{
    if (!w || !fn || w->n_channels >= K26TICK_CHANNEL_CAP) {
        return K26TICK_CHANNEL_NULL;
    }
    K26TickChan *c = &w->ch[w->n_channels];
    c->enabled = 1;
    c->hz = hz;
    c->step_dt_s = (hz > 0.0) ? (1.0 / hz) : 0.0;
    c->accum_s = 0.0;
    c->fn = fn;
    c->user = user;
    /* Copy name for diagnostics — bounded, NUL-terminated. */
    if (name) {
        size_t n = strlen(name);
        if (n >= sizeof(c->name)) n = sizeof(c->name) - 1;
        memcpy(c->name, name, n);
        c->name[n] = '\0';
    } else {
        c->name[0] = '\0';
    }
    if (c->step_dt_s > w->longest_step_dt_s) {
        w->longest_step_dt_s = c->step_dt_s;
    }
    K26TickChannel h = { (uint32_t)w->n_channels };
    w->n_channels++;
    return h;
}

void k26tick_set_enabled(K26TickWorld *w, K26TickChannel ch, int on)
{
    if (!w || ch.id >= w->n_channels) return;
    K26TickChan *c = &w->ch[ch.id];
    int was = c->enabled;
    c->enabled = on ? 1 : 0;
    /* Re-enable clears the accumulator so the channel doesn't burst
     * out a backlog of steps the moment it switches on. */
    if (!was && c->enabled) c->accum_s = 0.0;
}

size_t k26tick_channel_count(const K26TickWorld *w)
{
    return w ? w->n_channels : 0;
}

/* Recompute the world-wide longest step_dt by scanning the full
 * channel array. Called after a live rate change since the old
 * longest might have belonged to the changed channel. */
static void recompute_longest_(K26TickWorld *w)
{
    double longest = 0.0;
    for (size_t i = 0; i < w->n_channels; i++) {
        double sd = w->ch[i].step_dt_s;
        if (sd > longest) longest = sd;
    }
    w->longest_step_dt_s = longest;
}

int k26tick_channel_set_hz(K26TickWorld *w, K26TickChannel ch, double hz)
{
    if (!w || ch.id >= w->n_channels) return -1;
    /* Reject NaN, negative, or > 1 MHz. hz=0 is permitted (means
     * render-rate: fire once per advance with wallclock dt). The
     * upper bound mirrors the spiral-of-death guard's working range
     * — calls outside [0, 1e6] are programming errors, not valid
     * input. */
    if (!(hz >= 0.0) || hz > 1.0e6) return -1;
    K26TickChan *c = &w->ch[ch.id];
    c->hz        = hz;
    c->step_dt_s = (hz > 0.0) ? (1.0 / hz) : 0.0;
    /* Accumulator preserved: any partial-tick remainder under the
     * old rate continues to count toward the next tick at the new
     * rate. If the accumulator now exceeds the new step_dt,
     * k26tick_advance's next call will fire the catch-up ticks
     * naturally (subject to the spiral-of-death clamp). */
    recompute_longest_(w);
    return 0;
}

double k26tick_channel_get_hz(const K26TickWorld *w, K26TickChannel ch)
{
    if (!w || ch.id >= w->n_channels) return -1.0;
    return w->ch[ch.id].hz;
}
