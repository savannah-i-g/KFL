/* test_channel_rate.c - live channel rate change gate.
 *
 * Verifies k26tick_channel_set_hz: updates a channel's target hz
 * mid-run, recomputes step_dt, and preserves the accumulator so
 * the channel doesn't lose partial work or burst-catchup.
 *
 * Two coverage axes:
 *   Test 1: Cadence change. Register a 60 Hz channel. Advance the
 *     world 1 second; expect ~60 callbacks. Change the rate to
 *     30 Hz. Advance another second; expect ~30 callbacks.
 *   Test 2: Get-after-set round-trip. set_hz(120) followed by
 *     get_hz() returns 120. set_hz(0) for render-rate is permitted
 *     and get_hz returns 0. */
#include "k26tick.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int g_tick_count = 0;

static void counter_cb_(double step_dt_s, void *user)
{
    (void)step_dt_s;
    (void)user;
    g_tick_count++;
}

static int test_cadence_change(void)
{
    K26TickWorld *w = k26tick_open();
    assert(w);
    K26TickChannel ch = k26tick_add_channel(w, "phys", 60.0,
                                              counter_cb_, NULL);
    assert(ch.id != 0xFFFFFFFFu);

    g_tick_count = 0;
    /* Drive 1 second of wall time in 10ms increments — small enough
     * that 60 Hz fires at most 1× per advance, avoiding burst-catchup. */
    for (int i = 0; i < 100; i++) k26tick_advance(w, 0.01);
    int at_60hz = g_tick_count;
    fprintf(stderr, "test_cadence_change: 1s @ 60Hz → %d ticks "
                    "(expected ~60)\n", at_60hz);
    if (at_60hz < 58 || at_60hz > 62) {
        fprintf(stderr, "FAIL: 60Hz tick count %d out of [58, 62]\n",
                at_60hz);
        k26tick_close(w);
        return 1;
    }

    /* Live rate change → 30 Hz. */
    int rc = k26tick_channel_set_hz(w, ch, 30.0);
    if (rc != 0) {
        fprintf(stderr, "FAIL: set_hz returned %d\n", rc);
        k26tick_close(w);
        return 1;
    }
    if (k26tick_channel_get_hz(w, ch) != 30.0) {
        fprintf(stderr, "FAIL: get_hz returned %.3f after set_hz(30)\n",
                k26tick_channel_get_hz(w, ch));
        k26tick_close(w);
        return 1;
    }

    g_tick_count = 0;
    for (int i = 0; i < 100; i++) k26tick_advance(w, 0.01);
    int at_30hz = g_tick_count;
    fprintf(stderr, "test_cadence_change: 1s @ 30Hz → %d ticks "
                    "(expected ~30)\n", at_30hz);
    if (at_30hz < 28 || at_30hz > 32) {
        fprintf(stderr, "FAIL: 30Hz tick count %d out of [28, 32]\n",
                at_30hz);
        k26tick_close(w);
        return 1;
    }

    k26tick_close(w);
    return 0;
}

static int test_get_after_set(void)
{
    K26TickWorld *w = k26tick_open();
    assert(w);
    K26TickChannel ch = k26tick_add_channel(w, "render", 60.0,
                                              counter_cb_, NULL);
    assert(ch.id != 0xFFFFFFFFu);

    /* Initial rate. */
    if (k26tick_channel_get_hz(w, ch) != 60.0) {
        fprintf(stderr, "FAIL: initial hz != 60\n");
        k26tick_close(w);
        return 1;
    }

    /* Set + get round-trip. */
    assert(k26tick_channel_set_hz(w, ch, 120.0) == 0);
    if (k26tick_channel_get_hz(w, ch) != 120.0) {
        fprintf(stderr, "FAIL: get_hz != 120 after set\n");
        k26tick_close(w);
        return 1;
    }

    /* Render-rate (hz=0) is legal — channel becomes "fire once per
     * advance with wallclock dt". */
    assert(k26tick_channel_set_hz(w, ch, 0.0) == 0);
    if (k26tick_channel_get_hz(w, ch) != 0.0) {
        fprintf(stderr, "FAIL: get_hz != 0 after set_hz(0)\n");
        k26tick_close(w);
        return 1;
    }

    /* Negative hz is rejected. */
    int rc = k26tick_channel_set_hz(w, ch, -5.0);
    if (rc == 0) {
        fprintf(stderr, "FAIL: set_hz(-5) returned success\n");
        k26tick_close(w);
        return 1;
    }

    /* Null world returns -1. */
    if (k26tick_channel_set_hz(NULL, ch, 60.0) != -1) {
        fprintf(stderr, "FAIL: NULL world didn't return -1\n");
        k26tick_close(w);
        return 1;
    }

    k26tick_close(w);
    return 0;
}

int main(void)
{
    if (test_cadence_change()) return 1;
    if (test_get_after_set())  return 1;
    fprintf(stderr, "test_channel_rate: OK\n");
    return 0;
}
