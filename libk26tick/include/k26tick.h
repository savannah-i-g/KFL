/* libk26tick — fixed-step accumulator + named rate channels.
 *
 * A K26TickWorld holds N rate channels. Each channel has a target Hz
 * and a per-step callback. Per render frame the consumer calls
 * k26tick_advance(world, wallclock_dt_s); the world dispatches each
 * channel's callback as many times as needed (with that channel's
 * fixed step_dt), absorbing the leftover into the channel's
 * accumulator.
 *
 * Why an accumulator and not a naive sub-tick split: given the same
 * world seed and the same wallclock dt pattern, an accumulator
 * produces the same step sequence — replay-safe. The cost is one
 * residual double per channel.
 *
 * Spiral-of-death guard: each advance clamps wallclock_dt to
 * 5 × longest registered step_dt. A 500 ms render hitch can't queue
 * 60 phys steps that themselves take >100 ms each.
 *
 * Threading: one K26TickWorld is used by one thread. No locking. */
#ifndef K26TICK_H
#define K26TICK_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define K26TICK_LIB_VERSION "0.1.0"

typedef struct K26TickWorld K26TickWorld;
typedef void (*K26TickStepFn)(double step_dt_s, void *user);

typedef struct {
    uint32_t id;
} K26TickChannel;

#define K26TICK_CHANNEL_NULL ((K26TickChannel){0xFFFFFFFFu})

K26TickWorld *k26tick_open(void);
void          k26tick_close(K26TickWorld *w);

/* Register a fixed-rate channel. hz==0 means render-rate: the callback
 * fires once per advance with the actual wallclock dt (no accumulator).
 * Returns K26TICK_CHANNEL_NULL on out-of-slots (cap 32). */
K26TickChannel k26tick_add_channel(K26TickWorld *w, const char *name,
                                   double hz, K26TickStepFn fn, void *user);

/* Disable/re-enable a channel without removing it. Disabled channels
 * accumulate nothing and don't fire — accumulator resets on re-enable
 * so the channel doesn't burst-catchup. */
void k26tick_set_enabled(K26TickWorld *w, K26TickChannel ch, int on);

/* Live rate change. Updates the channel's target hz and recomputes
 * its step_dt accordingly. The accumulator is preserved (so the
 * channel doesn't lose partial work or burst-catchup), and the
 * world's longest_step_dt_s is recomputed across all enabled
 * channels. hz=0 switches the channel to render-rate (fire once
 * per advance with the wallclock dt). Returns 0 on success,
 * -1 on null world / bad handle. */
int k26tick_channel_set_hz(K26TickWorld *w, K26TickChannel ch, double hz);

/* Read the channel's current hz (0 for render-rate). Returns -1 on
 * null world / bad handle. */
double k26tick_channel_get_hz(const K26TickWorld *w, K26TickChannel ch);

/* Advance the world by wallclock_dt_s. Calls each channel's fn the
 * right number of times with its own step_dt. Returns interpolation
 * alpha [0,1] for the slowest fixed channel — useful for visual blend
 * between two physics states. Returns 0.0 if no fixed channels exist. */
double k26tick_advance(K26TickWorld *w, double wallclock_dt_s);

/* Monotonic-clock helper. Seconds since first call (or some fixed
 * baseline). Replaces nav_now_s. CLOCK_MONOTONIC backed. */
double k26tick_now_s(void);

/* Diagnostic — number of registered channels. */
size_t k26tick_channel_count(const K26TickWorld *w);

#ifdef __cplusplus
}
#endif

#endif /* K26TICK_H */
