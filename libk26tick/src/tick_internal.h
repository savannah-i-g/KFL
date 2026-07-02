/* tick_internal.h — private types shared across libk26tick TUs. */
#ifndef K26TICK_INTERNAL_H
#define K26TICK_INTERNAL_H

#include "k26tick.h"

#define K26TICK_CHANNEL_CAP 32

typedef struct {
    int           enabled;
    double        hz;          /* 0 = render-rate */
    double        step_dt_s;   /* 1/hz, or 0 for render-rate */
    double        accum_s;     /* unspent wallclock time */
    K26TickStepFn fn;
    void         *user;
    char          name[24];    /* diagnostic */
} K26TickChan;

struct K26TickWorld {
    K26TickChan ch[K26TICK_CHANNEL_CAP];
    size_t      n_channels;
    /* Longest registered step_dt — used to bound the spiral-of-death
     * clamp in k26tick_advance. Updated on add_channel. */
    double      longest_step_dt_s;
};

#endif /* K26TICK_INTERNAL_H */
