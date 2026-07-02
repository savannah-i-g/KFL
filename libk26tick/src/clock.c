/* clock.c — monotonic clock. CLOCK_MONOTONIC, second baseline lazily
 * established at first call. Replaces k26atlas's nav_now_s. */
#include "k26tick.h"

#include <time.h>

static double s_baseline_s = 0.0;
static int    s_have_baseline = 0;

static double mono_now_raw_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

double k26tick_now_s(void)
{
    if (!s_have_baseline) {
        s_baseline_s = mono_now_raw_s();
        s_have_baseline = 1;
        return 0.0;
    }
    return mono_now_raw_s() - s_baseline_s;
}
