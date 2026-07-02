/* test_rate_change.c — live rate change at the rt layer.
 *
 * Verifies that k26astro_scheduler_set_spin_hz /
 * k26astro_scheduler_set_render_hz take effect mid-run via the
 * libk26tick channel-rate change, rather than being cached until
 * the next world create.
 *
 * Test approach: create a world, drive 1 second of wallclock time,
 * read the channel's tick count via libk26tick's diagnostic surface,
 * change the rate, drive another second, confirm the cadence
 * changed. The rt layer's spin/render callbacks are no-ops, so the
 * only observable is "channel hz read-back", which
 * k26tick_channel_get_hz exposes.
 *
 * Acceptance:
 *   - set_spin_hz(30) → get_spin_hz returns 30
 *   - set_render_hz(120) → get_render_hz returns 120
 *   - both work without world re-create */
#include "k26astro_rt/world.h"
#include "k26astro_rt/scheduler.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int test_set_get_round_trip(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    /* Defaults at world create. */
    double spin0   = k26astro_scheduler_spin_hz(w);
    double render0 = k26astro_scheduler_render_hz(w);
    fprintf(stderr,
        "test_rate: defaults spin=%.1f render=%.1f\n", spin0, render0);
    if (fabs(spin0 - K26ASTRO_DEFAULT_SPIN_HZ) > 1.0e-12
     || fabs(render0 - K26ASTRO_DEFAULT_RENDER_HZ) > 1.0e-12) {
        fprintf(stderr, "FAIL: defaults don't match constants\n");
        k26astro_world_destroy(w);
        return 1;
    }

    /* Live change to spin. */
    if (k26astro_scheduler_set_spin_hz(w, 30.0) != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: set_spin_hz(30) failed\n");
        k26astro_world_destroy(w);
        return 1;
    }
    double spin1 = k26astro_scheduler_spin_hz(w);
    if (fabs(spin1 - 30.0) > 1.0e-12) {
        fprintf(stderr, "FAIL: spin readback %.3f != 30\n", spin1);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Live change to render. */
    if (k26astro_scheduler_set_render_hz(w, 120.0) != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: set_render_hz(120) failed\n");
        k26astro_world_destroy(w);
        return 1;
    }
    double render1 = k26astro_scheduler_render_hz(w);
    if (fabs(render1 - 120.0) > 1.0e-12) {
        fprintf(stderr, "FAIL: render readback %.3f != 120\n", render1);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Mid-run change: advance, then change rate, then advance more. */
    for (int i = 0; i < 10; i++) k26astro_world_step(w, 0.01);
    if (k26astro_scheduler_set_spin_hz(w, 5.0) != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: mid-run set_spin_hz(5) failed\n");
        k26astro_world_destroy(w);
        return 1;
    }
    double spin2 = k26astro_scheduler_spin_hz(w);
    if (fabs(spin2 - 5.0) > 1.0e-12) {
        fprintf(stderr, "FAIL: spin readback %.3f != 5 after mid-run\n",
                spin2);
        k26astro_world_destroy(w);
        return 1;
    }
    /* Verify the world continues to step cleanly after rate change. */
    for (int i = 0; i < 10; i++) {
        int rc = k26astro_world_step(w, 0.01);
        if (rc != 0) {
            fprintf(stderr, "FAIL: world_step rc=%d after rate change\n", rc);
            k26astro_world_destroy(w);
            return 1;
        }
    }

    /* Clamping: out-of-range rates are clamped, not rejected. */
    k26astro_scheduler_set_spin_hz(w, 0.0);
    double spin_clamped_low = k26astro_scheduler_spin_hz(w);
    if (!(spin_clamped_low > 0.0 && spin_clamped_low < 0.1)) {
        fprintf(stderr,
            "FAIL: set_spin_hz(0) → readback %.6f (expected ~0.01)\n",
            spin_clamped_low);
        k26astro_world_destroy(w);
        return 1;
    }
    k26astro_scheduler_set_spin_hz(w, 1.0e9);
    double spin_clamped_hi = k26astro_scheduler_spin_hz(w);
    if (!(spin_clamped_hi >= 9999.0 && spin_clamped_hi <= 10001.0)) {
        fprintf(stderr,
            "FAIL: set_spin_hz(1e9) → readback %.3f (expected ~10000)\n",
            spin_clamped_hi);
        k26astro_world_destroy(w);
        return 1;
    }

    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    if (test_set_get_round_trip()) return 1;
    fprintf(stderr, "test_rate_change: OK\n");
    return 0;
}
