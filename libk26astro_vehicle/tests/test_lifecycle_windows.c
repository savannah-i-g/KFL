/* test_lifecycle_windows.c — active / rails sliding mode-state gate.
 *
 * Exercises the lifecycle resolution rule:
 *   - No windows → active.
 *   - Single window: inside (including endpoints) → window's kind;
 *     outside → active default.
 *   - Overlapping windows: latest registration wins inside the overlap.
 *   - Disjoint windows: each respected.
 *   - Reversed window (t_end < t_start) rejected. */

#include "k26astro_vehicle/lifecycle.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"

#include <assert.h>
#include <stdio.h>

static K26AstroEpoch ep_plus_(K26AstroEpoch base, double delta_s)
{
    k26astro_epoch_add_seconds(&base, delta_s);
    return base;
}

int main(void)
{
    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);

    K26AstroEpoch t0  = k26astro_epoch_j2000_tt();
    K26AstroEpoch t5  = ep_plus_(t0,  5.0);
    K26AstroEpoch t10 = ep_plus_(t0, 10.0);
    K26AstroEpoch t15 = ep_plus_(t0, 15.0);
    K26AstroEpoch t20 = ep_plus_(t0, 20.0);
    K26AstroEpoch t25 = ep_plus_(t0, 25.0);
    K26AstroEpoch t30 = ep_plus_(t0, 30.0);
    K26AstroEpoch t35 = ep_plus_(t0, 35.0);
    K26AstroEpoch t40 = ep_plus_(t0, 40.0);

    /* (a) Default — no windows → active everywhere. */
    assert(k26astro_vehicle_lifecycle_window_count(v) == 0);
    assert(k26astro_vehicle_is_active_at(v, t10) == true);
    assert(k26astro_vehicle_is_active_at(v, t30) == true);

    /* (b) Single rails window [10, 20]. */
    int w0 = k26astro_vehicle_schedule_rails_window(v, t10, t20);
    assert(w0 == 0);
    assert(k26astro_vehicle_lifecycle_window_count(v) == 1);

    assert(k26astro_vehicle_is_active_at(v, t5)  == true);   /* before */
    assert(k26astro_vehicle_is_active_at(v, t10) == false);  /* start inclusive */
    assert(k26astro_vehicle_is_active_at(v, t15) == false);
    assert(k26astro_vehicle_is_active_at(v, t20) == false);  /* end inclusive */
    assert(k26astro_vehicle_is_active_at(v, t25) == true);   /* after */

    /* (c) Overlapping active window [15, 25] — registered later, so
     * latest-wins must give active inside the overlap. */
    int w1 = k26astro_vehicle_schedule_active_window(v, t15, t25);
    assert(w1 == 1);
    assert(k26astro_vehicle_is_active_at(v, t15) == true);   /* active overrides rails */
    assert(k26astro_vehicle_is_active_at(v, t20) == true);
    assert(k26astro_vehicle_is_active_at(v, t25) == true);
    /* t10 only contained by the rails window, so still rails. */
    assert(k26astro_vehicle_is_active_at(v, t10) == false);
    /* t5 outside everything → default active. */
    assert(k26astro_vehicle_is_active_at(v, t5)  == true);

    /* (d) Disjoint third rails window [30, 40]. */
    int w2 = k26astro_vehicle_schedule_rails_window(v, t30, t40);
    assert(w2 == 2);
    assert(k26astro_vehicle_is_active_at(v, t30) == false);
    assert(k26astro_vehicle_is_active_at(v, t35) == false);
    assert(k26astro_vehicle_is_active_at(v, t40) == false);
    /* Mid-overlap zone untouched. */
    assert(k26astro_vehicle_is_active_at(v, t25) == true);

    /* (e) Reversed window rejected. */
    assert(k26astro_vehicle_schedule_active_window(v, t40, t30) == -1);
    assert(k26astro_vehicle_schedule_rails_window (v, t40, t30) == -1);
    assert(k26astro_vehicle_lifecycle_window_count(v) == 3);

    /* (f) Zero-duration window [t, t] honoured. */
    int w3 = k26astro_vehicle_schedule_active_window(v, t10, t10);
    assert(w3 == 3);
    assert(k26astro_vehicle_is_active_at(v, t10) == true);   /* latest-wins */

    /* NULL safety: defaults to active. */
    assert(k26astro_vehicle_is_active_at(NULL, t0) == true);
    assert(k26astro_vehicle_lifecycle_window_count(NULL) == 0);
    assert(k26astro_vehicle_schedule_active_window(NULL, t0, t10) == -1);
    assert(k26astro_vehicle_schedule_rails_window (NULL, t0, t10) == -1);

    k26astro_vehicle_destroy(v);

    printf("test_lifecycle_windows: OK\n");
    return 0;
}
