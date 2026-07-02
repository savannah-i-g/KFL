/* lifecycle.c — active / rails sliding mode-state machine.
 *
 * Windows are stored in registration order. is_active_at scans
 * back-to-front so the latest-registered window containing t wins
 * on overlap. Outside all windows the default is active.
 *
 * Complexity: O(n_windows) per is_active_at call. n_windows is
 * typically small (single-digit per vehicle); a sorted-by-end +
 * binary-search optimisation is worth adding only when mission-
 * scale schedules surface. */

#include "k26astro_vehicle/lifecycle.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"
#include "vehicle_internal.h"

#include <stdlib.h>

static int append_window_(K26AstroVehicle *v,
                          K26AstroEpoch t_start, K26AstroEpoch t_end,
                          K26AstroLifecycleKind_ kind)
{
    if (!v) return -1;
    /* Reject reversed windows (t_end < t_start). */
    if (k26astro_epoch_diff_seconds(&t_end, &t_start) < 0.0) return -1;
    if (v->n_windows >= v->cap_windows) {
        int new_cap = v->cap_windows ? (v->cap_windows * 2) : 4;
        K26AstroVehicleLifecycleWindow_ *grown = realloc(
            v->windows,
            (size_t)new_cap * sizeof(K26AstroVehicleLifecycleWindow_));
        if (!grown) return -1;
        v->windows     = grown;
        v->cap_windows = new_cap;
    }
    v->windows[v->n_windows].t_start = t_start;
    v->windows[v->n_windows].t_end   = t_end;
    v->windows[v->n_windows].kind    = kind;
    return v->n_windows++;
}

int k26astro_vehicle_schedule_active_window(K26AstroVehicle *v,
                                            K26AstroEpoch t_start,
                                            K26AstroEpoch t_end)
{
    return append_window_(v, t_start, t_end, K26ASTRO_LC_ACTIVE);
}

int k26astro_vehicle_schedule_rails_window(K26AstroVehicle *v,
                                           K26AstroEpoch t_start,
                                           K26AstroEpoch t_end)
{
    return append_window_(v, t_start, t_end, K26ASTRO_LC_RAILS);
}

bool k26astro_vehicle_is_active_at(const K26AstroVehicle *v, K26AstroEpoch t)
{
    if (!v) return true;
    if (v->n_windows == 0) return true;
    /* Latest-registered wins on overlap: walk back-to-front. */
    for (int i = v->n_windows - 1; i >= 0; i--) {
        double dt_lo = k26astro_epoch_diff_seconds(&t, &v->windows[i].t_start);
        double dt_hi = k26astro_epoch_diff_seconds(&v->windows[i].t_end, &t);
        if (dt_lo >= 0.0 && dt_hi >= 0.0) {
            return v->windows[i].kind == K26ASTRO_LC_ACTIVE;
        }
    }
    /* t outside every window — default active. */
    return true;
}

int k26astro_vehicle_lifecycle_window_count(const K26AstroVehicle *v)
{
    return v ? v->n_windows : 0;
}
