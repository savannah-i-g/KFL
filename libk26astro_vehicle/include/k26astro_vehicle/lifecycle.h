/* k26astro_vehicle/lifecycle.h — active / rails windowed scheduler.
 *
 * A vehicle is in one of two lifecycle modes at any given epoch:
 *
 *   active  — the integrator pays per-substep dot_m callback cost
 *             for every engine cluster registered on the vehicle.
 *             Default mode when no windows are scheduled.
 *
 *   rails   — the integrator treats vehicle mass as constant and
 *             skips dot_m callbacks. Used during long quiescent
 *             cruise legs where engines are off (interplanetary
 *             cruise, Hohmann coast).
 *
 * Multiple windows of either kind may be scheduled and may overlap.
 * Resolution rule for the mode at epoch t:
 *
 *   1. No windows scheduled  → active. Small programs work without
 *      any scheduling.
 *   2. At least one window contains t (inclusive endpoints) → the
 *      latest-registered such window decides. (Sliding mode-state
 *      machine; higher-level mission-timeline surfaces compose on
 *      top.)
 *   3. t outside every window  → active. Gaps default to active;
 *      a quiescent gap must be explicitly scheduled as rails to
 *      opt out.
 *
 * Boundary convention: both endpoints inclusive. A zero-duration
 * window [t, t] containing t is honoured.
 *
 * Time-scale handling: k26astro_epoch_diff_seconds asserts equal
 * scales in debug builds. The caller is responsible for using a
 * consistent scale across schedule and is_active_at calls. */
#ifndef K26ASTRO_VEHICLE_LIFECYCLE_H
#define K26ASTRO_VEHICLE_LIFECYCLE_H

#include <stdbool.h>

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Schedule an active window. Returns the window index (>= 0) on
 * success, -1 on NULL vehicle, allocation failure, or reversed
 * window (t_end < t_start). The vehicle keeps windows in
 * registration order; latest registration wins on overlap. */
int k26astro_vehicle_schedule_active_window(K26AstroVehicle *v,
                                            K26AstroEpoch t_start,
                                            K26AstroEpoch t_end);

/* Schedule a rails window. Same semantics + return as
 * schedule_active_window; the kind field distinguishes. */
int k26astro_vehicle_schedule_rails_window (K26AstroVehicle *v,
                                            K26AstroEpoch t_start,
                                            K26AstroEpoch t_end);

/* Returns true if the vehicle is in active mode at epoch t per the
 * rule above; false otherwise. Returns true on NULL vehicle
 * (defensive: a caller without a vehicle handle is by definition
 * in the no-windows default). */
bool k26astro_vehicle_is_active_at(const K26AstroVehicle *v,
                                   K26AstroEpoch t);

/* Total number of scheduled windows of both kinds. Diagnostic. */
int  k26astro_vehicle_lifecycle_window_count(const K26AstroVehicle *v);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_LIFECYCLE_H */
