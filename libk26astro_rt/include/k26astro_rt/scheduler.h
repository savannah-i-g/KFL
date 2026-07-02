/* libk26astro_rt — multi-rate scheduler API.
 *
 * Wraps libk26tick. Three channels are registered at world create:
 *   orbit   (hz=0, fires once per advance with wallclock dt)
 *   spin    (default 10 Hz, attitude propagation)
 *   render  (default 60 Hz, scene update for libk26astro_render)
 *
 * The orbit channel always runs at "wallclock" rate; the inner
 * substep granularity is controlled by the active integrator
 * (k26astro_grav_advise_step). Spin and render are fixed-rate
 * accumulator channels per libk26tick semantics. */
#ifndef K26ASTRO_RT_SCHEDULER_H
#define K26ASTRO_RT_SCHEDULER_H

#include "k26astro_rt/world.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Default rates set at world create. */
#define K26ASTRO_DEFAULT_SPIN_HZ    10.0
#define K26ASTRO_DEFAULT_RENDER_HZ  60.0

/* Setters. Both clamp to [0.01, 10000.0] Hz to avoid spiral-of-death. */
int    k26astro_scheduler_set_spin_hz   (K26AstroWorld *world, double hz);
int    k26astro_scheduler_set_render_hz (K26AstroWorld *world, double hz);

/* Get current rates. */
double k26astro_scheduler_spin_hz   (const K26AstroWorld *world);
double k26astro_scheduler_render_hz (const K26AstroWorld *world);

/* Interpolation alpha [0,1] for the slowest fixed channel. Useful
 * for visual blend between two physics states. Forwarded from
 * libk26tick. */
double k26astro_scheduler_alpha     (const K26AstroWorld *world);

/* ---- Per-vehicle rate channel --------------------------------- *
 *
 * A sub-orbit-channel cadence pinned to one vehicle. The world's
 * orbit channel runs at wallclock rate; running it faster to absorb
 * RCS-pulse or engine-ignition transients would be wasteful for
 * every other body in the world. A vehicle rate channel lives
 * alongside orbit/spin/render at a caller-chosen Hz and dispatches
 * a per-vehicle step callback that is gated on the vehicle's
 * active-vs-rails state.
 *
 * Lifecycle: each channel allocates one libk26tick slot at create.
 * The vehicle handle is forward-declared (K26AstroVehicle defined
 * in libk26astro_vehicle); the channel never inspects vehicle
 * internals — it only passes the handle to the caller-supplied
 * predicate + step functions.
 *
 * Slot accounting: libk26tick caps total channels per world at 32
 * (orbit + spin + render consume 3; ≤29 vehicle channels per world
 * lifetime). Destroying a vehicle rate channel disables its slot
 * but does not reclaim the slot for new creation. */
struct K26AstroVehicle;  /* defined in libk26astro_vehicle */

typedef struct K26AstroVehicleRateChannel K26AstroVehicleRateChannel;

/* Predicate: returns non-zero when the vehicle should receive its
 * per-channel step (active mode). When the predicate returns 0
 * (rails mode, deactivated mission segment) the channel skips its
 * step and accumulates as if disabled. Passing NULL is equivalent
 * to a predicate that always returns 1. */
typedef int  (*K26AstroVehicleRatePredicateFn)(
                  const struct K26AstroVehicle *v, void *ctx);

/* Per-channel step. Called at the channel's fixed step_dt (1 / hz). */
typedef void (*K26AstroVehicleRateStepFn)     (
                  struct K26AstroVehicle *v, double dt, void *ctx);

/* Create a per-vehicle rate channel at `hz` (clamped to [0.01, 10000]
 * Hz per libk26tick). Returns NULL on null world / vehicle / step,
 * or when libk26tick slots are exhausted. */
K26AstroVehicleRateChannel *
k26astro_rt_vehicle_rate_channel_new(K26AstroWorld                  *world,
                                     struct K26AstroVehicle         *vehicle,
                                     double                          hz,
                                     K26AstroVehicleRatePredicateFn  is_active,
                                     K26AstroVehicleRateStepFn       step,
                                     void                           *ctx);

/* Disable the underlying libk26tick slot and free the wrapper.
 * No-op on NULL. */
void k26astro_rt_vehicle_rate_channel_destroy(
    K26AstroVehicleRateChannel *ch);

/* Live rate change. Returns 0 on success, non-zero on null channel
 * or libk26tick rejection. */
int  k26astro_rt_vehicle_rate_channel_set_hz(
    K26AstroVehicleRateChannel *ch, double hz);

/* Returns the channel's current Hz. 0 if NULL. */
double k26astro_rt_vehicle_rate_channel_hz(
    const K26AstroVehicleRateChannel *ch);

/* Diagnostic: returns 1 if the predicate currently says the vehicle
 * is active, 0 otherwise (or on NULL channel / NULL predicate
 * shortcuts to 1). */
int  k26astro_rt_vehicle_rate_channel_active(
    const K26AstroVehicleRateChannel *ch);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_SCHEDULER_H */
