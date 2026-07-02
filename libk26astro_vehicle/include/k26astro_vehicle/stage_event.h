/* k26astro_vehicle/stage_event.h — discrete-event timeline with
 * registration into libk26astro_grav's event-time root-finder.
 *
 * A stage event is an (epoch, mass_drop_kg, post-event diagonal
 * inertia) record. Scheduled events compose one-to-one with the
 * gravity integrator's event registry: each event becomes a
 * K26AstroGravEvent with
 *
 *   predicate(state, t_s, ctx)  : returns 1 when the current
 *                                  simulation time (seconds since
 *                                  J2000) is at or past the
 *                                  scheduled epoch, 0 otherwise.
 *                                  Monotonic across any candidate
 *                                  integrator substep.
 *
 *   handler  (state, t_s, ctx)  : runs after the bisector has
 *                                  advanced to t_event. Mutates the
 *                                  vehicle's basic mass, propagates
 *                                  the new mass to the bound body
 *                                  via k26astro_body_set_mass, and
 *                                  updates the diagonal entries of
 *                                  the Ext inertia tensor (off-
 *                                  diagonals preserved).
 *
 * The vehicle owns the timeline; the gravity state owns the
 * registry. Stage events are flushed to a specific gravity state via
 * k26astro_vehicle_register_stage_events_with — scheduling is
 * decoupled from gravity so per-vehicle unit tests can run without
 * a gravity state.
 *
 * Determinism: predicates are monotone in simulation time;
 * bisection is bit-stable on IEEE-754 (per libk26astro_grav's
 * advance_with_events.c). Multiple events at the same epoch fire in
 * registration order. */
#ifndef K26ASTRO_VEHICLE_STAGE_EVENT_H
#define K26ASTRO_VEHICLE_STAGE_EVENT_H

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"
#include "k26astro_grav/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Record a stage event in the vehicle's timeline. The timeline
 * stays sorted ascending by epoch on every insertion. Returns the
 * resulting sorted index (>= 0) on success, -1 on NULL vehicle or
 * allocation failure.
 *
 * mass_drop_kg is subtracted from basic_mass when the event fires
 * (clamped to >= 0). The post-event diagonal triple replaces the
 * diagonal entries of the Ext inertia tensor; off-diagonals are
 * preserved (a staged-substructure drop typically symmetrises the
 * remaining vehicle, and off-diagonal mutation across stage events
 * is a tank-distribution concern that this surface does not
 * currently model). */
int k26astro_vehicle_schedule_stage_event(K26AstroVehicle *v,
                                          K26AstroEpoch epoch,
                                          double mass_drop_kg,
                                          double ixx_after,
                                          double iyy_after,
                                          double izz_after);

/* Total number of scheduled events (including events that have
 * already fired). Returns 0 on NULL. */
int           k26astro_vehicle_event_count   (const K26AstroVehicle *v);

/* Returns the epoch of the event at sorted index `idx`. Returns the
 * J2000 TT epoch on out-of-range index or NULL vehicle. */
K26AstroEpoch k26astro_vehicle_event_epoch_at(const K26AstroVehicle *v, int idx);

/* Flush every scheduled stage event into `grav`'s K26AstroGravEvent
 * registry. Each registered event carries a vehicle-owned context
 * record so the handler can mutate vehicle state at fire time.
 *
 * Precondition: the vehicle must outlive `grav`. Context records
 * are vehicle-owned and freed by k26astro_vehicle_destroy; if
 * gravity-state callbacks fire after the vehicle has been
 * destroyed, they touch freed memory.
 *
 * Calling this twice with the same (vehicle, grav) pair will
 * double-register the events — this is a caller error and is not
 * currently guarded against. The expected sequence is: configure
 * stages, register once, integrate.
 *
 * Returns K26ASTRO_E_OK on success, K26ASTRO_E_NULL on NULL inputs,
 * K26ASTRO_E_ALLOC if the gravity registry cannot grow. */
int k26astro_vehicle_register_stage_events_with(K26AstroVehicle *v,
                                                K26AstroGravState *grav);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_STAGE_EVENT_H */
