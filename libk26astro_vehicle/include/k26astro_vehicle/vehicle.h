/* k26astro_vehicle/vehicle.h — spacecraft identity + composition root.
 *
 * A K26AstroVehicle is the identity of a spacecraft: dry mass plus
 * mass growth allowance, full 3x3 inertia tensor in the body frame,
 * centre-of-mass offset, optional bind to a K26AstroBody (the
 * position / velocity / attitude provider), composition slots for
 * subsystems (engine clusters, tanks, thermal networks, RCS quads,
 * EPS buses, attitude controllers, sensors), a discrete stage-event
 * timeline, and active/rails lifecycle windows.
 *
 * Subsystem libraries hang their per-vehicle state off the vehicle
 * handle. The vehicle owns typed opaque handles to subsystem state;
 * the subsystem libraries own the heap allocations themselves.
 * Subsystem-handle lifetime is mediated by a back-reference protocol
 * (k26astro_vehicle_invalidate_slot for subsystem-initiated
 * invalidation; a destroy-time notify shim — overridable per
 * subsystem via weak linkage — fires when the vehicle goes away).
 *
 * Mass-budget terms follow NASA-STD-1000 / GSFC-STD-1000RevH /
 * AIAA S-120A-2015 discipline:
 *
 *   basic_mass_kg    — current actual mass (NASA basic mass)
 *   mga_kg           — mass growth allowance (AIAA S-120A-2015)
 *   predicted_mass   — basic + mga
 *   mev_mass         — propellant-sizing mass (basic + mga + reserve)
 *
 * Full 3x3 inertia uses K26AstroAttitudeStateExt from libk26astro_body.
 * The diagonal setter zeroes off-diagonals and routes through the
 * same Ext path; vehicles never use the K26V3 inertia_diag path that
 * celestial bodies use.
 *
 * Stage events compile to K26AstroGravEvent registrations against a
 * K26AstroGravState (see stage_event.h). Events evolve mass and the
 * inertia tensor at integration-quality timing via the event-time
 * root-finder inside the gravity integrator.
 *
 * Active / rails lifecycle windows (lifecycle.h) gate per-substep
 * dot_m callback cost during quiescent cruise phases. Default mode
 * when no windows are scheduled is active.
 *
 * Determinism: single-threaded. The mass_accum field is a scalar
 * double, not atomic. A future parallel engine-cluster dispatch
 * would need atomic accumulation. */
#ifndef K26ASTRO_VEHICLE_VEHICLE_H
#define K26ASTRO_VEHICLE_VEHICLE_H

#include <stdint.h>
#include <stddef.h>

#include "k26astro_core/epoch.h"
#include "k26astro_body/body.h"
#include "k26astro_body/attitude.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Opaque vehicle handle ------------------------------------- */

/* The vehicle struct definition lives in vehicle.c; consumers see
 * only the typedef. K26AstroVehicle is forward-declared by
 * k26astro_rt/scheduler.h via `struct K26AstroVehicle;`. */
typedef struct K26AstroVehicle K26AstroVehicle;

/* ---- Subsystem-handle forward declarations -------------------- *
 *
 * The concrete struct definitions live in the owning subsystem
 * library. The vehicle treats each handle as an opaque pointer.
 * Subsystem libraries are expected to publish a generation accessor
 * for their handle type so the vehicle's lifetime protocol can
 * detect a stale registration. */
struct K26AstroEngineCluster;
struct K26AstroTank;
struct K26AstroThermalNetwork;
struct K26AstroRCSQuad;
struct K26AstroEPSBus;
struct K26AstroAttitudeCtrl;
struct K26AstroSensor;
struct K26AstroMissionTimeline;
struct K26AstroAeroTable;
struct K26AstroPayload;

/* ---- Construction --------------------------------------------- */

/* Allocate a zero-initialised vehicle. basic_mass = 0; mga = 0;
 * inertia = identity (1 kg.m^2 on diagonal); com_offset = 0;
 * no body bound; no slots populated; no events; no lifecycle
 * windows. Returns NULL on allocation failure.
 *
 * Caller pairs with k26astro_vehicle_destroy. */
K26AstroVehicle *k26astro_vehicle_new(void);

/* Free a vehicle. Walks every composition slot and invokes the
 * subsystem notify shim so live subsystems can invalidate any
 * back-reference they hold. Frees the Ext attitude state's torque
 * registry. Frees stage-event context records previously handed to
 * any K26AstroGravState registry. The bound K26AstroBody is not
 * freed — the world owns the body's lifetime.
 *
 * Safe on NULL. */
void k26astro_vehicle_destroy(K26AstroVehicle *v);

/* Bind a body. The body becomes the vehicle's position / velocity /
 * attitude provider. The vehicle does not take ownership; the caller
 * keeps the body alive for the vehicle's lifetime. Re-binding is
 * allowed and silently replaces the previous binding. Safe on NULL
 * inputs. */
void k26astro_vehicle_bind_body(K26AstroVehicle *v, K26AstroBody *b);

/* ---- Mass setters --------------------------------------------- */

/* Set the current NASA basic mass. Initial value 0. Negative input
 * clamps to 0. Safe on NULL. */
void k26astro_vehicle_set_dry_mass(K26AstroVehicle *v, double basic_mass_kg);

/* Set the AIAA S-120A-2015 mass growth allowance. Initial value 0.
 * Negative input clamps to 0. Safe on NULL. */
void k26astro_vehicle_set_mga(K26AstroVehicle *v, double mga_kg);

/* ---- Inertia setters ------------------------------------------ */

/* Set inertia from a principal-axis (diagonal) triple. Constructs
 * the equivalent K26M3 with zero off-diagonals and routes through
 * k26astro_attitude_update_inertia (which recomputes the inverse).
 * Negative principal moments are rejected silently (a physical
 * inertia tensor is symmetric positive-definite). Safe on NULL. */
void k26astro_vehicle_set_inertia_diag(K26AstroVehicle *v,
                                       double ixx, double iyy, double izz);

/* Set the full 3x3 inertia tensor (body frame, kg.m^2). Assumed
 * symmetric positive-definite; symmetry is not enforced. Recomputes
 * the inverse via k26astro_attitude_update_inertia. Safe on NULL.
 *
 * Singular inertia (|det| ~ 0) is accepted by attitude_update_inertia
 * which zeroes the inverse; the vehicle is then unusable for torque
 * steps until inertia is re-set with a non-singular tensor. */
void k26astro_vehicle_set_inertia_full(K26AstroVehicle *v, K26M3 inertia);

/* ---- COM setter ----------------------------------------------- */

/* Set the centre-of-mass offset in the vehicle's body frame. The
 * offset is the displacement from the body-frame origin to the COM
 * — relevant for thrust-vector-induced torques about the COM. The
 * value is stored verbatim and returned by k26astro_vehicle_com_at;
 * propellant-distribution-driven COM drift couples in when tank
 * subsystems land. Safe on NULL. */
void k26astro_vehicle_set_com_offset(K26AstroVehicle *v,
                                     double x, double y, double z);

/* ---- Composition slots ---------------------------------------- *
 *
 * Each setter records a typed opaque handle in the vehicle. Handles
 * are non-owning. Lifetime contract:
 *
 *   - A subsystem either lives at least as long as its registration
 *     on the vehicle, or notifies the vehicle via
 *     k26astro_vehicle_invalidate_slot before its heap allocation
 *     disappears.
 *   - Vehicle destruction walks every live slot and calls the
 *     per-subsystem notify shim
 *     (k26astro_vehicle_notify_subsystem_destroy_) so subsystem
 *     destructors can blank any back-reference they hold.
 *
 * Each add returns the slot index >= 0 on success, -1 on NULL inputs
 * or allocation failure. */
typedef enum {
    K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER   = 1,
    K26ASTRO_VEHICLE_SLOT_TANK             = 2,
    K26ASTRO_VEHICLE_SLOT_THERMAL_NETWORK  = 3,
    K26ASTRO_VEHICLE_SLOT_RCS_QUAD         = 4,
    K26ASTRO_VEHICLE_SLOT_EPS_BUS          = 5,
    K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL    = 6,
    K26ASTRO_VEHICLE_SLOT_SENSOR           = 7,
    K26ASTRO_VEHICLE_SLOT_MISSION_TIMELINE = 8,
    K26ASTRO_VEHICLE_SLOT_AERO_TABLE       = 9,
    K26ASTRO_VEHICLE_SLOT_PAYLOAD          = 10,
    K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE   = 11
} K26AstroVehicleSlotKind;

int  k26astro_vehicle_add_engine_cluster (K26AstroVehicle *v,
                                          struct K26AstroEngineCluster *cluster);
int  k26astro_vehicle_add_tank           (K26AstroVehicle *v,
                                          struct K26AstroTank *tank);
int  k26astro_vehicle_add_thermal_network(K26AstroVehicle *v,
                                          struct K26AstroThermalNetwork *net);
int  k26astro_vehicle_add_rcs_quad       (K26AstroVehicle *v,
                                          struct K26AstroRCSQuad *quad);
int  k26astro_vehicle_add_eps_bus        (K26AstroVehicle *v,
                                          struct K26AstroEPSBus *bus);

/* Singleton — only one attitude controller per vehicle. Replaces any
 * previous controller (the previous occupant's notify shim fires
 * before replacement). Returns 0 on success, -1 on NULL vehicle. A
 * NULL ctrl clears the slot. */
int  k26astro_vehicle_set_attitude_ctrl  (K26AstroVehicle *v,
                                          struct K26AstroAttitudeCtrl *ctrl);

int  k26astro_vehicle_add_sensor         (K26AstroVehicle *v,
                                          struct K26AstroSensor *s);

/* Singleton — only one mission timeline per vehicle. Returns 0 on
 * success, -1 on NULL inputs, -2 if a timeline is already attached
 * (the previous occupant is not replaced; the caller must
 * invalidate first). */
int  k26astro_vehicle_add_mission_timeline(K26AstroVehicle *v,
                                           struct K26AstroMissionTimeline *tl);

/* Singleton — one aerodynamic coefficient table per vehicle. The
 * table opaque carries the aerodynamic reference lengths (S_ref,
 * c_ref, b_ref) and the Mach-α-β interpolation grid. Mirrors the
 * attitude controller singleton: replacing a live aero table fires
 * the previous occupant's notify shim before the swap, allowing
 * stage events that change the aerodynamic profile (heat-shield
 * jettison, second-stage drop) to swap tables mid-simulation.
 * Returns 0 on success, -1 on NULL vehicle. A NULL table clears
 * the slot. */
int  k26astro_vehicle_set_aero_table     (K26AstroVehicle *v,
                                          struct K26AstroAeroTable *t);

/* Payload mount — generic attachment surface for subsystems the
 * vehicle core does not model itself (instrument packages,
 * experiments, data recorders, deployables). A vehicle carries any
 * number of payloads; the slot array preserves registration order
 * so deterministic iteration sees payloads in a stable sequence.
 * Each payload carries its own mount geometry internally; the slot
 * record holds only the opaque pointer. Teardown dispatch is
 * per-handle via K26AstroPayload.on_owner_destroy (see payload.h) —
 * no per-kind vehicle hooks are involved, so any number of
 * payload-providing libraries can coexist in one link. Returns the
 * slot index >= 0 on success, -1 on NULL inputs or allocation
 * failure. */
int  k26astro_vehicle_add_payload        (K26AstroVehicle *v,
                                          struct K26AstroPayload *p);

/* Singleton — at most one unique payload per vehicle (a
 * vehicle-wide state package, e.g. a mission data recorder).
 * Replaces any previous occupant; the previous occupant's notify
 * shim fires before replacement. Returns 0 on success, -1 on NULL
 * vehicle. A NULL handle clears the slot. */
int  k26astro_vehicle_set_unique_payload (K26AstroVehicle *v,
                                          struct K26AstroPayload *p);

/* Subsystem-initiated invalidation. A subsystem destructor calls
 * this with its (kind, handle) before freeing itself. The vehicle
 * marks the slot dead in place (no compaction; slot indices stay
 * stable). Safe on dead handles, NULL vehicle, kind mismatch. */
void k26astro_vehicle_invalidate_slot(K26AstroVehicle *v,
                                      K26AstroVehicleSlotKind kind,
                                      void *handle);

/* ---- Queries -------------------------------------------------- */

/* Current basic mass plus any committed delta from past
 * commit_mass_step calls. In the absence of registered engine-
 * cluster callbacks this equals basic_mass_kg. Returns 0.0 on
 * NULL. */
double k26astro_vehicle_mass_now      (const K26AstroVehicle *v);

/* basic_mass + mga_kg. Returns 0.0 on NULL. */
double k26astro_vehicle_predicted_mass(const K26AstroVehicle *v);

/* Propellant-sizing mass (MEV). Returns predicted_mass; an explicit
 * reserve term is not currently modelled. Returns 0.0 on NULL. */
double k26astro_vehicle_mev_mass      (const K26AstroVehicle *v);

/* Walks the sorted stage-event timeline and applies every event
 * with epoch <= t. Returns the resulting basic mass (clamped to
 * >= 0). Returns 0.0 on NULL. */
double k26astro_vehicle_mass_at       (const K26AstroVehicle *v,
                                       K26AstroEpoch t);

/* Returns the static com_offset. Propellant-distribution coupling
 * across registered tanks is not currently modelled — the value is
 * whatever the caller most recently set via _set_com_offset.
 * Returns {0, 0, 0} on NULL. */
K26V3 k26astro_vehicle_com_at         (const K26AstroVehicle *v,
                                       K26AstroEpoch t);

/* Walks the stage-event timeline applying each event's post-event
 * diagonal moments to a working K26M3 (initialised from the current
 * vehicle inertia tensor). Off-diagonals are preserved across stage
 * events. Returns the identity matrix on NULL. */
K26M3 k26astro_vehicle_inertia_at     (const K26AstroVehicle *v,
                                       K26AstroEpoch t);

/* ---- Internal accessors for cross-lib wiring ------------------ *
 *
 * Consumed by mass_step.c, stage_event.c, lifecycle.c, and by
 * subsystem libraries that build on the composition root. Not part
 * of the kflbi surface. */

/* Returns a pointer to the vehicle's Ext attitude state. NULL on
 * NULL vehicle. The returned pointer's lifetime is the vehicle's. */
K26AstroAttitudeStateExt *
k26astro_vehicle_attitude_ext(K26AstroVehicle *v);

/* Returns the bound body, or NULL if no body is bound. */
K26AstroBody *k26astro_vehicle_body(K26AstroVehicle *v);

/* Per-substep mass-flow accumulator. Engine-cluster thrust callbacks
 * add dot_m via _add (negative for mass leaving the vehicle);
 * k26astro_vehicle_commit_mass_step consumes the accumulated value
 * at substep close. Initial value 0. */
double k26astro_vehicle_mass_accum_get  (const K26AstroVehicle *v);
void   k26astro_vehicle_mass_accum_add  (K26AstroVehicle *v, double dot_m_kg_per_s);
void   k26astro_vehicle_mass_accum_clear(K26AstroVehicle *v);

/* Vehicle generation counter. Increments on destroy. Subsystems may
 * cache the generation at registration time and check it later to
 * detect a destroyed owner. */
uint64_t k26astro_vehicle_generation(const K26AstroVehicle *v);

/* Per-substep mass commit. Applies accumulated dot_m * dt to
 * basic_mass and propagates to the bound body via
 * k26astro_body_set_mass. Clears the accumulator. Safe on NULL. */
void k26astro_vehicle_commit_mass_step(K26AstroVehicle *v, double dt);

/* ---- Non-gravitational acceleration cache --------------------- *
 *
 * Most-recent inertial-frame non-gravitational acceleration on the
 * vehicle's bound body. Populated by a perturb-capture callback that
 * the runtime can install via the world's enable-imu-accel-cache
 * helper; consumed by inertial-measurement-unit sensors to derive a
 * body-frame felt-acceleration reading.
 *
 * Defaults to {0, 0, 0} on vehicle construction. The setter is
 * also a public API to allow direct population for fixtures or
 * tests that do not wire a runtime capture path.
 *
 * Sensor consumers should be aware that perturb callbacks registered
 * on the world's gravity state after the cache is enabled will not
 * contribute to the captured sum unless the cache enable call is
 * re-issued. */
K26V3 k26astro_vehicle_last_non_grav_accel_inertial(
    const K26AstroVehicle *v);

void k26astro_vehicle_set_last_non_grav_accel_inertial(
    K26AstroVehicle *v, K26V3 a);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_VEHICLE_H */
