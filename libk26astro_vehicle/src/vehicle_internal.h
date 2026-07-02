/* vehicle_internal.h — private struct layout shared between
 * vehicle.c, composition.c, stage_event.c, lifecycle.c, and
 * mass_step.c. Not installed; not part of the public ABI. */
#ifndef K26ASTRO_VEHICLE_INTERNAL_H
#define K26ASTRO_VEHICLE_INTERNAL_H

#include <stdint.h>

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"
#include "k26astro_body/attitude.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Per-slot record ------------------------------------------ *
 *
 * Slot indices are stable across the vehicle's life: invalidation
 * flips `live` to 0 without compacting the array. The `gen` field
 * is reserved for a subsystem-published generation counter; until
 * a subsystem library publishes a `_get_generation` accessor the
 * slot records 0 at registration time. */
typedef struct {
    void                      *handle;   /* opaque subsystem pointer */
    uint64_t                   gen;      /* subsystem generation snapshot */
    K26AstroVehicleSlotKind    kind;
    int                        live;     /* 1 = registered; 0 = invalidated */
} K26AstroVehicleSlot_;

/* ---- Stage-event record --------------------------------------- *
 *
 * Stored sorted ascending by epoch. The `fired` flag is set by the
 * gravity handler when it runs at the bisected event epoch. */
typedef struct {
    K26AstroEpoch  epoch;
    double         mass_drop_kg;
    double         ixx_after;
    double         iyy_after;
    double         izz_after;
    int            fired;
} K26AstroVehicleStageEvent_;

/* ---- Stage-event handler context ------------------------------ *
 *
 * Heap-allocated by k26astro_vehicle_register_stage_events_with;
 * the vehicle keeps a parallel pointer array so destroy can free
 * every context that was handed to a gravity registry. */
typedef struct {
    K26AstroVehicle             *vehicle;
    K26AstroVehicleStageEvent_  *event;
    double                       scheduled_t_s;
} K26AstroVehicleEventCtx_;

/* ---- Lifecycle window ----------------------------------------- *
 *
 * Stored in registration order; is_active_at scans back-to-front so
 * the latest-registered window containing t wins on overlap. */
typedef enum {
    K26ASTRO_LC_ACTIVE = 1,
    K26ASTRO_LC_RAILS  = 2
} K26AstroLifecycleKind_;

typedef struct {
    K26AstroEpoch           t_start;
    K26AstroEpoch           t_end;
    K26AstroLifecycleKind_  kind;
} K26AstroVehicleLifecycleWindow_;

/* ---- The vehicle struct --------------------------------------- */

struct K26AstroVehicle {
    /* Mass */
    double  basic_mass_kg;
    double  mga_kg;
    double  mass_accum;            /* per-substep dot_m kg/s accumulator */

    /* Inertia (full 3x3 via Ext attitude state — embedded by value) */
    K26AstroAttitudeStateExt attitude_ext;

    /* COM offset (body frame) */
    K26V3   com_offset;

    /* Bound body (non-owning) */
    K26AstroBody *body;

    /* Composition slots — flex array per kind for cache-friendly
     * per-kind iteration. The attitude controller is a singleton. */
    K26AstroVehicleSlot_ *slots_engine;        int n_engine,        cap_engine;
    K26AstroVehicleSlot_ *slots_tank;          int n_tank,          cap_tank;
    K26AstroVehicleSlot_ *slots_thermal;       int n_thermal,       cap_thermal;
    K26AstroVehicleSlot_ *slots_rcs;           int n_rcs,           cap_rcs;
    K26AstroVehicleSlot_ *slots_eps;           int n_eps,           cap_eps;
    K26AstroVehicleSlot_ *slots_sensor;        int n_sensor,        cap_sensor;
    K26AstroVehicleSlot_ *slots_payload;       int n_payload,       cap_payload;
    K26AstroVehicleSlot_  slot_attitude;       /* singleton */
    K26AstroVehicleSlot_  slot_mission;        /* singleton */
    K26AstroVehicleSlot_  slot_aero_table;     /* singleton */
    K26AstroVehicleSlot_  slot_payload_unique; /* singleton */

    /* Stage-event timeline (sorted) + parallel ctx pointer pool */
    K26AstroVehicleStageEvent_  *events;        int n_events,  cap_events;
    K26AstroVehicleEventCtx_   **event_ctxs;    int n_ctxs,    cap_ctxs;

    /* Lifecycle windows (registration order) */
    K26AstroVehicleLifecycleWindow_ *windows;
    int                              n_windows, cap_windows;

    /* Generation counter — increments on destroy. */
    uint64_t generation;

    /* Inertial-frame non-gravitational acceleration cache. Populated
     * by an opt-in perturb-capture callback installed via
     * k26astro_world_enable_imu_accel_cache. Read by inertial-
     * measurement-unit sensors via
     * k26astro_vehicle_last_non_grav_accel_inertial. */
    K26V3 last_non_grav_accel_inertial;
};

/* ---- Internal cross-file declarations ------------------------- */

/* Implemented in composition.c. Walks every live slot and invokes
 * the per-subsystem notify shim. Called by k26astro_vehicle_destroy. */
void k26astro_vehicle_notify_all_slots_(K26AstroVehicle *v);

/* Weak default shim implemented in composition.c. Subsystem
 * libraries — and the subsystem-outlive test — provide strong
 * overrides at link time. */
void k26astro_vehicle_notify_subsystem_destroy_(
    void *handle, K26AstroVehicleSlotKind kind);

/* Per-kind destroy hooks. Each subsystem library provides a strong
 * override for the kind it owns; libs that do not care about a kind
 * leave the weak no-op default in place. The dispatcher in
 * k26astro_vehicle_notify_subsystem_destroy_ routes by slot kind, so
 * multiple subsystem libraries can co-link without symbol conflicts. */
void k26astro_vehicle_on_engine_cluster_destroy_ (void *handle);
void k26astro_vehicle_on_tank_destroy_           (void *handle);
void k26astro_vehicle_on_thermal_network_destroy_(void *handle);
void k26astro_vehicle_on_rcs_quad_destroy_       (void *handle);
void k26astro_vehicle_on_eps_bus_destroy_        (void *handle);
void k26astro_vehicle_on_attitude_ctrl_destroy_  (void *handle);
void k26astro_vehicle_on_sensor_destroy_         (void *handle);
void k26astro_vehicle_on_mission_timeline_destroy_(void *handle);
void k26astro_vehicle_on_aero_table_destroy_     (void *handle);

/* PAYLOAD slots have no per-kind hooks: dispatch is per-handle via
 * K26AstroPayload.on_owner_destroy inside the default notify shim. */

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_INTERNAL_H */
