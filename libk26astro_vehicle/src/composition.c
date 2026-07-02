/* composition.c — typed-opaque slot management for vehicle
 * subsystems and the destroy-time notify protocol.
 *
 * Slot bookkeeping uses one flex array per slot kind for
 * cache-friendly per-kind iteration. Insertion grows the array
 * (doubling realloc from capacity 4); invalidation flips a `live`
 * flag in place so slot indices stay stable across the vehicle's
 * lifetime.
 *
 * Subsystem lifetime protocol:
 *
 *   - vehicle_add_<kind>(v, handle): record the (handle, kind)
 *     pair with live = 1 and a generation counter snapshot. The
 *     generation field is reserved for a per-subsystem accessor;
 *     until a subsystem library publishes one, the slot records 0.
 *
 *   - vehicle_destroy(v): walks every live slot and invokes the
 *     per-subsystem notify shim. The default no-op shim is defined
 *     here; a subsystem library overrides it with a strong
 *     definition; tests provide a strong override locally to
 *     verify the call-through.
 *
 *   - vehicle_invalidate_slot(v, kind, handle): subsystem
 *     destructors call this before freeing themselves, marking the
 *     slot dead so the vehicle no longer dereferences the freed
 *     pointer on subsequent walks.
 *
 * The dual-direction back-reference is what closes the use-after-
 * free window. The test_vehicle_subsystem_outlive gate exercises
 * the protocol with a stub subsystem inline. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/payload.h"
#include "vehicle_internal.h"

#include <stdlib.h>

/* ---- Weak default per-kind hooks ----------------------------- *
 *
 * A subsystem library overrides exactly the hook for the kind it
 * owns. Hooks for unowned kinds stay weak no-ops; no link conflict
 * arises when multiple subsystem libs coexist in the same binary
 * because each lib publishes a different symbol. */
__attribute__((weak)) void
k26astro_vehicle_on_engine_cluster_destroy_(void *h)  { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_tank_destroy_(void *h)            { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_thermal_network_destroy_(void *h) { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_rcs_quad_destroy_(void *h)        { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_eps_bus_destroy_(void *h)         { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_attitude_ctrl_destroy_(void *h)   { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_sensor_destroy_(void *h)          { (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_mission_timeline_destroy_(void *h){ (void)h; }
__attribute__((weak)) void
k26astro_vehicle_on_aero_table_destroy_(void *h)      { (void)h; }

/* PAYLOAD slots deliberately have no per-kind hooks. Payload
 * teardown dispatches per-handle through the on_owner_destroy
 * callback carried in the K26AstroPayload base (payload.h), so any
 * number of payload-providing libraries can co-link without symbol
 * conflicts and without touching this file. */
static void payload_notify_(void *handle)
{
    struct K26AstroPayload *p = (struct K26AstroPayload *)handle;
    if (!p) return;
    if (p->on_owner_destroy) {
        p->on_owner_destroy(p);
    }
    p->owner      = NULL;
    p->owner_slot = -1;
}

/* ---- Weak default notify shim -------------------------------- *
 *
 * Dispatcher: routes (handle, kind) to the appropriate per-kind
 * hook. The dispatcher stays weak so the subsystem-outlive test's
 * strong override of this symbol continues to shadow it (the
 * per-kind hooks are never reached in that test's call graph). */
__attribute__((weak))
void k26astro_vehicle_notify_subsystem_destroy_(void *handle,
                                                K26AstroVehicleSlotKind kind)
{
    switch (kind) {
    case K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER:
        k26astro_vehicle_on_engine_cluster_destroy_(handle);  break;
    case K26ASTRO_VEHICLE_SLOT_TANK:
        k26astro_vehicle_on_tank_destroy_(handle);            break;
    case K26ASTRO_VEHICLE_SLOT_THERMAL_NETWORK:
        k26astro_vehicle_on_thermal_network_destroy_(handle); break;
    case K26ASTRO_VEHICLE_SLOT_RCS_QUAD:
        k26astro_vehicle_on_rcs_quad_destroy_(handle);        break;
    case K26ASTRO_VEHICLE_SLOT_EPS_BUS:
        k26astro_vehicle_on_eps_bus_destroy_(handle);         break;
    case K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL:
        k26astro_vehicle_on_attitude_ctrl_destroy_(handle);   break;
    case K26ASTRO_VEHICLE_SLOT_SENSOR:
        k26astro_vehicle_on_sensor_destroy_(handle);          break;
    case K26ASTRO_VEHICLE_SLOT_MISSION_TIMELINE:
        k26astro_vehicle_on_mission_timeline_destroy_(handle); break;
    case K26ASTRO_VEHICLE_SLOT_AERO_TABLE:
        k26astro_vehicle_on_aero_table_destroy_(handle);      break;
    case K26ASTRO_VEHICLE_SLOT_PAYLOAD:
    case K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE:
        payload_notify_(handle);                              break;
    }
}

/* ---- Walk every live slot ------------------------------------- */

static void notify_array_(K26AstroVehicleSlot_ *arr, int n,
                          K26AstroVehicleSlotKind kind)
{
    for (int i = 0; i < n; i++) {
        if (arr[i].live) {
            k26astro_vehicle_notify_subsystem_destroy_(arr[i].handle, kind);
        }
    }
}

void k26astro_vehicle_notify_all_slots_(K26AstroVehicle *v)
{
    if (!v) return;
    notify_array_(v->slots_engine,  v->n_engine,
                  K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER);
    notify_array_(v->slots_tank,    v->n_tank,
                  K26ASTRO_VEHICLE_SLOT_TANK);
    notify_array_(v->slots_thermal, v->n_thermal,
                  K26ASTRO_VEHICLE_SLOT_THERMAL_NETWORK);
    notify_array_(v->slots_rcs,     v->n_rcs,
                  K26ASTRO_VEHICLE_SLOT_RCS_QUAD);
    notify_array_(v->slots_eps,     v->n_eps,
                  K26ASTRO_VEHICLE_SLOT_EPS_BUS);
    notify_array_(v->slots_sensor,  v->n_sensor,
                  K26ASTRO_VEHICLE_SLOT_SENSOR);
    notify_array_(v->slots_payload, v->n_payload,
                  K26ASTRO_VEHICLE_SLOT_PAYLOAD);
    if (v->slot_attitude.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_attitude.handle,
            K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL);
    }
    if (v->slot_mission.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_mission.handle,
            K26ASTRO_VEHICLE_SLOT_MISSION_TIMELINE);
    }
    if (v->slot_aero_table.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_aero_table.handle,
            K26ASTRO_VEHICLE_SLOT_AERO_TABLE);
    }
    if (v->slot_payload_unique.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_payload_unique.handle,
            K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE);
    }
}

/* ---- Generic slot append helper ------------------------------- */

static int slot_array_append_(K26AstroVehicleSlot_ **arr,
                              int *n, int *cap,
                              void *handle,
                              K26AstroVehicleSlotKind kind)
{
    if (*n >= *cap) {
        int new_cap = *cap ? (*cap * 2) : 4;
        K26AstroVehicleSlot_ *grown = (K26AstroVehicleSlot_ *)realloc(
            *arr, (size_t)new_cap * sizeof(K26AstroVehicleSlot_));
        if (!grown) return -1;
        *arr = grown;
        *cap = new_cap;
    }
    (*arr)[*n].handle = handle;
    (*arr)[*n].gen    = 0;
    (*arr)[*n].kind   = kind;
    (*arr)[*n].live   = 1;
    return (*n)++;
}

/* ---- Per-kind add functions ----------------------------------- */

int k26astro_vehicle_add_engine_cluster(K26AstroVehicle *v,
                                        struct K26AstroEngineCluster *c)
{
    if (!v || !c) return -1;
    return slot_array_append_(&v->slots_engine, &v->n_engine, &v->cap_engine,
                              c, K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER);
}

int k26astro_vehicle_add_tank(K26AstroVehicle *v,
                              struct K26AstroTank *t)
{
    if (!v || !t) return -1;
    return slot_array_append_(&v->slots_tank, &v->n_tank, &v->cap_tank,
                              t, K26ASTRO_VEHICLE_SLOT_TANK);
}

int k26astro_vehicle_add_thermal_network(K26AstroVehicle *v,
                                         struct K26AstroThermalNetwork *n)
{
    if (!v || !n) return -1;
    return slot_array_append_(&v->slots_thermal, &v->n_thermal, &v->cap_thermal,
                              n, K26ASTRO_VEHICLE_SLOT_THERMAL_NETWORK);
}

int k26astro_vehicle_add_rcs_quad(K26AstroVehicle *v,
                                  struct K26AstroRCSQuad *q)
{
    if (!v || !q) return -1;
    return slot_array_append_(&v->slots_rcs, &v->n_rcs, &v->cap_rcs,
                              q, K26ASTRO_VEHICLE_SLOT_RCS_QUAD);
}

int k26astro_vehicle_add_eps_bus(K26AstroVehicle *v,
                                 struct K26AstroEPSBus *b)
{
    if (!v || !b) return -1;
    return slot_array_append_(&v->slots_eps, &v->n_eps, &v->cap_eps,
                              b, K26ASTRO_VEHICLE_SLOT_EPS_BUS);
}

int k26astro_vehicle_add_sensor(K26AstroVehicle *v,
                                struct K26AstroSensor *s)
{
    if (!v || !s) return -1;
    return slot_array_append_(&v->slots_sensor, &v->n_sensor, &v->cap_sensor,
                              s, K26ASTRO_VEHICLE_SLOT_SENSOR);
}

int k26astro_vehicle_add_payload(K26AstroVehicle *v,
                                 struct K26AstroPayload *p)
{
    if (!v || !p) return -1;
    return slot_array_append_(&v->slots_payload, &v->n_payload,
                              &v->cap_payload,
                              p, K26ASTRO_VEHICLE_SLOT_PAYLOAD);
}

int k26astro_vehicle_set_unique_payload(K26AstroVehicle *v,
                                        struct K26AstroPayload *p)
{
    if (!v) return -1;
    /* Notify any previous occupant before replacement — mirrors the
     * attitude-controller singleton's replace-with-notify semantics. */
    if (v->slot_payload_unique.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_payload_unique.handle,
            K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE);
    }
    v->slot_payload_unique.handle = p;
    v->slot_payload_unique.kind   = K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE;
    v->slot_payload_unique.gen    = 0;
    v->slot_payload_unique.live   = (p != NULL) ? 1 : 0;
    return 0;
}

/* The mission-timeline back-link setter is defined inside
 * libk26astro_mission as a non-static helper; the vehicle has no
 * dependency on the mission library, so the call is wired by a
 * forward declaration here. The mission lib is the only place that
 * publishes a strong definition. */
extern void k26astro_mission_timeline_set_bound_vehicle_(
    struct K26AstroMissionTimeline *tl,
    struct K26AstroVehicle         *v);
__attribute__((weak)) void
k26astro_mission_timeline_set_bound_vehicle_(
    struct K26AstroMissionTimeline *tl,
    struct K26AstroVehicle         *v)
{
    (void)tl; (void)v;
}

int k26astro_vehicle_add_mission_timeline(K26AstroVehicle *v,
                                          struct K26AstroMissionTimeline *tl)
{
    if (!v || !tl) return -1;
    if (v->slot_mission.live) return -2;
    v->slot_mission.handle = tl;
    v->slot_mission.kind   = K26ASTRO_VEHICLE_SLOT_MISSION_TIMELINE;
    v->slot_mission.gen    = 0;
    v->slot_mission.live   = 1;
    k26astro_mission_timeline_set_bound_vehicle_(tl, v);
    return 0;
}

int k26astro_vehicle_set_attitude_ctrl(K26AstroVehicle *v,
                                       struct K26AstroAttitudeCtrl *ctrl)
{
    if (!v) return -1;
    /* Notify any previous occupant before replacement. */
    if (v->slot_attitude.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_attitude.handle,
            K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL);
    }
    v->slot_attitude.handle = ctrl;
    v->slot_attitude.kind   = K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL;
    v->slot_attitude.gen    = 0;
    v->slot_attitude.live   = (ctrl != NULL) ? 1 : 0;
    return 0;
}

int k26astro_vehicle_set_aero_table(K26AstroVehicle *v,
                                    struct K26AstroAeroTable *t)
{
    if (!v) return -1;
    /* Notify any previous occupant before replacement. Stage events
     * that change the aerodynamic profile (heat-shield jettison,
     * second-stage drop) swap tables mid-simulation through this
     * same code path. */
    if (v->slot_aero_table.live) {
        k26astro_vehicle_notify_subsystem_destroy_(
            v->slot_aero_table.handle,
            K26ASTRO_VEHICLE_SLOT_AERO_TABLE);
    }
    v->slot_aero_table.handle = t;
    v->slot_aero_table.kind   = K26ASTRO_VEHICLE_SLOT_AERO_TABLE;
    v->slot_aero_table.gen    = 0;
    v->slot_aero_table.live   = (t != NULL) ? 1 : 0;
    return 0;
}

/* ---- Subsystem-initiated invalidation ------------------------- */

static void invalidate_in_(K26AstroVehicleSlot_ *arr, int n, void *handle)
{
    for (int i = 0; i < n; i++) {
        if (arr[i].live && arr[i].handle == handle) {
            arr[i].live   = 0;
            arr[i].handle = NULL;
            return;
        }
    }
}

void k26astro_vehicle_invalidate_slot(K26AstroVehicle *v,
                                      K26AstroVehicleSlotKind kind,
                                      void *handle)
{
    if (!v) return;
    switch (kind) {
        case K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER:
            invalidate_in_(v->slots_engine, v->n_engine, handle); break;
        case K26ASTRO_VEHICLE_SLOT_TANK:
            invalidate_in_(v->slots_tank, v->n_tank, handle); break;
        case K26ASTRO_VEHICLE_SLOT_THERMAL_NETWORK:
            invalidate_in_(v->slots_thermal, v->n_thermal, handle); break;
        case K26ASTRO_VEHICLE_SLOT_RCS_QUAD:
            invalidate_in_(v->slots_rcs, v->n_rcs, handle); break;
        case K26ASTRO_VEHICLE_SLOT_EPS_BUS:
            invalidate_in_(v->slots_eps, v->n_eps, handle); break;
        case K26ASTRO_VEHICLE_SLOT_SENSOR:
            invalidate_in_(v->slots_sensor, v->n_sensor, handle); break;
        case K26ASTRO_VEHICLE_SLOT_ATTITUDE_CTRL:
            if (v->slot_attitude.live && v->slot_attitude.handle == handle) {
                v->slot_attitude.live   = 0;
                v->slot_attitude.handle = NULL;
            }
            break;
        case K26ASTRO_VEHICLE_SLOT_MISSION_TIMELINE:
            if (v->slot_mission.live && v->slot_mission.handle == handle) {
                v->slot_mission.live   = 0;
                v->slot_mission.handle = NULL;
            }
            break;
        case K26ASTRO_VEHICLE_SLOT_AERO_TABLE:
            if (v->slot_aero_table.live && v->slot_aero_table.handle == handle) {
                v->slot_aero_table.live   = 0;
                v->slot_aero_table.handle = NULL;
            }
            break;
        case K26ASTRO_VEHICLE_SLOT_PAYLOAD:
            invalidate_in_(v->slots_payload, v->n_payload, handle); break;
        case K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE:
            if (v->slot_payload_unique.live && v->slot_payload_unique.handle == handle) {
                v->slot_payload_unique.live   = 0;
                v->slot_payload_unique.handle = NULL;
            }
            break;
    }
}
