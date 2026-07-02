/* test_vehicle_subsystem_outlive.c — lifetime-protocol gate.
 *
 * Validates the dual-direction back-reference protocol that closes
 * the use-after-free window between vehicle and subsystem handles.
 *
 * Scenarios:
 *   (1) Vehicle destroyed before subsystem: the destroy-time notify
 *       shim must fire on the live slot and the stub subsystem must
 *       observe its owner pointer cleared.
 *   (2) Subsystem destroyed before vehicle: subsystem calls
 *       k26astro_vehicle_invalidate_slot before freeing itself; the
 *       vehicle's subsequent destroy must skip the dead slot
 *       without dereferencing the freed handle.
 *   (3) Both alive, vehicle destroyed first: equivalent to (1) but
 *       confirms the post-destroy subsystem cleanup is clean.
 *
 * The test defines a stub K26AstroEngineCluster struct inline and
 * provides a strong override of
 * k26astro_vehicle_notify_subsystem_destroy_ — the weak default in
 * composition.c is shadowed at link time. Run under ASan/LSan/UBSan. */

#include "k26astro_vehicle/vehicle.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Stub subsystem — the real struct lives in the owning subsystem
 * library; this is just enough state to exercise the protocol. */
struct K26AstroEngineCluster {
    uint64_t          generation;
    K26AstroVehicle  *owner;
    uint64_t          owner_gen_at_register;
    int               destructor_ran;
};

/* Track notify-shim invocations for assertions below. */
static int notify_calls_;
static int notify_last_kind_;

/* Strong override of the weak default in composition.c. */
void k26astro_vehicle_notify_subsystem_destroy_(void *handle,
                                                K26AstroVehicleSlotKind kind)
{
    notify_calls_++;
    notify_last_kind_ = (int)kind;
    if (kind != K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER) return;
    struct K26AstroEngineCluster *c = (struct K26AstroEngineCluster *)handle;
    if (c) c->owner = NULL;
}

int main(void)
{
    /* ---- Scenario 1: vehicle destroyed before subsystem. ------ */
    notify_calls_     = 0;
    notify_last_kind_ = 0;

    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);
    struct K26AstroEngineCluster *c = calloc(1, sizeof(*c));
    assert(c != NULL);
    c->owner                 = v;
    c->owner_gen_at_register = k26astro_vehicle_generation(v);

    int slot = k26astro_vehicle_add_engine_cluster(v, c);
    assert(slot >= 0);

    /* Tear down the vehicle. The notify shim fires on the live slot. */
    k26astro_vehicle_destroy(v);
    assert(notify_calls_     == 1);
    assert(notify_last_kind_ == K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER);
    assert(c->owner == NULL);
    free(c);

    /* ---- Scenario 2: subsystem destroyed first. -------------- */
    notify_calls_ = 0;
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    struct K26AstroEngineCluster *c2 = calloc(1, sizeof(*c2));
    c2->owner = v2;
    int slot2 = k26astro_vehicle_add_engine_cluster(v2, c2);
    assert(slot2 >= 0);

    /* Subsystem-initiated invalidation before freeing. */
    k26astro_vehicle_invalidate_slot(v2,
        K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER, c2);
    free(c2);

    /* Vehicle destruction must skip the dead slot — no notify call,
     * no dereference of the freed pointer. */
    k26astro_vehicle_destroy(v2);
    assert(notify_calls_ == 0);

    /* ---- Scenario 3: both alive, vehicle destroyed first. ---- */
    notify_calls_ = 0;
    K26AstroVehicle *v3 = k26astro_vehicle_new();
    struct K26AstroEngineCluster *c3 = calloc(1, sizeof(*c3));
    c3->owner = v3;
    assert(k26astro_vehicle_add_engine_cluster(v3, c3) >= 0);

    k26astro_vehicle_destroy(v3);
    assert(notify_calls_ == 1);
    assert(c3->owner == NULL);
    /* Subsystem cleanup post-vehicle — clean, owner already nulled. */
    free(c3);

    /* ---- Scenario 4: invalidate on a kind / handle that isn't
     *                  registered — must be a no-op. -------------- */
    K26AstroVehicle *v4 = k26astro_vehicle_new();
    struct K26AstroEngineCluster *fake = (struct K26AstroEngineCluster *)0xDEADBEEF;
    k26astro_vehicle_invalidate_slot(v4,
        K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER, fake);
    /* No registered slots so destroy walks zero entries. */
    notify_calls_ = 0;
    k26astro_vehicle_destroy(v4);
    assert(notify_calls_ == 0);

    /* ---- Scenario 5: add_engine_cluster with NULL handle -> -1 -- */
    K26AstroVehicle *v5 = k26astro_vehicle_new();
    assert(k26astro_vehicle_add_engine_cluster(v5, NULL) == -1);
    /* NULL vehicle: also -1. */
    assert(k26astro_vehicle_add_engine_cluster(NULL, NULL) == -1);
    k26astro_vehicle_destroy(v5);

    /* ---- Scenario 6: NULL inputs to invalidate_slot are safe -- */
    k26astro_vehicle_invalidate_slot(NULL,
        K26ASTRO_VEHICLE_SLOT_ENGINE_CLUSTER, NULL);

    printf("test_vehicle_subsystem_outlive: OK\n");
    return 0;
}
