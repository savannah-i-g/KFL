/* test_vehicle_payload_slots.c — lifetime-protocol gate for the two
 * generic PAYLOAD composition slots (array + unique singleton).
 *
 * Exercises the same dual-direction back-reference protocol as
 * test_vehicle_subsystem_outlive does for engine_cluster, but through
 * the per-handle dispatch path: payloads carry their own
 * on_owner_destroy callback in the K26AstroPayload base, and the
 * DEFAULT notify shim routes to it — so unlike the per-kind hook
 * tests, this test deliberately does NOT override
 * k26astro_vehicle_notify_subsystem_destroy_. The unique slot's
 * replace-with-notify semantics mirror attitude_ctrl. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/payload.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* ---- Concrete payload stub ------------------------------------- */

#define TEST_PAYLOAD_KIND_A 0x54500001u
#define TEST_PAYLOAD_KIND_B 0x54500002u

struct TestPayload {
    K26AstroPayload base;   /* first member; cast contract */
    int             owner_destroy_ran;
};

static int callback_total_;

static void on_owner_destroy_(K26AstroPayload *p)
{
    callback_total_++;
    ((struct TestPayload *)p)->owner_destroy_ran = 1;
}

static struct TestPayload *payload_new_(uint32_t kind, int with_callback)
{
    struct TestPayload *t = calloc(1, sizeof(*t));
    t->base.kind             = kind;
    t->base.owner_slot       = -1;
    t->base.on_owner_destroy = with_callback ? on_owner_destroy_ : NULL;
    return t;
}

/* ---- Array slot -------------------------------------------------- */

static void test_payload_array_slot_(void)
{
    /* Two payloads on one vehicle, destroy vehicle first. Both
     * callbacks fire; both stubs see their owner cleared. */
    callback_total_ = 0;
    K26AstroVehicle *v = k26astro_vehicle_new();
    struct TestPayload *p1 = payload_new_(TEST_PAYLOAD_KIND_A, 1);
    struct TestPayload *p2 = payload_new_(TEST_PAYLOAD_KIND_B, 1);
    assert(k26astro_payload_bind(&p1->base, v) == 0);
    assert(k26astro_payload_bind(&p2->base, v) == 0);
    assert(p1->base.owner == v && p1->base.owner_slot == 0);
    assert(p2->base.owner == v && p2->base.owner_slot == 1);
    assert(k26astro_payload_kind(&p1->base) == TEST_PAYLOAD_KIND_A);

    k26astro_vehicle_destroy(v);
    assert(callback_total_ == 2);
    assert(p1->owner_destroy_ran && p1->base.owner == NULL);
    assert(p2->owner_destroy_ran && p2->base.owner == NULL);
    free(p1);
    free(p2);

    /* NULL callback: owner still blanked, nothing else happens. */
    callback_total_ = 0;
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    struct TestPayload *p3 = payload_new_(TEST_PAYLOAD_KIND_A, 0);
    assert(k26astro_payload_bind(&p3->base, v2) == 0);
    k26astro_vehicle_destroy(v2);
    assert(callback_total_ == 0);
    assert(p3->base.owner == NULL);
    free(p3);

    /* Payload-side unlink skips the notify path. */
    callback_total_ = 0;
    K26AstroVehicle *v3 = k26astro_vehicle_new();
    struct TestPayload *p4 = payload_new_(TEST_PAYLOAD_KIND_A, 1);
    struct TestPayload *p5 = payload_new_(TEST_PAYLOAD_KIND_B, 1);
    assert(k26astro_payload_bind(&p4->base, v3) == 0);
    assert(k26astro_payload_bind(&p5->base, v3) == 0);
    k26astro_payload_unlink_(&p4->base);
    assert(p4->base.owner == NULL && p4->base.owner_slot == -1);
    free(p4);
    k26astro_vehicle_destroy(v3);
    assert(callback_total_ == 1);           /* only p5 */
    assert(p5->owner_destroy_ran);
    free(p5);

    /* NULL handle / NULL vehicle return -1. */
    K26AstroVehicle *v4 = k26astro_vehicle_new();
    assert(k26astro_vehicle_add_payload(v4, NULL) == -1);
    assert(k26astro_vehicle_add_payload(NULL, NULL) == -1);
    assert(k26astro_payload_bind(NULL, v4) == -1);
    k26astro_vehicle_destroy(v4);
}

/* ---- Unique singleton slot -------------------------------------- */

static void test_payload_unique_slot_(void)
{
    /* Replace-with-notify: setting a second payload fires the
     * previous occupant's callback before swap. */
    callback_total_ = 0;
    K26AstroVehicle *v = k26astro_vehicle_new();
    struct TestPayload *u1 = payload_new_(TEST_PAYLOAD_KIND_A, 1);
    struct TestPayload *u2 = payload_new_(TEST_PAYLOAD_KIND_A, 1);

    assert(k26astro_payload_bind_unique(&u1->base, v) == 0);
    assert(u1->base.owner_slot == K26ASTRO_PAYLOAD_SLOT_UNIQUE);
    assert(callback_total_ == 0);

    assert(k26astro_payload_bind_unique(&u2->base, v) == 0);
    assert(callback_total_ == 1);
    assert(u1->owner_destroy_ran && u1->base.owner == NULL);
    free(u1);

    /* Vehicle destroy notifies the live occupant. */
    k26astro_vehicle_destroy(v);
    assert(callback_total_ == 2);
    assert(u2->base.owner == NULL);
    free(u2);

    /* NULL handle evicts the occupant (with notify) and clears. */
    callback_total_ = 0;
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    struct TestPayload *u3 = payload_new_(TEST_PAYLOAD_KIND_A, 1);
    assert(k26astro_payload_bind_unique(&u3->base, v2) == 0);
    assert(k26astro_vehicle_set_unique_payload(v2, NULL) == 0);
    assert(callback_total_ == 1);
    assert(u3->base.owner == NULL);
    free(u3);
    /* No live occupant — destroy skips the notify. */
    k26astro_vehicle_destroy(v2);
    assert(callback_total_ == 1);

    /* Payload-initiated unlink on the singleton. */
    callback_total_ = 0;
    K26AstroVehicle *v3 = k26astro_vehicle_new();
    struct TestPayload *u4 = payload_new_(TEST_PAYLOAD_KIND_A, 1);
    assert(k26astro_payload_bind_unique(&u4->base, v3) == 0);
    k26astro_payload_unlink_(&u4->base);
    free(u4);
    k26astro_vehicle_destroy(v3);
    assert(callback_total_ == 0);

    /* NULL vehicle returns -1. */
    assert(k26astro_vehicle_set_unique_payload(NULL, NULL) == -1);
    assert(k26astro_payload_bind_unique(NULL, NULL) == -1);
}

int main(void)
{
    test_payload_array_slot_();
    test_payload_unique_slot_();
    printf("test_vehicle_payload_slots: OK\n");
    return 0;
}
