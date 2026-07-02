/* payload.h — generic attached-payload base for the vehicle's
 * PAYLOAD composition slots.
 *
 * A payload is any subsystem the vehicle core does not model
 * itself: an instrument package, an experiment, a data recorder,
 * a deployable. External libraries define a concrete struct
 * beginning with a K26AstroPayload field, set the kind tag and
 * (optionally) the on_owner_destroy callback at construction,
 * and bind the payload into either the PAYLOAD array slot
 * (order-preserving, any number) or the PAYLOAD_UNIQUE singleton
 * slot (replace-with-notify, at most one).
 *
 * Dispatch is per-handle via the function pointer carried in the
 * base — no per-kind vehicle hooks exist for payloads, so any
 * number of payload-providing libraries can coexist in one link
 * without strong-definition conflicts.
 *
 * Kind tags are an open namespace owned by the payload-providing
 * library; 0 is reserved as "none". Libraries that need to
 * distinguish their payloads from others' should choose tags with
 * a library-specific high byte and document them.
 *
 * Lifetime contract (same dual-direction protocol as the other
 * composition slots):
 *
 *   - vehicle_destroy: walks the PAYLOAD slots, calls each live
 *     payload's on_owner_destroy(payload) if set, then blanks the
 *     payload's owner back-reference.
 *   - payload-side destroy: the concrete `_destroy` function calls
 *     k26astro_payload_unlink_(p) before freeing its storage —
 *     this nulls the vehicle slot in place so the vehicle's
 *     subsequent destroy walk skips the dead payload.
 *
 * The K26AstroPayload base is the SAME pointer the vehicle slot
 * holds; a concrete struct may be cast to it via its first member.
 */
#ifndef K26ASTRO_VEHICLE_PAYLOAD_H
#define K26ASTRO_VEHICLE_PAYLOAD_H

#include <stdint.h>

#include "k26astro_vehicle/vehicle.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct K26AstroPayload K26AstroPayload;

/* Concrete struct definition. Payload-providing libraries embed
 * this as their first field and downcast back via the
 * (K26AstroPayload *) pointer. Owner back-reference + slot index
 * are maintained by the bind / unlink helpers below.
 *
 * owner_slot: index >= 0 while bound into the PAYLOAD array;
 * K26ASTRO_PAYLOAD_SLOT_UNIQUE while occupying the singleton;
 * -1 while unattached. */
#define K26ASTRO_PAYLOAD_SLOT_UNIQUE (-2)

struct K26AstroPayload {
    uint32_t                  kind;      /* open tag; 0 = none */
    struct K26AstroVehicle   *owner;
    int                       owner_slot;
    uint64_t                  generation;

    /* Called by the vehicle's PAYLOAD-slot teardown when the
     * vehicle is destroyed (or a unique payload is replaced) while
     * this payload is still attached. May be NULL; the owner
     * back-reference is blanked by the vehicle either way. */
    void (*on_owner_destroy)(K26AstroPayload *p);
};

/* Read the kind tag. Safe on NULL (returns 0). */
uint32_t k26astro_payload_kind(const K26AstroPayload *p);

/* Bind a payload into the vehicle's PAYLOAD array slot. Sets the
 * back-reference + records the slot index. Returns 0 on success,
 * -1 on NULL inputs or slot-array allocation failure. */
int k26astro_payload_bind(K26AstroPayload *p, struct K26AstroVehicle *v);

/* Bind a payload into the vehicle's PAYLOAD_UNIQUE singleton slot.
 * Any previous occupant is notified (its on_owner_destroy fires,
 * its owner back-reference is blanked) before replacement.
 * Returns 0 on success, -1 on NULL inputs. */
int k26astro_payload_bind_unique(K26AstroPayload *p, struct K26AstroVehicle *v);

/* Unlink a payload from its owning vehicle. Called from the
 * concrete `_destroy` function before freeing. Nulls the vehicle
 * slot in place so the vehicle's subsequent destroy walk skips
 * the dead payload. Safe on NULL or unattached. */
void k26astro_payload_unlink_(K26AstroPayload *p);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_PAYLOAD_H */
