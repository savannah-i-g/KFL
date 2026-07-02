/* payload.c — bind / unlink helpers for the generic PAYLOAD
 * composition slots. Dispatch on vehicle teardown is per-handle
 * (K26AstroPayload.on_owner_destroy) and lives in composition.c's
 * default notify shim; no per-kind vehicle hooks exist for
 * payloads. */

#include "k26astro_vehicle/payload.h"
#include "k26astro_vehicle/vehicle.h"

#include <stddef.h>

uint32_t k26astro_payload_kind(const K26AstroPayload *p)
{
    return p ? p->kind : 0u;
}

int k26astro_payload_bind(K26AstroPayload *p, struct K26AstroVehicle *v)
{
    if (!p || !v) return -1;
    int slot = k26astro_vehicle_add_payload(v, p);
    if (slot < 0) return -1;
    p->owner      = v;
    p->owner_slot = slot;
    return 0;
}

int k26astro_payload_bind_unique(K26AstroPayload *p, struct K26AstroVehicle *v)
{
    if (!p || !v) return -1;
    if (k26astro_vehicle_set_unique_payload(v, p) < 0) return -1;
    p->owner      = v;
    p->owner_slot = K26ASTRO_PAYLOAD_SLOT_UNIQUE;
    return 0;
}

void k26astro_payload_unlink_(K26AstroPayload *p)
{
    if (!p || !p->owner) return;
    k26astro_vehicle_invalidate_slot(p->owner,
        p->owner_slot == K26ASTRO_PAYLOAD_SLOT_UNIQUE
            ? K26ASTRO_VEHICLE_SLOT_PAYLOAD_UNIQUE
            : K26ASTRO_VEHICLE_SLOT_PAYLOAD,
        p);
    p->owner      = NULL;
    p->owner_slot = -1;
}
