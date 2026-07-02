/* event_registry.c — K26AstroGravEvent registration + lifetime.
 *
 * Mirrors the perturbation-list pattern in perturb_registry.c: lazy
 * heap allocation, doubling realloc from an initial capacity of 4.
 *
 * The full struct K26AstroEventList layout lives in grav_state.c so
 * the destroy path can free the inner buffer without an additional
 * accessor. This file owns the public register / clear / set-tol
 * surface. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"

#include <stdlib.h>

/* Layout shared with grav_state.c. */
struct K26AstroEventList {
    K26AstroGravEvent *events;
    int                count;
    int                capacity;
};

int k26astro_grav_register_event(K26AstroGravState *state,
                                 K26AstroGravEvent event)
{
    if (!state || !event.predicate || !event.handler) return K26ASTRO_E_NULL;
    if (!state->events) {
        state->events = calloc(1, sizeof(struct K26AstroEventList));
        if (!state->events) return K26ASTRO_E_ALLOC;
    }
    struct K26AstroEventList *p = state->events;
    if (p->count >= p->capacity) {
        int new_cap = p->capacity ? p->capacity * 2 : 4;
        K26AstroGravEvent *new_events = realloc(p->events,
            (size_t)new_cap * sizeof(K26AstroGravEvent));
        if (!new_events) return K26ASTRO_E_ALLOC;
        p->events   = new_events;
        p->capacity = new_cap;
    }
    p->events[p->count++] = event;
    return K26ASTRO_E_OK;
}

void k26astro_grav_clear_events(K26AstroGravState *state)
{
    if (!state || !state->events) return;
    free(state->events->events);
    free(state->events);
    state->events = NULL;
}

int k26astro_grav_set_event_tol_s(K26AstroGravState *state, double tol_s)
{
    if (!state) return K26ASTRO_E_NULL;
    if (!(tol_s > 0.0)) return K26ASTRO_E_BAD_INPUT;
    state->event_tol_s = tol_s;
    return K26ASTRO_E_OK;
}
