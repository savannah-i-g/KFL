/* grav_step_internal.h — private dispatch helper.
 *
 * The public k26astro_grav_step routes through the event-time
 * wrapper when state->events has registered events; otherwise it
 * dispatches the inner integrator directly. The wrapper in
 * advance_with_events.c uses this helper as its sub-step entry
 * point. */
#ifndef K26ASTRO_GRAV_STEP_INTERNAL_H
#define K26ASTRO_GRAV_STEP_INTERNAL_H

#include "k26astro_grav/grav.h"

/* Dispatch the configured integrator for a single step of duration
 * `dt`. Mirrors the original k26astro_grav_step body before the
 * event-wrapper split. Does NOT consult the event registry. */
int k26astro_grav_step_inner_dispatch(K26AstroGravState *state, double dt);

/* Step with event-time root-finding. Called by k26astro_grav_step
 * when the event registry is non-empty. */
int k26astro_grav_step_with_events(K26AstroGravState *state, double dt);

#endif /* K26ASTRO_GRAV_STEP_INTERNAL_H */
