/* stage_event.c — vehicle stage-event timeline with registration
 * into the libk26astro_grav event-time root-finder.
 *
 * The timeline is sorted ascending by epoch. Insertions binary-
 * search the insertion point (linear walk; n_events is small) and
 * memmove the tail.
 *
 * Epoch-to-seconds conversion mirrors libk26astro_grav's
 * advance_with_events.c exactly:
 *
 *   t_s = (double)e.days_since_J2000 * 86400.0 + e.seconds_of_day
 *
 * The same arithmetic both sides yields bit-exact predicate
 * comparison against the integrator's working time, so the
 * predicate is well-behaved on the candidate substep interval.
 *
 * Handler semantics:
 *
 *   1. The gravity bisector restores the pre-substep snapshot.
 *   2. It advances the inner integrator to t_event.
 *   3. The handler runs with the integrator state at t_event;
 *      mutations to bodies[].mass / .gm survive into the next
 *      outer substep's fresh snapshot.
 *
 * The handler therefore:
 *   - decrements vehicle basic_mass,
 *   - propagates the new mass to the bound body via
 *     k26astro_body_set_mass,
 *   - updates the Ext attitude state's inertia diagonal via
 *     k26astro_attitude_update_inertia (off-diagonals preserved). */

#include "k26astro_vehicle/stage_event.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26astro_body/attitude.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_core/epoch.h"
#include "vehicle_internal.h"

#include <stdlib.h>
#include <string.h>

/* ---- Epoch <-> seconds helper -------------------------------- */

static double epoch_to_s_(K26AstroEpoch e)
{
    return (double)e.days_since_J2000 * 86400.0 + e.seconds_of_day;
}

static int epoch_less_(K26AstroEpoch a, K26AstroEpoch b)
{
    return epoch_to_s_(a) < epoch_to_s_(b);
}

/* ---- Predicate + handler ------------------------------------- */

/* Predicate is one-shot: once the handler has set fired=1 the
 * predicate stays on the post-event side. Without this, the
 * bisector's slight undershoot (handler runs at t_event = scheduled
 * - epsilon) lets the recursive substep see another sign change
 * across [scheduled - epsilon, scheduled] and the same event would
 * fire repeatedly until basic_mass underflows. */
static int predicate_(const K26AstroGravState *state, double t, void *ctx)
{
    (void)state;
    const K26AstroVehicleEventCtx_ *ec = (const K26AstroVehicleEventCtx_ *)ctx;
    if (ec->event->fired) return 1;
    return (t >= ec->scheduled_t_s) ? 1 : 0;
}

static int handler_(K26AstroGravState *state, double t, void *ctx)
{
    (void)state;
    (void)t;
    K26AstroVehicleEventCtx_   *ec = (K26AstroVehicleEventCtx_ *)ctx;
    K26AstroVehicle            *v  = ec->vehicle;
    K26AstroVehicleStageEvent_ *ev = ec->event;

    /* Mass step. */
    double new_mass = v->basic_mass_kg - ev->mass_drop_kg;
    if (new_mass < 0.0) new_mass = 0.0;
    v->basic_mass_kg = new_mass;

    /* Propagate to the bound body so the gravity integrator sees
     * the new GM at the next substep. */
    if (v->body) {
        k26astro_body_set_mass(v->body, new_mass);
    }

    /* Inertia diagonal update — off-diagonals preserved from the
     * current tensor. */
    K26M3 I = v->attitude_ext.inertia;
    I.m[0][0] = ev->ixx_after;
    I.m[1][1] = ev->iyy_after;
    I.m[2][2] = ev->izz_after;
    k26astro_attitude_update_inertia(&v->attitude_ext, I);

    ev->fired = 1;
    return K26ASTRO_E_OK;
}

/* ---- Schedule ------------------------------------------------- */

int k26astro_vehicle_schedule_stage_event(K26AstroVehicle *v,
                                          K26AstroEpoch epoch,
                                          double mass_drop_kg,
                                          double ixx_after,
                                          double iyy_after,
                                          double izz_after)
{
    if (!v) return -1;
    if (v->n_events >= v->cap_events) {
        int new_cap = v->cap_events ? (v->cap_events * 2) : 4;
        K26AstroVehicleStageEvent_ *grown = (K26AstroVehicleStageEvent_ *)realloc(
            v->events,
            (size_t)new_cap * sizeof(K26AstroVehicleStageEvent_));
        if (!grown) return -1;
        v->events     = grown;
        v->cap_events = new_cap;
    }
    /* Insertion-point search — linear; events lists are short. */
    int i = 0;
    while (i < v->n_events && epoch_less_(v->events[i].epoch, epoch)) i++;
    if (i < v->n_events) {
        memmove(&v->events[i + 1], &v->events[i],
                (size_t)(v->n_events - i) * sizeof(K26AstroVehicleStageEvent_));
    }
    v->events[i].epoch        = epoch;
    v->events[i].mass_drop_kg = (mass_drop_kg < 0.0) ? 0.0 : mass_drop_kg;
    v->events[i].ixx_after    = ixx_after;
    v->events[i].iyy_after    = iyy_after;
    v->events[i].izz_after    = izz_after;
    v->events[i].fired        = 0;
    v->n_events++;
    return i;
}

int k26astro_vehicle_event_count(const K26AstroVehicle *v)
{
    return v ? v->n_events : 0;
}

K26AstroEpoch k26astro_vehicle_event_epoch_at(const K26AstroVehicle *v, int idx)
{
    if (!v || idx < 0 || idx >= v->n_events) {
        return k26astro_epoch_j2000_tt();
    }
    return v->events[idx].epoch;
}

/* ---- Flush to a gravity-state registry ------------------------ */

int k26astro_vehicle_register_stage_events_with(K26AstroVehicle *v,
                                                K26AstroGravState *grav)
{
    if (!v || !grav) return K26ASTRO_E_NULL;

    /* Grow the parallel ctx pool to fit the new events. */
    int needed = v->n_ctxs + v->n_events;
    if (needed > v->cap_ctxs) {
        int new_cap = v->cap_ctxs ? v->cap_ctxs : 4;
        while (new_cap < needed) new_cap *= 2;
        K26AstroVehicleEventCtx_ **grown = (K26AstroVehicleEventCtx_ **)realloc(
            v->event_ctxs,
            (size_t)new_cap * sizeof(K26AstroVehicleEventCtx_ *));
        if (!grown) return K26ASTRO_E_ALLOC;
        v->event_ctxs = grown;
        v->cap_ctxs   = new_cap;
    }

    for (int i = 0; i < v->n_events; i++) {
        K26AstroVehicleEventCtx_ *ec = (K26AstroVehicleEventCtx_ *)malloc(sizeof(*ec));
        if (!ec) return K26ASTRO_E_ALLOC;
        ec->vehicle       = v;
        ec->event         = &v->events[i];
        ec->scheduled_t_s = epoch_to_s_(v->events[i].epoch);

        /* Track for cleanup on destroy. */
        v->event_ctxs[v->n_ctxs++] = ec;

        K26AstroGravEvent ge = { predicate_, handler_, ec };
        int rc = k26astro_grav_register_event(grav, ge);
        if (rc != K26ASTRO_E_OK) return rc;
    }
    return K26ASTRO_E_OK;
}
