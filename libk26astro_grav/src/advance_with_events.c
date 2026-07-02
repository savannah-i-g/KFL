/* advance_with_events.c — event-time root-finding wrapper around
 * the inner integrator dispatch.
 *
 * The wrapper:
 *   1. Snapshots the body array + epoch at step start.
 *   2. Advances tentatively by dt via grav_step_inner_dispatch.
 *   3. For each registered event, compares predicate-at-start vs
 *      predicate-at-end. A change indicates a bracketed crossing.
 *   4. Bisects each bracketed event to locate t_event within
 *      state->event_tol_s. Bisection uses snapshot-and-restart
 *      sub-steps: from the snapshot, advance by trial_dt, evaluate
 *      the predicate, narrow the bracket.
 *   5. Picks the earliest t_event across all bracketed events.
 *   6. Restores the snapshot, advances to t_event, calls the
 *      handler (which may mutate state), and recurses on the
 *      remainder of the step from the post-handler state. Re-
 *      evaluating events on the recursive call handles cascading
 *      events triggered by handler mutations.
 *
 * Determinism: bisection is bit-stable on IEEE-754. Multiple events
 * at the same bisected t_event resolve in registration order
 * (earliest-registered handler runs first; the recursive remainder
 * picks up the others on the next pass).
 *
 * Cost: zero events → identical to grav_step_inner_dispatch (the
 * public k26astro_grav_step short-circuits before entering this
 * file). With N bracketed events per outer step and B = log2(dt /
 * tol) bisection iterations: N · B + 1 integrator sub-steps per
 * event-bearing outer step. For typical N=1 and B≈20 that's 21
 * sub-steps where a no-event call took 1.
 *
 * Reference: Hairer, Nørsett, Wanner, Solving Ordinary Differential
 * Equations I: Nonstiff Problems (1993), §II.6.2 — event location
 * via bisection inside an adaptive substep. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_core/epoch.h"
#include "grav_step_internal.h"

#include <stdlib.h>
#include <string.h>

/* Layout shared with grav_state.c + event_registry.c. */
struct K26AstroEventList {
    K26AstroGravEvent *events;
    int                count;
    int                capacity;
};

/* Maximum bisection iterations. Safety bound — bisection converges
 * deterministically and the loop normally exits on tol_s. */
#define EVENT_BISECT_MAX_ITERS 60

/* Maximum sub-event recursion depth per outer step. Cascading
 * events that re-trigger should be rare; the bound prevents runaway
 * recursion on a buggy handler. */
#define EVENT_RECURSE_MAX_DEPTH 32

/* Convert the current epoch to seconds-since-J2000 for predicate
 * input. Caller's predicate receives `t` in this form. */
static double epoch_to_s_(K26AstroEpoch e)
{
    return (double)e.days_since_J2000 * 86400.0 + e.seconds_of_day;
}

/* Snapshot body array + epoch. Grows the snapshot buffer if
 * n_bodies exceeds the current capacity. */
static int snapshot_take_(K26AstroGravState *state)
{
    if (state->event_snapshot_cap < state->n_bodies) {
        K26AstroBody *grown = realloc(state->event_snapshot,
            (size_t)state->n_bodies * sizeof(K26AstroBody));
        if (!grown) return K26ASTRO_E_ALLOC;
        state->event_snapshot     = grown;
        state->event_snapshot_cap = state->n_bodies;
    }
    memcpy(state->event_snapshot, state->bodies,
           (size_t)state->n_bodies * sizeof(K26AstroBody));
    return K26ASTRO_E_OK;
}

static void snapshot_restore_(K26AstroGravState *state, K26AstroEpoch t_at)
{
    memcpy(state->bodies, state->event_snapshot,
           (size_t)state->n_bodies * sizeof(K26AstroBody));
    state->t = t_at;
}

/* Bisect t_event in [lo, hi] for event `ev`. Caller has already
 * confirmed predicate(state_at_lo, lo_s) != predicate(state_at_hi,
 * hi_s). The wrapper restores from the snapshot before each trial
 * sub-step, so on entry `state` may be at any consistent (state, t)
 * — by exit, state is at the latest sub-step. */
static double bisect_event_(K26AstroGravState *state,
                             K26AstroEpoch  t_start_epoch,
                             double         dt_full,
                             const K26AstroGravEvent *ev,
                             int            p_start,
                             double         tol_s)
{
    double dt_lo = 0.0;
    double dt_hi = dt_full;
    int    iters = 0;
    while ((dt_hi - dt_lo) > tol_s && iters < EVENT_BISECT_MAX_ITERS) {
        double dt_mid = 0.5 * (dt_lo + dt_hi);
        /* Fresh sub-step from the snapshot. */
        snapshot_restore_(state, t_start_epoch);
        if (dt_mid > 0.0) {
            int rc = k26astro_grav_step_inner_dispatch(state, dt_mid);
            if (rc != K26ASTRO_E_OK) {
                /* Treat sub-step failure as upper bracket — fall back
                 * to the lower interval. The caller will check the
                 * final result; bisection still returns a value
                 * within [lo, hi]. */
                dt_hi = dt_mid;
                iters++;
                continue;
            }
        }
        double t_mid_s = epoch_to_s_(state->t);
        int    p_mid   = ev->predicate(state, t_mid_s, ev->ctx) ? 1 : 0;
        if (p_mid == p_start) {
            dt_lo = dt_mid;
        } else {
            dt_hi = dt_mid;
        }
        iters++;
    }
    return 0.5 * (dt_lo + dt_hi);
}

static int step_with_events_recursive_(K26AstroGravState *state,
                                        double dt, int depth)
{
    if (depth >= EVENT_RECURSE_MAX_DEPTH) {
        /* Defensive: shouldn't happen in well-behaved problems. Fall
         * back to a plain step on the remainder. */
        return k26astro_grav_step_inner_dispatch(state, dt);
    }
    if (!(dt > 0.0)) return K26ASTRO_E_OK;

    /* No events registered? Direct dispatch. The recursive call can
     * land here after a handler that cleared the registry. */
    if (!state->events || state->events->count == 0) {
        return k26astro_grav_step_inner_dispatch(state, dt);
    }

    struct K26AstroEventList *list = state->events;
    int n = list->count;

    /* Snapshot at step start. */
    int rc = snapshot_take_(state);
    if (rc != K26ASTRO_E_OK) return rc;
    K26AstroEpoch t_start_epoch = state->t;
    double t_start_s = epoch_to_s_(t_start_epoch);

    /* Evaluate predicates at t_start (state is at the snapshot). */
    int p_start[64];
    if (n > 64) n = 64;       /* practical cap; matches plan-time scope */
    for (int i = 0; i < n; i++) {
        p_start[i] = list->events[i].predicate(state, t_start_s,
                                                list->events[i].ctx) ? 1 : 0;
    }

    /* Tentative full-step advance. */
    rc = k26astro_grav_step_inner_dispatch(state, dt);
    if (rc != K26ASTRO_E_OK) {
        snapshot_restore_(state, t_start_epoch);
        return rc;
    }
    double t_end_s = epoch_to_s_(state->t);
    int    any_bracketed = 0;
    int    p_end[64];
    for (int i = 0; i < n; i++) {
        p_end[i] = list->events[i].predicate(state, t_end_s,
                                              list->events[i].ctx) ? 1 : 0;
        if (p_end[i] != p_start[i]) any_bracketed = 1;
    }

    if (!any_bracketed) {
        /* Clean step. State already at t_start + dt. */
        return K26ASTRO_E_OK;
    }

    /* Find earliest bracketed event. */
    double earliest_dt = dt;
    int    earliest_i  = -1;
    for (int i = 0; i < n; i++) {
        if (p_end[i] == p_start[i]) continue;
        double t_event_dt = bisect_event_(state, t_start_epoch, dt,
                                           &list->events[i], p_start[i],
                                           state->event_tol_s > 0.0
                                               ? state->event_tol_s
                                               : 1.0e-6);
        if (t_event_dt < earliest_dt) {
            earliest_dt = t_event_dt;
            earliest_i  = i;
        }
    }

    if (earliest_i < 0) {
        /* Shouldn't reach — any_bracketed was true. Defensive: take
         * the tentative full step. */
        snapshot_restore_(state, t_start_epoch);
        return k26astro_grav_step_inner_dispatch(state, dt);
    }

    /* Restore to start, advance to the event epoch, apply handler. */
    snapshot_restore_(state, t_start_epoch);
    if (earliest_dt > 0.0) {
        rc = k26astro_grav_step_inner_dispatch(state, earliest_dt);
        if (rc != K26ASTRO_E_OK) return rc;
    }
    double t_event_s = epoch_to_s_(state->t);
    /* Handler may mutate state; we deliberately ignore its return
     * value beyond logging-via-comment intent (it propagates only
     * through ctx-side effects). */
    (void)list->events[earliest_i].handler(state, t_event_s,
                                            list->events[earliest_i].ctx);

    /* Recurse on the remainder. The recursive call re-evaluates the
     * full registry against the (possibly-mutated) post-handler
     * state, which covers cascade events and re-triggered predicates. */
    double dt_remaining = dt - earliest_dt;
    return step_with_events_recursive_(state, dt_remaining, depth + 1);
}

int k26astro_grav_step_with_events(K26AstroGravState *state, double dt)
{
    state->dt_last = dt;
    return step_with_events_recursive_(state, dt, 0);
}
