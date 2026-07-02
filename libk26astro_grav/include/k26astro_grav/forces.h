/* k26astro_grav/forces.h — N-body force evaluation.
 *
 * The view abstraction lets a future Barnes-Hut + SOA pass swap the
 * inner-loop implementation without touching integrator drivers.
 *
 * Memory layout is currently AOS (an array of K26AstroBody). The
 * view carries a const pointer + count; the integrator pipeline
 * never reads body->pos directly from inside a force loop; it
 * goes through the view. A future SOA migration would change the
 * view to { const double *xs, *ys, *zs, *ms; int n; } with loop
 * bodies updated in force_direct.c only.
 *
 * Force accumulation in portable mode switches to Neumaier
 * compensated summation (libk26astro_core/sum.h) when
 * n_bodies ≥ K26ASTRO_COMPENSATED_SUM_THRESHOLD (= 16). Below the
 * threshold, plain summation is bit-identical anyway.
 */
#ifndef K26ASTRO_GRAV_FORCES_H
#define K26ASTRO_GRAV_FORCES_H

#include "k26astro_body/body.h"
#include "k26astro_grav/grav.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct K26AstroGravView {
    const K26AstroBody *bodies;
    int n;
} K26AstroGravView;

/* Build a view over a body array. */
static inline K26AstroGravView k26astro_grav_view(const K26AstroBody *bodies,
                                                   int n)
{
    K26AstroGravView v = { bodies, n };
    return v;
}

/* Direct O(N²) N-body acceleration. Writes per-body accel into
 * `accel_out` (caller-allocated, length n). Acceleration is the
 * gradient of Newtonian point-mass potential summed over all other
 * bodies. */
void k26astro_grav_force_direct(const K26AstroGravView *view,
                                K26V3 *accel_out);

/* Direct N² with Plummer softening: r² → r² + ε². Use for clusters
 * with close encounters; ε is the per-body softening radius. */
void k26astro_grav_force_direct_softened(const K26AstroGravView *view,
                                          double softening,
                                          K26V3 *accel_out);

/* Total acceleration: direct N² (softened if state->softening > 0)
 * + all enabled perturbations + user-registered perturbs.
 *
 * The integrator drivers all call this single entry point so any
 * future perturbation registration shows up automatically.
 *
 * If state->mercurius is non-NULL, the direct-N² pair loop applies
 * the per-pair K weight specified in the context (see
 * mercurius.h). Perturbations are NOT weighted (J2/SRP/GR are
 * additive corrections that don't participate in the MERCURIUS
 * splitting). */
void k26astro_grav_accel_total(const K26AstroGravState *state,
                                K26V3 *accel_out);

/* ---- MERCURIUS pair-by-pair force decomposition ----------------- */

/* Per-pair K weight for the Rein-Tamayo 2019 force split.
 * `k_weight` is K(y_ij) per the quintic smoothstep: 1 = entirely
 * near (IAS15), 0 = entirely far (WH), 0..1 in the transition. */
typedef struct {
    int    i, j;
    double k_weight;
} K26AstroPairWeight;

typedef enum {
    /* Outer integrator: applies (1-K)·F_ij per encounter pair plus
     * full F_ij for all non-encounter pairs. */
    K26ASTRO_MERCURIUS_FAR  = 1,
    /* Inner integrator: applies K·F_ij per encounter pair plus
     * zero for all non-encounter pairs. */
    K26ASTRO_MERCURIUS_NEAR = 2
} K26AstroMercuriusMode;

/* MERCURIUS orchestration context. Caller sets state->mercurius to
 * point at one of these for the duration of a step, then clears.
 * `pair_weights` is a caller-owned array; the integrator does not
 * mutate it. */
struct K26AstroMercuriusContext {
    K26AstroMercuriusMode      mode;
    const K26AstroPairWeight  *pair_weights;
    int                        n_pair_weights;
};

/* Pair-by-pair weighted direct N². Caller-side decomposition
 * primitive (does not consult state->mercurius). Emits both a_far
 * and a_near arrays in one pass over the pair list, useful for
 * tests + the MERCURIUS orchestrator's continuity diagnostics.
 *
 * Semantics:
 *   - For each (i, j) with i < j:
 *       F_ij = full Newtonian pair force
 *       K_ij = pair_weights[k].k_weight if (i,j) present, else 0
 *       a_far[i]  += (1 - K_ij) · F_ij / m_i;  a_far[j]  -= (1 - K_ij) · F_ij / m_j
 *       a_near[i] +=      K_ij · F_ij / m_i;   a_near[j] -=      K_ij · F_ij / m_j
 *   - For non-encounter pairs (K_ij = 0): full force on a_far, zero
 *     on a_near.
 *
 * `out_far` and/or `out_near` may be NULL to skip that side. */
void k26astro_grav_force_direct_weighted(
        const K26AstroGravView    *view,
        const K26AstroPairWeight  *pair_weights, int n_pair_weights,
        K26V3 *out_far, K26V3 *out_near);

/* ---- Event-time root-finding registry --------------------------- *
 *
 * Discrete events (stage drop, parachute deploy, mode transition,
 * engine ignition / cutoff at a precise epoch, plane crossing, etc.)
 * fire mid-step inside the integrator. Each registered event
 * exposes:
 *
 *   - `predicate(state, t, ctx)` : 0 if the event has not fired by
 *     `t`, 1 if it has. `t` is the current simulation time in
 *     seconds since J2000 (consistent across substeps). The
 *     predicate must be monotonic across the candidate interval —
 *     at most one sign change per step.
 *
 *   - `handler(state, t, ctx)` : called at the bisected event epoch
 *     with the state already advanced to `t`. The handler MAY
 *     mutate the state (mass step, inertia step, body removal,
 *     velocity discontinuity); it is the only place where such
 *     mutations are sanctioned mid-step.
 *
 * k26astro_grav_step bisects within the integrator step to locate
 * t_event ∈ [t_start, t_start + dt] (within state->event_tol_s,
 * default 1e-6 s) and advances the remainder of the step with the
 * post-event state. The integrator's full order is preserved across
 * the discontinuity because the substep length is reduced
 * gracefully on either side of the event boundary.
 *
 * Determinism: bisection is bit-stable on IEEE-754. Multiple events
 * at the same epoch resolve in registration order. */
typedef int  (*K26AstroGravEventPredicateFn)(
                  const K26AstroGravState *state, double t, void *ctx);
typedef int  (*K26AstroGravEventHandlerFn)  (
                  K26AstroGravState *state,       double t, void *ctx);

typedef struct {
    K26AstroGravEventPredicateFn predicate;
    K26AstroGravEventHandlerFn   handler;
    void                        *ctx;
} K26AstroGravEvent;

/* Register an event. Lifetime of `event.ctx` is the caller's
 * responsibility. Returns K26ASTRO_E_OK on success,
 * K26ASTRO_E_NULL / K26ASTRO_E_ALLOC on failure. */
int  k26astro_grav_register_event(K26AstroGravState *state,
                                   K26AstroGravEvent event);

/* Drop all registered events and free the registry buffer. No-op
 * if no events were ever registered. */
void k26astro_grav_clear_events  (K26AstroGravState *state);

/* Override the event-time bisection tolerance (seconds). Default
 * 1e-6 s set at k26astro_grav_state_init. Values <= 0 are rejected
 * and the existing tolerance retained. */
int  k26astro_grav_set_event_tol_s(K26AstroGravState *state,
                                    double tol_s);

#ifdef __cplusplus
}
#endif

#endif
