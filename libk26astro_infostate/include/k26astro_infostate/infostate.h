/* k26astro_infostate/infostate.h — retarded-time observer machinery
 * for participants.
 *
 * An infostate is the per-observer record of what the observer
 * "knows" about each tracked target — bounded by the finite speed
 * of light. For two participants A and B separated by range R, A's
 * picture of B is delayed by R/c; tactical decisions made at A's
 * clock-time t must be made against B's state at t - R/c, not at t.
 *
 * The lib provides:
 *
 *   - A per-observer opaque (K26AstroInfostate) that attaches to a
 *     K26AstroVehicle's PAYLOAD_UNIQUE singleton slot (see
 *     k26astro_infostate_attach).
 *
 *   - A per-target history ring buffer that the world step pushes
 *     each target's (t, x, v) state into once per integration tick.
 *
 *   - A retarded-time fixed-point solver
 *     (k26astro_infostate_observe) that iterates
 *
 *         |x_target(t - τ) - x_observer(t)| = c · τ
 *
 *     until |Δτ| falls below tolerance. Newton convergence is
 *     typically reached in 2-3 iterations for engagement-bubble
 *     ranges (mutual ranges 10^4–10^6 m → light-times 30 μs–3 ms);
 *     larger ranges (Earth-Jupiter at 5 AU → 33 minutes) converge
 *     in 3-4 iterations because the target's velocity is small
 *     relative to c.
 *
 * The retarded-time history of each target is reconstructed by
 * linear interpolation between the two history samples bracketing
 * (t - τ). Linear interpolation matches the determinism guarantees
 * of the underlying integrator's per-tick state output without
 * adding higher-order extrapolation error.
 *
 * Determinism: the lib is pure-deterministic — no RNG, no
 * thread-local state. PORTABLE and REFERENCED modes produce
 * bit-identical retarded-time states given identical inputs.
 *
 * Celestial-body targets: the EPHEM modality tag is reserved for
 * targets whose state comes from libk26astro_ephem. The infostate
 * lib does not re-implement ephemeris evaluation — consumers
 * should observe celestial bodies through libk26astro_rt's
 * observer modes (which run their own light-time iteration) and
 * reserve this table for vehicle targets.
 *
 * Scope limitations (stated, not deferred):
 *   - Linear interpolation between history samples; higher-order
 *     interpolation is not provided. Tick rates fast enough that
 *     target velocity changes negligibly between samples are the
 *     expected operating regime.
 *   - History ring buffer has a fixed capacity at construction;
 *     overflow drops the oldest sample. Resolution of retarded
 *     times older than (capacity × mean tick interval) returns
 *     observation.valid = 0.
 *   - Aberration and frame-dragging effects are not modelled;
 *     consumers needing post-Newtonian corrections should pair
 *     this lib with libk26astro_rt's APPARENT observer mode for
 *     celestial-body targets.
 *   - Target pointers are identity-only and tracks are never
 *     retired: a destroyed target's history lingers, and an
 *     allocator that reuses the address for a new vehicle would
 *     silently alias the old history onto it. Do not reuse
 *     per-target history across target lifetimes; a
 *     `_forget_target` surface is the intended future remedy.
 */
#ifndef K26ASTRO_INFOSTATE_INFOSTATE_H
#define K26ASTRO_INFOSTATE_INFOSTATE_H

#include <stdint.h>

#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations — the vehicle and body opaques are defined
 * in their owning libraries. */
struct K26AstroVehicle;

/* ---- Opaque handle -------------------------------------------- */

typedef struct K26AstroInfostate K26AstroInfostate;

/* ---- Observation modality ------------------------------------- *
 *
 * Records which sensor produced the track. Tactical libs use this
 * for downstream attribution (e.g. a radar-derived track has a
 * different uncertainty model than an IR-derived one). EPHEM is
 * reserved for celestial-body targets whose state comes from
 * libk26astro_ephem rather than a sensor reading. */
typedef enum {
    K26ASTRO_INFOSTATE_MODALITY_NONE   = 0,
    K26ASTRO_INFOSTATE_MODALITY_IR     = 1,
    K26ASTRO_INFOSTATE_MODALITY_RADAR  = 2,
    K26ASTRO_INFOSTATE_MODALITY_LIDAR  = 3,
    K26ASTRO_INFOSTATE_MODALITY_EPHEM  = 4
} K26AstroInfostateModality;

/* ---- Observation snapshot ------------------------------------- *
 *
 * Returned by k26astro_infostate_observe and
 * k26astro_infostate_latest. Position and velocity are in the
 * world's inertial reference frame at t_retarded (the emission
 * time at the target, in coordinate time); the observer's clock
 * reads t_observer at the time of reception. */
typedef struct {
    K26AstroEpoch              t_observer;   /* clock at observer */
    K26AstroEpoch              t_retarded;   /* emission time at target */
    K26V3                      position;     /* target r at t_retarded */
    K26V3                      velocity;     /* target v at t_retarded */
    double                     range_m;      /* observer-target distance */
    double                     age_s;        /* t_observer - t_retarded */
    K26AstroInfostateModality  modality;
    int                        iters;        /* fixed-point iterations used */
    int                        valid;        /* 1 if an observation was
                                              * produced, 0 otherwise */
} K26AstroInfostateObservation;

/* ---- Construction --------------------------------------------- */

/* Allocate an infostate bound to an observer vehicle. The observer
 * vehicle must have a bound K26AstroBody (used as the observer's
 * position/velocity source); calling without a bound body returns
 * NULL.
 *
 * history_capacity caps the per-target ring-buffer length. A
 * conservative default of 1024 samples covers engagement bubbles
 * up to ~17 minutes of past state at 1 Hz tick (or 17 seconds at
 * 60 Hz tick). Targets that need deeper history can request more
 * at the cost of linear memory. Minimum value 2 (need at least
 * two samples to interpolate); smaller values clamp to 2.
 *
 * Returns NULL on NULL observer, observer without a bound body, or
 * allocation failure. */
K26AstroInfostate *k26astro_infostate_new(struct K26AstroVehicle *observer,
                                          int history_capacity);

/* Bind the infostate into its observer vehicle's PAYLOAD_UNIQUE
 * singleton slot. Any previous unique payload on the observer is
 * evicted (its notify callback fires). Returns 0 on success, -1 on
 * NULL infostate or missing observer. */
int k26astro_infostate_attach(K26AstroInfostate *s);

/* Free the infostate and unlink from the observer's payload slot.
 * The vehicle slot is nulled in place before the heap allocation is
 * released so the back-reference is closed. Safe on NULL. */
void k26astro_infostate_destroy(K26AstroInfostate *s);

/* Returns the observer vehicle the infostate is bound to. NULL on
 * NULL input. */
struct K26AstroVehicle *
k26astro_infostate_observer(K26AstroInfostate *s);

/* ---- Target history push -------------------------------------- *
 *
 * Push a target's current state into the per-target history ring
 * buffer. Called once per integration tick per (observer, target)
 * pair by the world step loop. Pushing an epoch older than the
 * most-recent entry is a logic error and is silently dropped (the
 * ring buffer expects monotone-increasing epochs); pushing the
 * same epoch twice replaces the previous entry.
 *
 * The first push for a previously-unseen target allocates a new
 * target slot inside the infostate. Slot count is capped at
 * K26ASTRO_INFOSTATE_MAX_TARGETS — see infostate_consts.h.
 * Pushing beyond the cap is a silent no-op.
 *
 * Safe on NULL inputs. */
void k26astro_infostate_target_push(K26AstroInfostate *s,
                                    const struct K26AstroVehicle *target,
                                    K26AstroEpoch t,
                                    K26V3 position_inertial,
                                    K26V3 velocity_inertial);

/* ---- Retarded-time observation -------------------------------- *
 *
 * Resolve the target's state at the observer's clock-time t, using
 * the light-time fixed-point iteration described in the file
 * header. The observer's position is read from its bound body's
 * current state at t (the world is expected to have stepped the
 * observer to t before this call).
 *
 * On convergence, the returned observation has .valid = 1 and the
 * .position/.velocity/.t_retarded fields populated. If the target's
 * history doesn't extend far enough back to cover t - τ (the ring
 * buffer's oldest sample is newer than t - τ), the observation
 * returns .valid = 0 and the caller should treat the target as not
 * yet visible at the requested clock time.
 *
 * The modality parameter is recorded into the observation snapshot
 * unchanged — the lib does not apply per-modality processing in
 * v0.1; that responsibility lies with the sensor-side libk26astro_
 * detect surface that calls this. */
K26AstroInfostateObservation
k26astro_infostate_observe(K26AstroInfostate *s,
                           const struct K26AstroVehicle *target,
                           K26AstroEpoch t,
                           K26AstroInfostateModality modality);

/* ---- Inspection ----------------------------------------------- */

/* Returns the latest pushed sample for `target` as an observation
 * snapshot (modality NONE, age_s=0, range_m=0). Useful for replay
 * and debugging without recomputing the retarded-time iteration.
 * .valid = 0 if target has never been pushed. */
K26AstroInfostateObservation
k26astro_infostate_latest(const K26AstroInfostate *s,
                          const struct K26AstroVehicle *target);

/* Number of history samples currently populated for `target`.
 * Returns 0 if the target has never been pushed or s is NULL. */
int k26astro_infostate_history_length(const K26AstroInfostate *s,
                                      const struct K26AstroVehicle *target);

/* History capacity the infostate was constructed with. */
int k26astro_infostate_history_capacity(const K26AstroInfostate *s);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_INFOSTATE_INFOSTATE_H */
