/* infostate_internal.h — private struct layout shared across the
 * lib's .c files. Not installed; not part of the public ABI. */
#ifndef K26ASTRO_INFOSTATE_INTERNAL_H
#define K26ASTRO_INFOSTATE_INTERNAL_H

#include <stdint.h>

#include "k26astro_infostate/infostate.h"
#include "k26astro_infostate/infostate_consts.h"
#include "k26astro_core/epoch.h"
#include "k26astro_vehicle/payload.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Per-target history entry --------------------------------- */
typedef struct {
    K26AstroEpoch  t;
    K26V3          position;
    K26V3          velocity;
} K26AstroInfostateHistoryPt_;

/* ---- Per-target track record ---------------------------------- *
 *
 * The history buffer is a ring; .head indexes the next write slot
 * and .count tracks the populated length up to capacity. The
 * latest sample lives at (head - 1) mod capacity; the oldest lives
 * at (head - count) mod capacity. */
typedef struct {
    const struct K26AstroVehicle *target;     /* identity + non-owning */
    K26AstroInfostateHistoryPt_  *history;    /* capacity-sized ring */
    int                           capacity;
    int                           head;
    int                           count;
    int                           live;       /* 0 if slot freed */
} K26AstroInfostateTrack_;

/* ---- Infostate struct ----------------------------------------- */
struct K26AstroInfostate {
    K26AstroPayload          base;            /* vehicle PAYLOAD_UNIQUE slot; first member */
    struct K26AstroVehicle  *observer;        /* non-owning */
    int                      history_capacity;

    /* Per-target table; sparse — count is the number of unique
     * targets ever pushed (capped at K26ASTRO_INFOSTATE_MAX_TARGETS).
     * Target lookup is linear scan; the table is small enough that
     * binary indexing would not pay back. */
    K26AstroInfostateTrack_  tracks[K26ASTRO_INFOSTATE_MAX_TARGETS];
    int                      n_tracks;
};

/* ---- Internal helpers (cross-file) ---------------------------- */

/* Locates the track record for `target`, or NULL if not yet
 * present. */
K26AstroInfostateTrack_ *
k26astro_infostate_find_track_(K26AstroInfostate *s,
                               const struct K26AstroVehicle *target);

/* Locates or appends; returns NULL if the target table is full. */
K26AstroInfostateTrack_ *
k26astro_infostate_emplace_track_(K26AstroInfostate *s,
                                  const struct K26AstroVehicle *target);

/* Returns the latest sample for `track` (most-recently-pushed),
 * with `out_valid = 1`; returns zero-initialised + `out_valid = 0`
 * if the track is empty. */
K26AstroInfostateHistoryPt_
k26astro_infostate_track_latest_(const K26AstroInfostateTrack_ *track,
                                 int *out_valid);

/* Returns the oldest sample retained in the ring. */
K26AstroInfostateHistoryPt_
k26astro_infostate_track_oldest_(const K26AstroInfostateTrack_ *track,
                                 int *out_valid);

/* Linear-interpolates the track to clock-time t. Returns
 * `out_valid = 1` and populates position/velocity if t falls
 * inside the retained-history window; returns `out_valid = 0`
 * otherwise. Extrapolation forward of the latest sample is
 * explicitly rejected (caller should push then re-query). */
K26AstroInfostateHistoryPt_
k26astro_infostate_track_at_(const K26AstroInfostateTrack_ *track,
                             K26AstroEpoch t,
                             int *out_valid);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_INFOSTATE_INTERNAL_H */
