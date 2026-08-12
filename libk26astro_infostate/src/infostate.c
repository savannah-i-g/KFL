/* infostate.c — opaque construction, destruction, accessors, and
 * the target-track table machinery.
 *
 * The fixed-point retarded-time solver lives in retarded_time.c;
 * the vehicle-slot composition wiring is in this file. */

#include "k26astro_infostate/infostate.h"
#include "k26astro_infostate/infostate_consts.h"
#include "infostate_internal.h"

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#include <stdlib.h>
#include <string.h>

/* Vehicle-teardown callback carried in the payload base. Blanks the
 * observer back-reference so the retarded-time solver never reads a
 * freed vehicle after the owning vehicle is destroyed. */
static void infostate_on_owner_destroy_(K26AstroPayload *p)
{
    K26AstroInfostate *s = (K26AstroInfostate *)p;
    s->observer = NULL;
}

/* ---- Construction --------------------------------------------- */

K26AstroInfostate *k26astro_infostate_new(struct K26AstroVehicle *observer,
                                          int history_capacity)
{
    if (!observer) return NULL;
    /* The observer must have a bound body — the retarded-time
     * solver reads the observer's position through it. */
    if (!k26astro_vehicle_body(observer)) return NULL;

    if (history_capacity < K26ASTRO_INFOSTATE_MIN_HISTORY_CAPACITY) {
        history_capacity = K26ASTRO_INFOSTATE_MIN_HISTORY_CAPACITY;
    }

    K26AstroInfostate *s = (K26AstroInfostate *)calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->base.kind             = K26ASTRO_DEFENSE_KIND_INFOSTATE;
    s->base.owner_slot       = -1;
    s->base.on_owner_destroy = infostate_on_owner_destroy_;
    s->observer         = observer;
    s->history_capacity = history_capacity;
    s->n_tracks         = 0;
    return s;
}

int k26astro_infostate_attach(K26AstroInfostate *s)
{
    if (!s || !s->observer) return -1;
    return k26astro_payload_bind_unique(&s->base, s->observer);
}

void k26astro_infostate_destroy(K26AstroInfostate *s)
{
    if (!s) return;

    /* Unlink from the observer's PAYLOAD_UNIQUE slot so the
     * back-reference closes before the heap allocation disappears. */
    k26astro_payload_unlink_(&s->base);

    for (int i = 0; i < s->n_tracks; i++) {
        free(s->tracks[i].history);
    }
    free(s);
}

struct K26AstroVehicle *k26astro_infostate_observer(K26AstroInfostate *s)
{
    return s ? s->observer : NULL;
}

int k26astro_infostate_history_capacity(const K26AstroInfostate *s)
{
    return s ? s->history_capacity : 0;
}

/* ---- Target-track table machinery ----------------------------- */

K26AstroInfostateTrack_ *
k26astro_infostate_find_track_(K26AstroInfostate *s,
                               const struct K26AstroVehicle *target)
{
    if (!s || !target) return NULL;
    for (int i = 0; i < s->n_tracks; i++) {
        if (s->tracks[i].live && s->tracks[i].target == target) {
            return &s->tracks[i];
        }
    }
    return NULL;
}

K26AstroInfostateTrack_ *
k26astro_infostate_emplace_track_(K26AstroInfostate *s,
                                  const struct K26AstroVehicle *target)
{
    K26AstroInfostateTrack_ *t = k26astro_infostate_find_track_(s, target);
    if (t) return t;
    if (s->n_tracks >= K26ASTRO_INFOSTATE_MAX_TARGETS) return NULL;

    t = &s->tracks[s->n_tracks];
    t->target   = target;
    t->history  = (K26AstroInfostateHistoryPt_ *)calloc(
                      (size_t)s->history_capacity,
                      sizeof(K26AstroInfostateHistoryPt_));
    if (!t->history) return NULL;
    t->capacity = s->history_capacity;
    t->head     = 0;
    t->count    = 0;
    t->live     = 1;
    s->n_tracks++;
    return t;
}

/* ---- Push --------------------------------------------------- */

void k26astro_infostate_target_push(K26AstroInfostate *s,
                                    const struct K26AstroVehicle *target,
                                    K26AstroEpoch t,
                                    K26V3 position_inertial,
                                    K26V3 velocity_inertial)
{
    if (!s || !target) return;

    K26AstroInfostateTrack_ *trk = k26astro_infostate_emplace_track_(s, target);
    if (!trk) return;  /* table full — silent drop */

    /* Drop pushes older than the most-recent sample. The history
     * buffer expects monotone epochs; an out-of-order push is a
     * caller bug — silently ignored rather than masked. */
    if (trk->count > 0) {
        int latest = (trk->head - 1 + trk->capacity) % trk->capacity;
        double dt = k26astro_epoch_diff_seconds(&t, &trk->history[latest].t);
        if (dt < 0.0) return;
        if (dt == 0.0) {
            /* Replace the latest sample in place. */
            trk->history[latest].position = position_inertial;
            trk->history[latest].velocity = velocity_inertial;
            return;
        }
    }

    trk->history[trk->head].t        = t;
    trk->history[trk->head].position = position_inertial;
    trk->history[trk->head].velocity = velocity_inertial;
    trk->head = (trk->head + 1) % trk->capacity;
    if (trk->count < trk->capacity) trk->count++;
}

/* ---- Inspection ----------------------------------------------- */

K26AstroInfostateHistoryPt_
k26astro_infostate_track_latest_(const K26AstroInfostateTrack_ *track,
                                 int *out_valid)
{
    K26AstroInfostateHistoryPt_ pt;
    memset(&pt, 0, sizeof(pt));
    if (!track || track->count == 0) {
        if (out_valid) *out_valid = 0;
        return pt;
    }
    int idx = (track->head - 1 + track->capacity) % track->capacity;
    pt = track->history[idx];
    if (out_valid) *out_valid = 1;
    return pt;
}

K26AstroInfostateHistoryPt_
k26astro_infostate_track_oldest_(const K26AstroInfostateTrack_ *track,
                                 int *out_valid)
{
    K26AstroInfostateHistoryPt_ pt;
    memset(&pt, 0, sizeof(pt));
    if (!track || track->count == 0) {
        if (out_valid) *out_valid = 0;
        return pt;
    }
    int idx = (track->head - track->count + track->capacity) % track->capacity;
    pt = track->history[idx];
    if (out_valid) *out_valid = 1;
    return pt;
}

int k26astro_infostate_history_length(const K26AstroInfostate *s,
                                      const struct K26AstroVehicle *target)
{
    if (!s || !target) return 0;
    /* find_track_ is non-const; the cast is safe — the lookup
     * doesn't mutate state. */
    K26AstroInfostateTrack_ *t = k26astro_infostate_find_track_(
        (K26AstroInfostate *)s, target);
    return t ? t->count : 0;
}

K26AstroInfostateObservation
k26astro_infostate_latest(const K26AstroInfostate *s,
                          const struct K26AstroVehicle *target)
{
    K26AstroInfostateObservation obs;
    memset(&obs, 0, sizeof(obs));

    if (!s || !target) return obs;
    K26AstroInfostateTrack_ *trk = k26astro_infostate_find_track_(
        (K26AstroInfostate *)s, target);
    if (!trk) return obs;

    int v = 0;
    K26AstroInfostateHistoryPt_ pt =
        k26astro_infostate_track_latest_(trk, &v);
    if (!v) return obs;

    obs.t_observer = pt.t;
    obs.t_retarded = pt.t;
    obs.position   = pt.position;
    obs.velocity   = pt.velocity;
    obs.range_m    = 0.0;
    obs.age_s      = 0.0;
    obs.modality   = K26ASTRO_INFOSTATE_MODALITY_NONE;
    obs.iters      = 0;
    obs.valid      = 1;
    return obs;
}

