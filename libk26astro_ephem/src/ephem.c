/* ephem.c — K26AstroEphem query API over a K26AstroSpk + per-body
 * LRU cache for sub-millisecond repeat queries.
 *
 * The cache is keyed by NAIF body id; each body owns a tiny LRU of
 * (record-mid-epoch, segment-pointer, record-pointer) tuples. A
 * query that lands inside the same record as the previous query
 * skips the segment lookup + binary search entirely — important for
 * light-time iteration, which calls back into the ephemeris with an
 * epoch that's only nanoseconds off the previous one.
 */
#include "k26astro_ephem/ephem.h"
#include "k26astro_ephem/spk.h"
#include "k26astro_ephem/cheby.h"

#include "k26astro_core/consts.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* ---- Cache shape ------------------------------------------------ */

#define K26ASTRO_EPHEM_BODIES_MAX        16
#define K26ASTRO_EPHEM_CACHE_PER_BODY     8

typedef struct {
    int                       naif_id;
    /* Most-recently used at slot 0. */
    struct {
        const K26AstroSpkSegment *seg;
        int                       record_idx;
        double                    mid_et;
        double                    half_interval;
    } slots[K26ASTRO_EPHEM_CACHE_PER_BODY];
    int n_slots;
} BodyCache;

struct K26AstroEphem {
    K26AstroSpk *spk;
    BodyCache    cache[K26ASTRO_EPHEM_BODIES_MAX];
    int          n_cache;
};

/* ---- NAIF id ↔ friendly name table ----------------------------- */
typedef struct { int id; const char *name; } NameEntry;
static const NameEntry BODY_NAMES[] = {
    { 0,   "ssb" },        /* solar-system barycentre */
    { 1,   "mercury_b" },  /* mercury barycentre */
    { 2,   "venus_b" },
    { 3,   "earthmoon_b" }, /* Earth-Moon barycentre */
    { 4,   "mars_b" },
    { 5,   "jupiter_b" },
    { 6,   "saturn_b" },
    { 7,   "uranus_b" },
    { 8,   "neptune_b" },
    { 9,   "pluto_b" },
    { 10,  "sun" },
    { 199, "mercury" },
    { 299, "venus" },
    { 301, "moon" },
    { 399, "earth" },
    { 499, "mars" },
    { 599, "jupiter" },
    { 699, "saturn" },
    { 799, "uranus" },
    { 899, "neptune" },
    { 999, "pluto" },
};
static const int N_BODY_NAMES = (int)(sizeof(BODY_NAMES) / sizeof(BODY_NAMES[0]));

int k26astro_ephem_lookup_name(const char *name)
{
    if (!name) return -1;
    for (int i = 0; i < N_BODY_NAMES; i++) {
        if (strcasecmp(BODY_NAMES[i].name, name) == 0) return BODY_NAMES[i].id;
    }
    return -1;
}

const char *k26astro_ephem_id_name(int naif_id)
{
    for (int i = 0; i < N_BODY_NAMES; i++) {
        if (BODY_NAMES[i].id == naif_id) return BODY_NAMES[i].name;
    }
    return NULL;
}

/* ---- Lifecycle -------------------------------------------------- */

K26AstroEphem *k26astro_ephem_load(const char *path)
{
    if (!path) return NULL;
    K26AstroSpk *spk = k26astro_spk_open(path);
    if (!spk) return NULL;
    K26AstroEphem *e = (K26AstroEphem *)calloc(1, sizeof(*e));
    if (!e) { k26astro_spk_close(spk); return NULL; }
    e->spk = spk;
    return e;
}

K26AstroEphem *k26astro_ephem_load_default(void)
{
    return k26astro_ephem_load("/usr/share/k26astro/ephem/de441_inner.spk");
}

K26AstroEphem *k26astro_ephem_load_default_full(void)
{
    /* Probe the full-kernel install path first; fall back to the
     * inner-planet subset so consumers that only need Sun + inner
     * planets still see a valid ephem handle even if the optional
     * full-data APK is not installed. */
    K26AstroEphem *e =
        k26astro_ephem_load("/usr/share/k26astro/ephem/de441_full.spk");
    if (e) return e;
    return k26astro_ephem_load_default();
}

K26AstroEphem *k26astro_ephem_load_default_full_strict(void)
{
    return k26astro_ephem_load("/usr/share/k26astro/ephem/de441_full.spk");
}

int k26astro_ephem_has_outer_planets(const K26AstroEphem *e)
{
    if (!e || !e->spk) return 0;
    int n = k26astro_spk_n_segments(e->spk);
    for (int i = 0; i < n; i++) {
        const K26AstroSpkSegment *seg = k26astro_spk_segment(e->spk, i);
        if (!seg) continue;
        /* Jupiter barycentre (5) or Jupiter (599). The presence of
         * either signals the kernel extends beyond the inner-planet
         * subset; outer planets (Saturn, Uranus, Neptune) come with
         * Jupiter in every DE-series full kernel. */
        if (seg->target_body == 5 || seg->target_body == 599) return 1;
    }
    return 0;
}

void k26astro_ephem_close(K26AstroEphem *e)
{
    if (!e) return;
    if (e->spk) k26astro_spk_close(e->spk);
    free(e);
}

/* ---- Cache lookup ----------------------------------------------- */

static BodyCache *find_or_make_cache_(K26AstroEphem *e, int naif_id)
{
    for (int i = 0; i < e->n_cache; i++) {
        if (e->cache[i].naif_id == naif_id) return &e->cache[i];
    }
    if (e->n_cache >= K26ASTRO_EPHEM_BODIES_MAX) {
        /* Reuse the LRU body slot — for v0.1 we just bump the last
         * one. Reordering on hit lives in the per-record path below,
         * so this is rarely exercised in practice. */
        BodyCache *bc = &e->cache[K26ASTRO_EPHEM_BODIES_MAX - 1];
        memset(bc, 0, sizeof(*bc));
        bc->naif_id = naif_id;
        return bc;
    }
    BodyCache *bc = &e->cache[e->n_cache++];
    memset(bc, 0, sizeof(*bc));
    bc->naif_id = naif_id;
    return bc;
}

/* Returns a (seg, record_idx) pair for the given epoch. Hits the
 * cache on repeat queries inside the same record. */
static int locate_record_(K26AstroEphem *e, int naif_id, double et,
                          const K26AstroSpkSegment **out_seg,
                          int *out_record_idx)
{
    BodyCache *bc = find_or_make_cache_(e, naif_id);
    /* Linear scan over the body's cache (max 8 slots). */
    for (int i = 0; i < bc->n_slots; i++) {
        double mid = bc->slots[i].mid_et;
        double r   = bc->slots[i].half_interval;
        if (et >= mid - r && et <= mid + r) {
            /* Hit. Move this slot to MRU position. */
            if (i != 0) {
                __typeof__(bc->slots[0]) tmp = bc->slots[i];
                for (int j = i; j > 0; j--) bc->slots[j] = bc->slots[j - 1];
                bc->slots[0] = tmp;
            }
            *out_seg        = bc->slots[0].seg;
            *out_record_idx = bc->slots[0].record_idx;
            return 0;
        }
    }
    /* Miss: do the full lookup against the spk index. */
    const K26AstroSpkSegment *s = k26astro_spk_find_segment(e->spk, naif_id, et);
    if (!s) return 1;
    int idx = (int)((et - s->init_et) / s->interval_seconds);
    if (idx < 0) idx = 0;
    if (idx >= s->n_records) idx = s->n_records - 1;
    const double *rec = s->records + (size_t)idx * (size_t)s->record_size_doubles;
    /* Install at MRU slot, evicting the LRU. */
    int n = bc->n_slots;
    if (n < K26ASTRO_EPHEM_CACHE_PER_BODY) bc->n_slots = ++n;
    for (int j = n - 1; j > 0; j--) bc->slots[j] = bc->slots[j - 1];
    bc->slots[0].seg            = s;
    bc->slots[0].record_idx     = idx;
    bc->slots[0].mid_et         = rec[0];
    bc->slots[0].half_interval  = rec[1];
    *out_seg        = s;
    *out_record_idx = idx;
    return 0;
}

/* ---- Epoch → ET (TDB seconds past J2000) ------------------------ */
static double epoch_to_et_tdb_(const K26AstroEpoch *t)
{
    K26AstroEpoch tdb = *t;
    /* Convert to TDB if not already. */
    if (tdb.scale != K26A_TS_TDB) k26astro_epoch_convert(&tdb, K26A_TS_TDB);
    return (double)tdb.days_since_J2000 * 86400.0 + tdb.seconds_of_day;
}

/* ---- SPK observer-chain walker --------------------------------- */
/* DE-series SPK kernels store each body relative to its centre, not
 * relative to the Solar System Barycentre (SSB). For DE441 the inner-
 * planet chains are at most three hops (e.g. Earth 399 → EMB 3 → SSB
 * 0). The single-segment reader was treating every record as already
 * SSB-relative — that masked Earth at EMB-frame ICs into the
 * integrator and produced 1e16 m residuals over 100 years.
 *
 * This walker resolves a NAIF body id to SSB by summing position +
 * velocity at each hop. Composition is strictly linear (vector
 * addition in a fixed inertial frame). Terminates on center_body=0
 * (SSB) or when depth exceeds K26_EPHEM_CHAIN_MAX_DEPTH — six is
 * plenty for any DE-series chain. */
#define K26_EPHEM_CHAIN_MAX_DEPTH 6

static int resolve_to_ssb_(K26AstroEphem *e, int naif_id, double et,
                           K26V3 *out_pos_m, K26V3 *out_vel)
{
    K26V3 pos_m = { 0.0, 0.0, 0.0 };
    K26V3 vel   = { 0.0, 0.0, 0.0 };
    int current = naif_id;
    for (int depth = 0; depth < K26_EPHEM_CHAIN_MAX_DEPTH; depth++) {
        if (current == 0) {
            *out_pos_m = pos_m;
            *out_vel   = vel;
            return 0;
        }
        const K26AstroSpkSegment *seg = NULL;
        int idx = 0;
        if (locate_record_(e, current, et, &seg, &idx) != 0) return 1;
        const double *rec = seg->records
                          + (size_t)idx * (size_t)seg->record_size_doubles;
        double mid    = rec[0];
        double radius = rec[1];
        if (radius == 0.0) return 1;
        double s_norm = (et - mid) / radius;
        if (s_norm >  1.0) s_norm =  1.0;
        if (s_norm < -1.0) s_norm = -1.0;
        int K = seg->coeffs_per_axis;
        const double *cx = rec + 2;
        const double *cy = cx + K;
        const double *cz = cy + K;
        double xk, yk, zk, dx, dy, dz;
        k26astro_cheby_eval_both(cx, K, s_norm, &xk, &dx);
        k26astro_cheby_eval_both(cy, K, s_norm, &yk, &dy);
        k26astro_cheby_eval_both(cz, K, s_norm, &zk, &dz);
        double inv_r = 1.0 / radius;
        /* SPK distances are km; physical velocity (m/s) is
         * (d_km/ds_unit) * (1/radius) * 1000. */
        pos_m.x += xk * 1000.0;
        pos_m.y += yk * 1000.0;
        pos_m.z += zk * 1000.0;
        vel.x   += dx * inv_r * 1000.0;
        vel.y   += dy * inv_r * 1000.0;
        vel.z   += dz * inv_r * 1000.0;
        current = seg->center_body;
    }
    /* Chain too deep — malformed kernel or cycle. */
    return 2;
}

/* ---- Pos / state ------------------------------------------------ */

K26AstroPosTagged k26astro_ephem_body_pos(K26AstroEphem      *e,
                                          int                 naif_id,
                                          const K26AstroEpoch *t)
{
    K26AstroPosTagged out;
    out.p        = k26astro_pos_zero();
    out.frame_id = K26A_FRAME_INVALID;
    if (!e || !t) return out;
    double et = epoch_to_et_tdb_(t);
    K26V3 pos_m, vel;
    if (resolve_to_ssb_(e, naif_id, et, &pos_m, &vel) != 0) return out;
    out.p        = k26astro_pos_from_m(pos_m.x, pos_m.y, pos_m.z);
    out.frame_id = K26A_FRAME_ICRF;
    return out;
}

K26AstroStateXV k26astro_ephem_body_state(K26AstroEphem      *e,
                                          int                 naif_id,
                                          const K26AstroEpoch *t)
{
    K26AstroStateXV out;
    out.pos = k26astro_pos_zero();
    out.vel.x = out.vel.y = out.vel.z = 0.0;
    if (!e || !t) return out;
    double et = epoch_to_et_tdb_(t);
    K26V3 pos_m, vel;
    if (resolve_to_ssb_(e, naif_id, et, &pos_m, &vel) != 0) return out;
    out.pos = k26astro_pos_from_m(pos_m.x, pos_m.y, pos_m.z);
    out.vel = vel;
    return out;
}

/* ---- Light-time correction ------------------------------------ */

K26AstroPosTagged
k26astro_ephem_observe(K26AstroEphem        *e,
                       int                   target_body,
                       const K26AstroPos    *observer_pos,
                       const K26AstroEpoch  *t_obs,
                       int                   max_iter)
{
    K26AstroPosTagged miss;
    miss.p = k26astro_pos_zero();
    miss.frame_id = K26A_FRAME_INVALID;
    if (!e || !observer_pos || !t_obs) return miss;

    K26AstroEpoch t_emit = *t_obs;
    /* Iterate t_emit = t_obs - |r_target(t_emit) - r_observer| / c. */
    K26AstroPosTagged last = k26astro_ephem_body_pos(e, target_body, &t_emit);
    if (last.frame_id == K26A_FRAME_INVALID) return miss;

    if (max_iter < 0) max_iter = 0;
    for (int i = 0; i < max_iter; i++) {
        K26V3 sep = k26astro_pos_sub(&last.p, observer_pos);
        double dist_sq = sep.x * sep.x + sep.y * sep.y + sep.z * sep.z;
        if (dist_sq <= 0.0) break;
        double tau = sqrt(dist_sq) / K26A_C;
        K26AstroEpoch new_t = *t_obs;
        k26astro_epoch_add_seconds(&new_t, -tau);
        /* Convergence test in seconds. */
        if (i > 0) {
            double dt = (double)(new_t.days_since_J2000 - t_emit.days_since_J2000) * 86400.0
                      + (new_t.seconds_of_day - t_emit.seconds_of_day);
            if (dt > -1e-12 && dt < 1e-12) {
                t_emit = new_t;
                break;
            }
        }
        t_emit = new_t;
        last = k26astro_ephem_body_pos(e, target_body, &t_emit);
        if (last.frame_id == K26A_FRAME_INVALID) return miss;
    }

    return last;
}
