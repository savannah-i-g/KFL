/* k26astro_ephem/ephem.h — high-level ephemeris query API.
 *
 * Wraps an opaque K26AstroEphem (which owns an mmap'd K26AstroSpk
 * plus a per-body LRU cache for repeat queries) and exposes a clean
 * K26AstroPos / K26AstroEpoch interface.
 *
 * All positions are returned in the J2000 / ICRF frame in metres —
 * the SPK file's native distance unit is kilometres; we scale
 * inside the query path so callers never see km.
 *
 * Epochs are converted from the caller's K26AstroEpoch (any scale)
 * to TDB internally; SPK ephemeris time is TDB. */
#ifndef K26ASTRO_EPHEM_EPHEM_H
#define K26ASTRO_EPHEM_EPHEM_H

#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/frame.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct K26AstroEphem K26AstroEphem;

/* ---- Load / close ----------------------------------------------- */
K26AstroEphem *k26astro_ephem_load(const char *path);

/* Load the default DE441 inner-planet kernel installed by the K26
 * ISO (/usr/share/k26astro/ephem/de441_inner.spk). Returns NULL if
 * the kernel isn't installed. */
K26AstroEphem *k26astro_ephem_load_default(void);

/* Load the full DE441 kernel installed by the optional
 * libk26astro_ephem-data-full APK
 * (/usr/share/k26astro/ephem/de441_full.spk, ~3.3 GB covering
 * 1550-2650 with Sun + all planets + Moon). Falls back to the
 * inner-planet kernel if the full kernel is absent. Consumers
 * that require Jovian / outer-planet coverage should call
 * `_load_default_full_strict` instead, then test the handle for
 * NULL — or call this and gate on `_has_outer_planets`. Returns
 * NULL only if neither kernel is installed. */
K26AstroEphem *k26astro_ephem_load_default_full(void);

/* Load the full DE441 kernel without the inner-only fallback.
 * Returns NULL if the full kernel is not installed. Consumers
 * needing Jovian or outer-planet coverage should prefer this
 * variant over `_load_default_full` so an inner-only fallback
 * does not silently mask missing coverage. */
K26AstroEphem *k26astro_ephem_load_default_full_strict(void);

/* Returns 1 if the loaded kernel covers Jupiter (NAIF 5 or 599)
 * via any of its segments, else 0. Scans the kernel's segment
 * catalogue; no per-epoch evaluation. Use this to gate consumers
 * that require outer-planet coverage when the kernel was loaded
 * via the tolerant `_load_default_full` variant. */
int k26astro_ephem_has_outer_planets(const K26AstroEphem *e);

void           k26astro_ephem_close(K26AstroEphem *e);

/* ---- Body identity --------------------------------------------- */
/* NAIF id lookup. Names are case-insensitive: "earth", "EARTH",
 * "Earth" all match. Returns -1 if the name isn't known. */
int         k26astro_ephem_lookup_name(const char *name);
const char *k26astro_ephem_id_name(int naif_id);

/* ---- Position query -------------------------------------------- */
typedef struct {
    K26AstroPos pos;
    K26V3       vel;   /* metres per second */
} K26AstroStateXV;

/* Position of `naif_id` at epoch `t`. Returns a position tagged
 * ICRF. On query failure (no segment covers the epoch, or unknown
 * body), the returned position has frame_id == K26A_FRAME_INVALID
 * — callers should test this rather than trust the components. */
K26AstroPosTagged k26astro_ephem_body_pos(K26AstroEphem      *e,
                                          int                 naif_id,
                                          const K26AstroEpoch *t);

/* Position + velocity. Velocity comes from the analytic Chebyshev
 * derivative, so it's bit-equivalent across platforms. */
K26AstroStateXV k26astro_ephem_body_state(K26AstroEphem      *e,
                                          int                 naif_id,
                                          const K26AstroEpoch *t);

/* ---- Light-time-corrected observation -------------------------- */
/* Apparent position of `target_body` as seen from `observer_pos`
 * at `t_obs`. Solves |r_target(t_obs - τ) - r_observer(t_obs)| = c τ
 * by fixed-point iteration. Converges in 2-3 iterations for inner
 * solar-system distances. Returns the apparent position tagged ICRF.
 *
 * max_iter = 0 → no correction (position at t_obs exactly,
 *   useful as the "geometric" reference).
 * max_iter = 1 → one Newton step (rough light-time correction).
 * Typical: max_iter = 4 — converges to sub-nanosecond on inner SS. */
K26AstroPosTagged
k26astro_ephem_observe(K26AstroEphem        *e,
                       int                   target_body,
                       const K26AstroPos    *observer_pos,
                       const K26AstroEpoch  *t_obs,
                       int                   max_iter);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_EPHEM_EPHEM_H */
