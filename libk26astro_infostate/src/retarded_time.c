/* retarded_time.c — light-time fixed-point solver and ring-buffer
 * temporal interpolation.
 *
 * Solves |x_target(t - τ) - x_observer(t)| = c · τ for τ, where the
 * target's state at (t - τ) is reconstructed by linear
 * interpolation between the two history samples that bracket it.
 * Convergence criterion: |τ_{k+1} - τ_k| < K26ASTRO_INFOSTATE_FIXED_POINT_TOL_S.
 *
 * Iteration cap K26ASTRO_INFOSTATE_FIXED_POINT_MAX_ITERS guards
 * pathological geometries (target approaching at ~c) but is not
 * reached in any physically realistic engagement-bubble regime —
 * 2-3 iterations is typical.
 *
 * Linear interpolation rationale: the per-tick state output of
 * libk26astro_grav's integrator is the canonical sampled
 * representation of the target's trajectory in REFERENCED-mode
 * replay. Higher-order interpolation between samples would
 * introduce error not represented in the op-log; linear interp
 * keeps the retarded-time observation a deterministic function of
 * the recorded samples plus the call-site arguments. */

#include "k26astro_infostate/infostate.h"
#include "k26astro_infostate/infostate_consts.h"
#include "infostate_internal.h"

#include "k26astro_core/consts.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26m3d.h"

#include <math.h>
#include <string.h>

/* ---- Track-at-clock-time helper ------------------------------- *
 *
 * Returns the linearly-interpolated history sample at clock time t.
 * Sets out_valid = 1 if t falls inside [oldest, latest], 0 if t is
 * outside the retained-history window. Backward extrapolation is
 * rejected (the buffer's oldest sample is the earliest knowable
 * state); forward extrapolation is also rejected (the integrator
 * has not advanced the target past `latest` yet).
 *
 * The search is linear; the ring is small (typically <=1024) and
 * the bracketing sample is usually near the latest end (small
 * light-times), so a backwards-from-head linear scan beats binary
 * search for the common case. */
K26AstroInfostateHistoryPt_
k26astro_infostate_track_at_(const K26AstroInfostateTrack_ *track,
                             K26AstroEpoch t,
                             int *out_valid)
{
    K26AstroInfostateHistoryPt_ result;
    memset(&result, 0, sizeof(result));
    if (out_valid) *out_valid = 0;

    if (!track || track->count < 1) return result;

    /* If only one sample, exact match is the only valid answer. */
    if (track->count == 1) {
        int idx = (track->head - 1 + track->capacity) % track->capacity;
        double dt_s = k26astro_epoch_diff_seconds(&t, &track->history[idx].t);
        if (dt_s == 0.0) {
            if (out_valid) *out_valid = 1;
            return track->history[idx];
        }
        return result;
    }

    /* Latest sample bound. */
    int latest = (track->head - 1 + track->capacity) % track->capacity;
    double dt_latest = k26astro_epoch_diff_seconds(&t, &track->history[latest].t);
    if (dt_latest > 0.0) return result;             /* forward extrap rejected */
    if (dt_latest == 0.0) {
        if (out_valid) *out_valid = 1;
        return track->history[latest];
    }

    /* Oldest sample bound. */
    int oldest = (track->head - track->count + track->capacity) % track->capacity;
    double dt_oldest = k26astro_epoch_diff_seconds(&t, &track->history[oldest].t);
    if (dt_oldest < 0.0) return result;             /* backward extrap rejected */
    if (dt_oldest == 0.0) {
        if (out_valid) *out_valid = 1;
        return track->history[oldest];
    }

    /* Bracket: walk back from latest until we find the two
     * neighbours straddling t. */
    int hi = latest;
    int lo = (latest - 1 + track->capacity) % track->capacity;
    int steps_back = 1;
    while (steps_back < track->count) {
        double dt_lo = k26astro_epoch_diff_seconds(&t, &track->history[lo].t);
        if (dt_lo >= 0.0) break;
        hi = lo;
        lo = (lo - 1 + track->capacity) % track->capacity;
        steps_back++;
    }

    /* Linear interp between lo and hi at t. */
    double dt_lo = k26astro_epoch_diff_seconds(&t, &track->history[lo].t);
    double dt_hi = k26astro_epoch_diff_seconds(&track->history[hi].t,
                                               &track->history[lo].t);
    if (dt_hi <= 0.0) return result;                /* degenerate — should not happen */
    double alpha = dt_lo / dt_hi;
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    K26V3 r = k26m3d_v3_add(
        k26m3d_v3_scale(track->history[lo].position, 1.0 - alpha),
        k26m3d_v3_scale(track->history[hi].position, alpha));
    K26V3 v = k26m3d_v3_add(
        k26m3d_v3_scale(track->history[lo].velocity, 1.0 - alpha),
        k26m3d_v3_scale(track->history[hi].velocity, alpha));

    result.t        = t;
    result.position = r;
    result.velocity = v;
    if (out_valid) *out_valid = 1;
    return result;
}

/* ---- Observer position read ----------------------------------- */

static K26V3 observer_position_at_(const K26AstroVehicle *observer)
{
    K26V3 zero = { 0.0, 0.0, 0.0 };
    if (!observer) return zero;
    K26AstroBody *b = k26astro_vehicle_body((K26AstroVehicle *)observer);
    if (!b) return zero;

    /* The body's pos is a K26AstroPos (sector-grid). Subtract the
     * canonical zero to get an inertial-frame K26V3 in metres. */
    K26AstroPos origin = k26astro_pos_zero();
    return k26astro_pos_sub(&b->pos, &origin);
}

/* ---- The fixed-point retarded-time iteration ------------------ */

K26AstroInfostateObservation
k26astro_infostate_observe(K26AstroInfostate *s,
                           const struct K26AstroVehicle *target,
                           K26AstroEpoch t,
                           K26AstroInfostateModality modality)
{
    K26AstroInfostateObservation obs;
    memset(&obs, 0, sizeof(obs));
    obs.t_observer = t;
    obs.modality   = modality;

    if (!s || !target) return obs;
    if (!s->observer) return obs;  /* stale: vehicle-destroy hook nulled observer */

    K26AstroInfostateTrack_ *trk = k26astro_infostate_find_track_(s, target);
    if (!trk || trk->count < 1) return obs;

    const K26V3 x_obs = observer_position_at_(s->observer);
    const double c = K26A_C;

    /* Seed τ from the latest history sample's distance. */
    int v0 = 0;
    K26AstroInfostateHistoryPt_ pt0 =
        k26astro_infostate_track_latest_(trk, &v0);
    if (!v0) return obs;

    double r0 = k26m3d_v3_len(k26m3d_v3_sub(pt0.position, x_obs));
    double tau = r0 / c;

    int iters = 0;
    double tau_prev = 0.0;
    K26AstroInfostateHistoryPt_ pt_retarded = pt0;

    while (iters < K26ASTRO_INFOSTATE_FIXED_POINT_MAX_ITERS) {
        K26AstroEpoch t_ret = t;
        k26astro_epoch_add_seconds(&t_ret, -tau);

        int valid = 0;
        K26AstroInfostateHistoryPt_ pt =
            k26astro_infostate_track_at_(trk, t_ret, &valid);
        if (!valid) {
            /* History does not extend back to t - τ. The caller
             * sees .valid = 0 and treats the target as not visible
             * yet. */
            return obs;
        }
        pt_retarded = pt;

        double range = k26m3d_v3_len(k26m3d_v3_sub(pt.position, x_obs));
        double tau_new = range / c;
        double dtau = tau_new - tau;
        if (dtau < 0.0) dtau = -dtau;

        iters++;
        tau_prev = tau;
        tau      = tau_new;

        if (dtau < K26ASTRO_INFOSTATE_FIXED_POINT_TOL_S) break;
        (void)tau_prev;  /* reserved for damped iteration if needed */
    }

    K26AstroEpoch t_ret = t;
    k26astro_epoch_add_seconds(&t_ret, -tau);

    obs.t_retarded = t_ret;
    obs.position   = pt_retarded.position;
    obs.velocity   = pt_retarded.velocity;
    obs.range_m    = k26m3d_v3_len(k26m3d_v3_sub(pt_retarded.position, x_obs));
    obs.age_s      = tau;
    obs.iters      = iters;
    obs.valid      = 1;
    return obs;
}
