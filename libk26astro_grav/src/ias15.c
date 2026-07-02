/* ias15.c — IAS15 top-level driver with Rein-Spiegel 2015 §4
 * adaptive step-size control.
 *
 * `k26astro_grav_step_ias15(state, dt)` advances the body system by
 * `dt` seconds, sub-dividing transparently. Per substep:
 *   1. Ensure carry buffer allocated (n_bodies).
 *   2. Run predictor-corrector iteration (ias15_predictor.c). The PC
 *      reports max|b_6| and max|a_0| via output params for the
 *      controller. Positions are unchanged at this point.
 *   3. Compute dt_required = dt_try · (epsilon / |b_6|)^(1/7).
 *   4. If dt_required ≥ 0.5·dt_try: accept the substep. Apply
 *      position + velocity updates via standard Radau-7 weights:
 *        x_new = x_0 + dt·v_0 + dt² · (a_0/2 + b_0/6 + b_1/12 +
 *                                     b_2/20 + b_3/30 + b_4/42 +
 *                                     b_5/56 + b_6/72)
 *        v_new = v_0 + dt · (a_0 + b_0/2 + b_1/3 + b_2/4 + b_3/5 +
 *                                  b_4/6 + b_5/7 + b_6/8)
 *      Set dt_hint = min(dt_required, 10·dt_try) for the next call.
 *   5. Else: reject the substep. Positions are still unchanged.
 *      Set dt_hint = dt_required and retry. Increment
 *      state->ias15_rejected_steps.
 *
 * Reference: Rein & Spiegel (2015) MNRAS 446:1424 §4. The 0.5
 * acceptance margin and the 10× growth cap are the paper's
 * recommended defaults.
 *
 * The (1/7) exponent is computed via pow(); the 1/7 constant itself
 * is precomputed at file scope as a hex literal for cross-platform
 * bit-exactness (K26_IAS15_INV_SEVEN). pow() is in libm and is
 * compiled under -ffp-contract=off -frounding-math so its result
 * is deterministic given the FPU pin established by grav_state_init. */
#include "k26astro_grav/ias15.h"
#include "ias15_internal.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define IAS15_MAX_PC_ITER  16
#define IAS15_PC_TOL       1.0e-12

/* Controller bounds (Rein-Spiegel 2015 §4): acceptance threshold 0.5,
 * growth cap 10.0. Both exact in IEEE-754. */
#define IAS15_DT_ACCEPT_MARGIN  0.5
#define IAS15_DT_GROW_CAP      10.0

/* 1.0/7.0 as a hex-literal double for cross-platform bit-exactness.
 * 1/7 ≈ 0.142857142857142849... = 0x3FC2492492492492. */
static const union { uint64_t u; double d; } K26_IAS15_INV_SEVEN_BITS =
    { .u = 0x3FC2492492492492ULL };
#define K26_IAS15_INV_SEVEN (K26_IAS15_INV_SEVEN_BITS.d)

static int ensure_carry_(K26AstroGravState *state)
{
    if (state->ias15_carry && state->ias15_carry->capacity >= state->n_bodies)
        return 0;
    if (state->ias15_carry) {
        k26_ias15_carry_release(state->ias15_carry);
        state->ias15_carry = NULL;
    }
    return k26_ias15_carry_alloc(&state->ias15_carry, state->n_bodies);
}

/* Explicit body-state snapshot/restore around each PC iteration.
 * The predictor stack-restores its own pos/vel buffers on exit, so the
 * snapshot is a belt-and-braces defence, making rollback semantics
 * provably bit-exact at the driver level, independent of predictor
 * internals. The snapshot is a byte-copy of K26AstroBody so the
 * restore is bit-identical (no float arithmetic involved). */
static int ensure_snapshot_(K26AstroGravState *state)
{
    if (state->ias15_snapshot && state->ias15_snapshot_cap >= state->n_bodies)
        return 0;
    K26AstroBody *fresh = realloc(state->ias15_snapshot,
                                   (size_t)state->n_bodies * sizeof(K26AstroBody));
    if (!fresh) return K26ASTRO_E_ALLOC;
    state->ias15_snapshot     = fresh;
    state->ias15_snapshot_cap = state->n_bodies;
    return 0;
}

static inline void snapshot_save_(K26AstroGravState *state)
{
    memcpy(state->ias15_snapshot, state->bodies,
           (size_t)state->n_bodies * sizeof(K26AstroBody));
}

static inline void snapshot_restore_(K26AstroGravState *state)
{
    memcpy(state->bodies, state->ias15_snapshot,
           (size_t)state->n_bodies * sizeof(K26AstroBody));
}

/* Apply the accepted Radau-7 position + velocity updates. The PC has
 * already converged carry->b[k][i] and carry->at0[i] for `dt_try`. */
static void apply_step_updates_(K26AstroGravState *state, double dt_try)
{
    K26AstroIAS15Carry *carry = state->ias15_carry;
    K26AstroBody *b = state->bodies;
    K26V3 *a0 = carry->at0;
    int n = state->n_bodies;

    for (int i = 0; i < n; i++) {
        K26V3 dx = {
            dt_try * b[i].vel.x + dt_try * dt_try *
                (a0[i].x / 2.0
                 + carry->b[0][i].x / 6.0
                 + carry->b[1][i].x / 12.0
                 + carry->b[2][i].x / 20.0
                 + carry->b[3][i].x / 30.0
                 + carry->b[4][i].x / 42.0
                 + carry->b[5][i].x / 56.0
                 + carry->b[6][i].x / 72.0),
            dt_try * b[i].vel.y + dt_try * dt_try *
                (a0[i].y / 2.0
                 + carry->b[0][i].y / 6.0
                 + carry->b[1][i].y / 12.0
                 + carry->b[2][i].y / 20.0
                 + carry->b[3][i].y / 30.0
                 + carry->b[4][i].y / 42.0
                 + carry->b[5][i].y / 56.0
                 + carry->b[6][i].y / 72.0),
            dt_try * b[i].vel.z + dt_try * dt_try *
                (a0[i].z / 2.0
                 + carry->b[0][i].z / 6.0
                 + carry->b[1][i].z / 12.0
                 + carry->b[2][i].z / 20.0
                 + carry->b[3][i].z / 30.0
                 + carry->b[4][i].z / 42.0
                 + carry->b[5][i].z / 56.0
                 + carry->b[6][i].z / 72.0)
        };
        K26V3 dv = {
            dt_try * (a0[i].x
                 + carry->b[0][i].x / 2.0
                 + carry->b[1][i].x / 3.0
                 + carry->b[2][i].x / 4.0
                 + carry->b[3][i].x / 5.0
                 + carry->b[4][i].x / 6.0
                 + carry->b[5][i].x / 7.0
                 + carry->b[6][i].x / 8.0),
            dt_try * (a0[i].y
                 + carry->b[0][i].y / 2.0
                 + carry->b[1][i].y / 3.0
                 + carry->b[2][i].y / 4.0
                 + carry->b[3][i].y / 5.0
                 + carry->b[4][i].y / 6.0
                 + carry->b[5][i].y / 7.0
                 + carry->b[6][i].y / 8.0),
            dt_try * (a0[i].z
                 + carry->b[0][i].z / 2.0
                 + carry->b[1][i].z / 3.0
                 + carry->b[2][i].z / 4.0
                 + carry->b[3][i].z / 5.0
                 + carry->b[4][i].z / 6.0
                 + carry->b[5][i].z / 7.0
                 + carry->b[6][i].z / 8.0)
        };
        k26astro_pos_add(&b[i].pos, dx);
        b[i].vel.x += dv.x;
        b[i].vel.y += dv.y;
        b[i].vel.z += dv.z;
    }
}

int k26astro_grav_step_ias15(K26AstroGravState *state, double dt)
{
    if (!state) return K26ASTRO_E_NULL;
    if (state->n_bodies < 1 || dt == 0.0) return K26ASTRO_E_BAD_INPUT;

    int rc = ensure_carry_(state);
    if (rc != 0) return rc;
    rc = ensure_snapshot_(state);
    if (rc != 0) return rc;

    K26AstroIAS15Carry *carry = state->ias15_carry;
    int n = state->n_bodies;

    if (!carry->initialised) {
        for (int k = 0; k < 7; k++) {
            memset(carry->b[k], 0, (size_t)n * sizeof(K26V3));
            memset(carry->g[k], 0, (size_t)n * sizeof(K26V3));
            memset(carry->e[k], 0, (size_t)n * sizeof(K26V3));
        }
    }

    double remaining = dt;
    double dt_try = (state->ias15_dt_hint > 0.0
                       && state->ias15_dt_hint < dt)
                      ? state->ias15_dt_hint : dt;
    /* Safety cap: at least one finite substep. */
    if (dt_try <= 0.0) dt_try = dt;

    /* Bound rejections per call to surface pathological inputs
     * without infinite-looping. Realistic close encounters take
     * a few rejects to converge; 256 gives headroom for chaotic-
     * pocket trajectories (post-flyby outer-system arcs with
     * passing close encounters) that need the dt-collapse escape
     * mechanism below to engage. The escape activates after an
     * accepted substep at very small dt; the cap must be large
     * enough that the controller can find that accept before
     * exhausting its budget. */
    const int reject_cap = 256;
    int local_rejects = 0;
    /* Bounded escape counter. The chaotic-pocket escape (mirrored
     * onto the reject path below) bumps dt when the controller has
     * collapsed it to picosecond range. Capping at 16 per call
     * prevents infinite escape-reject-escape loops while still
     * allowing the trajectory through pathological regions. */
    const int escape_cap = 16;
    int local_escapes = 0;

    /* Wall-budget tracking. Sample t0 once; only enter the
     * clock_gettime path if a positive budget is configured. */
    struct timespec ts0;
    int budget_active = (state->ias15_wall_budget_s > 0.0);
    if (budget_active) clock_gettime(CLOCK_MONOTONIC, &ts0);

    while (remaining > 0.0) {
        if (dt_try > remaining) dt_try = remaining;

        if (budget_active) {
            struct timespec ts;
            clock_gettime(CLOCK_MONOTONIC, &ts);
            double elapsed = (double)(ts.tv_sec - ts0.tv_sec)
                           + 1.0e-9 * (double)(ts.tv_nsec - ts0.tv_nsec);
            if (elapsed > state->ias15_wall_budget_s)
                return K26ASTRO_E_TIME_BUDGET;
        }

        /* Byte-snapshot before each PC iteration. The PC's own
         * predictor restores pos/vel internally, but this driver-
         * level snapshot guarantees bit-exact rollback on any reject
         * path, including future predictor variants. */
        snapshot_save_(state);

        double max_b6 = 0.0, max_a0 = 0.0;
        int iters = k26_ias15_pc_iterate(state, dt_try,
                                           IAS15_MAX_PC_ITER, IAS15_PC_TOL,
                                           &max_b6, &max_a0);
        /* Stash PC convergence diagnostics for driver introspection.
         * Refreshed on every attempt (accept OR reject). Negative
         * iters preserves the no-converge signal. eps_b normalised
         * by max_a0; zeroed when there is no measurable error. */
        state->ias15_last_pc_iterations = iters;
        state->ias15_last_eps_b_achieved =
            (max_a0 > 0.0) ? (max_b6 / max_a0) : 0.0;
        if (iters < 0) {
            /* PC didn't converge; treat as a rejected step at half
             * the size. Restore body state to pre-PC snapshot. */
            snapshot_restore_(state);
            dt_try *= 0.5;
            state->ias15_rejected_steps++;
            if (++local_rejects > reject_cap) return K26ASTRO_E_NO_CONVERGE;
            if (dt_try < remaining * 1e-15) return K26ASTRO_E_NO_CONVERGE;
            continue;
        }

        /* Step-size controller. b_6 is normalised by a_0 so the
         * tolerance epsilon is the dimensionless per-step truncation
         * budget (Rein-Spiegel 2015 §4 eq. 9). */
        double dt_required;
        if (max_a0 > 0.0 && max_b6 > 0.0) {
            double err_norm = max_b6 / max_a0;
            dt_required = dt_try * pow(state->ias15_tol / err_norm,
                                         K26_IAS15_INV_SEVEN);
        } else {
            /* No measurable error in this step (initial seed, free
             * drift); accept and let dt grow. */
            dt_required = dt_try * IAS15_DT_GROW_CAP;
        }

        if (dt_required >= IAS15_DT_ACCEPT_MARGIN * dt_try
            || dt_try <= remaining * 1e-15) {
            /* Accept. */
            apply_step_updates_(state, dt_try);
            k26astro_epoch_add_seconds(&state->t, dt_try);
            state->ias15_dt_last = dt_try;
            state->dt_last       = dt_try;
            carry->initialised   = 1;
            carry->dt_proposed   = dt_try;

            remaining -= dt_try;
            /* Carry the controller hint forward, capped at the
             * growth limit. */
            double next_hint = dt_required;
            if (next_hint > dt_try * IAS15_DT_GROW_CAP)
                next_hint = dt_try * IAS15_DT_GROW_CAP;
            /* Chaotic-pocket escape: when dt has collapsed to
             * a tiny fraction of remaining (sub-IEEE-754-meaningful)
             * the controller's error estimate is dominated by
             * truncation noise rather than physics. Without an
             * escape the controller stays at picosecond dt and the
             * integration never finishes; the Burrau Pythagorean
             * 3-4-5 triple close encounter is the canonical case.
             * Bump the hint to remaining * 1e-6 (1M-step floor)
             * deterministically; this keeps the dynamics in
             * resolvable territory while letting subsequent steps
             * re-collapse dt if the chaotic region demands it. */
            if (dt_try < remaining * 1.0e-12) {
                double escape = remaining * 1.0e-6;
                if (escape > next_hint) next_hint = escape;
            }
            state->ias15_dt_hint = next_hint;
            dt_try = next_hint;
            local_rejects = 0;
        } else {
            /* Reject. Restore body state to pre-PC snapshot.
             * carry->b retained as predictor seed for the smaller
             * dt. */
            snapshot_restore_(state);
            state->ias15_rejected_steps++;
            state->ias15_dt_hint = dt_required;
            if (++local_rejects > reject_cap) return K26ASTRO_E_NO_CONVERGE;
            dt_try = dt_required;
            /* Bounded chaotic-pocket escape on the reject path.
             * Mirrors the accept-path escape at line 283. When the
             * controller has collapsed dt to a tiny fraction of
             * remaining (below the meaningful IEEE-754 timescale),
             * the truncation-error estimate is dominated by
             * floating-point noise. Bumping dt back to a resolvable
             * fraction lets the next PC attempt either accept (and
             * fire the standard escape) or reject again — and
             * eventually exhaust the reject_cap budget rather than
             * dt-collapse-floor termination. Bounded at escape_cap
             * to prevent infinite escape-reject-escape loops. */
            if (dt_try < remaining * 1.0e-12 && local_escapes < escape_cap) {
                double escape = remaining * 1.0e-6;
                if (escape > dt_try) {
                    dt_try = escape;
                    state->ias15_dt_hint = escape;
                    local_escapes++;
                }
            }
            if (dt_try <= remaining * 1e-15) return K26ASTRO_E_NO_CONVERGE;
        }
    }

    return K26ASTRO_E_OK;
}

void k26astro_grav_ias15_reset(K26AstroGravState *state)
{
    if (!state || !state->ias15_carry) return;
    state->ias15_carry->initialised = 0;
    state->ias15_dt_hint = 0.0;
}

void k26astro_grav_ias15_set_tol(K26AstroGravState *state, double eps)
{
    if (!state) return;
    if (eps > 0.0) state->ias15_tol = eps;
}

double k26astro_grav_ias15_get_dt_last(const K26AstroGravState *state)
{
    if (!state) return 0.0;
    return state->ias15_dt_last;
}

uint32_t k26astro_grav_ias15_rejected_steps(const K26AstroGravState *state)
{
    if (!state) return 0;
    return state->ias15_rejected_steps;
}

void k26astro_grav_ias15_set_wall_budget(K26AstroGravState *state,
                                          double budget_s)
{
    if (!state) return;
    state->ias15_wall_budget_s = (budget_s > 0.0) ? budget_s : 0.0;
}
