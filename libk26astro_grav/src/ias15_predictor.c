/* ias15_predictor.c — Gauss-Radau predictor-corrector iteration.
 *
 * For each of the 7 Radau substeps within a single IAS15 step:
 *   1. Predict body positions at h_n·dt using the current b-
 *      coefficient estimate.
 *   2. Evaluate the full acceleration there.
 *   3. Convert (a_n - a_0) into Newton divided differences g_n.
 *   4. Recompute the b-coefficient row from g via the c matrix.
 *
 * Iterate until |Δb_6| / |a_0| < tol; typically 2-6 iterations
 * for inner-SS scale problems, more during close encounters.
 *
 * The b-coefficients carry from one step to the next as the
 * predictor seed (predict the next step's b's by extrapolation in
 * dt; implemented in ias15.c's driver). */
#include "ias15_internal.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_grav/forces.h"
#include "k26astro_core/pos.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Position predictor at substep h:
 *   x(h) = x_0 + h·dt·v_0
 *        + (h·dt)² · (a_0/2 + h/6·b_0 + h²/12·b_1 + h³/20·b_2 + h⁴/30·b_3
 *                     + h⁵/42·b_4 + h⁶/56·b_5 + h⁷/72·b_6)
 *
 * For numerical stability we factor h·dt out:
 *   x(h) = x_0 + h·dt·(v_0 + h·dt/2·a_0 + h²·dt²/6·b_0 + ...)
 *
 * The exposed accel_total in K26AstroGravState computes accelerations
 * at body[i].pos, so we temporarily overwrite each body's K26AstroPos
 * with the predicted position; restore at end of iteration. */
static void predict_positions_(const K26AstroGravState *state, double dt,
                                double h,
                                const K26V3 *x0, const K26V3 *v0,
                                const K26V3 *a0, K26V3 *const b[7],
                                K26V3 *r_out)
{
    int n = state->n_bodies;
    double h_dt = h * dt;
    double h2 = h * h;
    double h3 = h2 * h, h4 = h2 * h2, h5 = h4 * h, h6 = h3 * h3, h7 = h6 * h;

    for (int i = 0; i < n; i++) {
        K26V3 acc = {
            a0[i].x / 2.0 + h/6.0*b[0][i].x + h2/12.0*b[1][i].x
                + h3/20.0*b[2][i].x + h4/30.0*b[3][i].x
                + h5/42.0*b[4][i].x + h6/56.0*b[5][i].x + h7/72.0*b[6][i].x,
            a0[i].y / 2.0 + h/6.0*b[0][i].y + h2/12.0*b[1][i].y
                + h3/20.0*b[2][i].y + h4/30.0*b[3][i].y
                + h5/42.0*b[4][i].y + h6/56.0*b[5][i].y + h7/72.0*b[6][i].y,
            a0[i].z / 2.0 + h/6.0*b[0][i].z + h2/12.0*b[1][i].z
                + h3/20.0*b[2][i].z + h4/30.0*b[3][i].z
                + h5/42.0*b[4][i].z + h6/56.0*b[5][i].z + h7/72.0*b[6][i].z
        };
        r_out[i].x = x0[i].x + h_dt * (v0[i].x + h_dt * acc.x);
        r_out[i].y = x0[i].y + h_dt * (v0[i].y + h_dt * acc.y);
        r_out[i].z = x0[i].z + h_dt * (v0[i].z + h_dt * acc.z);
    }
}

/* (apply_predicted_pos_ / restore_pos_ were factored out; the
 * predictor inlines the position swap and restore inside the substep
 * loop, since each substep needs its own swap with a different r_pred.) */

/* Update g_k from (a_n - a_0) using Newton's divided differences:
 *
 *   g_0 = (a_1 - a_0) · r_{1,0}
 *   g_k = ((((a_{k+1} - a_0) · r_{k+1,0} - g_0) · r_{k+1,1}
 *           - g_1) · r_{k+1,2} - g_2) · r_{k+1,3} - …
 *
 * This is the standard Newton-form divided-difference accumulator. */
static void update_g_(K26V3 *const g[7], const K26V3 *a_at_sub,
                       const K26V3 *a0, int n,
                       int sub, const double *rr)
{
    /* Compute g[sub-1] from (a_at_sub - a0) and previous g's. */
    /* k = sub - 1 (0-indexed). */
    int k = sub - 1;
    /* Pack-index helper for rr: rr[(i*(i-1))/2 + j] = 1/(h_i - h_j). */
    for (int i = 0; i < n; i++) {
        K26V3 diff = {
            a_at_sub[i].x - a0[i].x,
            a_at_sub[i].y - a0[i].y,
            a_at_sub[i].z - a0[i].z
        };
        /* Apply the divided-difference recursion for this row. */
        K26V3 acc = diff;
        for (int j = 0; j < k; j++) {
            int rr_idx = (sub * (sub - 1)) / 2 + j;
            double r = rr[rr_idx];
            acc.x = acc.x * r - g[j][i].x;
            acc.y = acc.y * r - g[j][i].y;
            acc.z = acc.z * r - g[j][i].z;
        }
        int rr_idx = (sub * (sub - 1)) / 2 + k;
        double r_last = rr[rr_idx];
        g[k][i].x = acc.x * r_last;
        g[k][i].y = acc.y * r_last;
        g[k][i].z = acc.z * r_last;
    }
}

/* Recompute b's from g's: b[k] = g[k] + sum_{j>k} c_{j,k}·g[j].
 *
 * The c matrix is the Newton-form elementary symmetric polynomials
 * (signed). Indexed as c_lt[(k(k-1))/2 + i] with k > i. We need
 * c_{j,k} for j > k. */
static void update_b_from_g_(K26V3 *const b[7], K26V3 *const g[7], int n,
                              const double *c)
{
    /* For each i = 0..6:  b[i] = g[i] + c_{i+1,i}*g[i+1] + c_{i+2,i}*g[i+2] + ... */
    for (int i = 0; i < 7; i++) {
        for (int bi = 0; bi < n; bi++) {
            b[i][bi] = g[i][bi];
            for (int j = i + 1; j < 7; j++) {
                int idx = (j * (j - 1)) / 2 + i;
                b[i][bi].x += c[idx] * g[j][bi].x;
                b[i][bi].y += c[idx] * g[j][bi].y;
                b[i][bi].z += c[idx] * g[j][bi].z;
            }
        }
    }
}

int k26_ias15_pc_iterate(K26AstroGravState *state, double dt,
                          int max_iter, double tol,
                          double *out_max_b6, double *out_max_a0)
{
    K26AstroIAS15Carry *carry = state->ias15_carry;
    int n = state->n_bodies;
    if (!carry || n < 1) return -1;

    k26_ias15_init_matrices();
    const double *c  = k26_ias15_c();
    const double *rr = k26_ias15_rr();

    /* Snapshot initial state. */
    K26V3 *x0 = malloc((size_t)n * sizeof(K26V3));
    K26V3 *v0 = malloc((size_t)n * sizeof(K26V3));
    K26V3 *a0 = malloc((size_t)n * sizeof(K26V3));
    K26AstroPos *pos_saved = malloc((size_t)n * sizeof(K26AstroPos));
    K26V3 *a_n = malloc((size_t)n * sizeof(K26V3));
    K26V3 *r_pred = malloc((size_t)n * sizeof(K26V3));
    /* α-fail-safe scratch: per-body active flag + characteristic
     * length |x_i - COM|. */
    char *body_active = malloc((size_t)n * sizeof(char));
    if (!x0 || !v0 || !a0 || !pos_saved || !a_n || !r_pred || !body_active) {
        free(x0); free(v0); free(a0); free(pos_saved); free(a_n); free(r_pred);
        free(body_active);
        return -1;
    }

    /* Initial position-in-metres + velocity snapshot for the predictor. */
    K26AstroPos origin = state->bodies[0].pos;
    for (int i = 0; i < n; i++) {
        K26V3 d = k26astro_pos_sub(&state->bodies[i].pos, &origin);
        x0[i] = d;
        v0[i] = state->bodies[i].vel;
    }

    /* COM-relative α-fail-safe (Rein & Spiegel 2015 §3.6).
     *
     * Compute system COM in body-0-relative metres, then characterise
     * each body by |x_i - COM| (its distance from the system barycentre).
     * This avoids the "body 0 privileged" trap of using body-0-relative
     * coords directly; body 0 always has x_i=0 and would be fail-safe
     * immune.
     *
     * Mark a body "inactive" (excluded from the max|a_0|/max|b_6|
     * maxima used by the controller's eq. 11) iff its motion during
     * the step is much smaller than its COM-distance, i.e. its
     * contribution to b_6 is FP-noise dominated rather than physical:
     *
     *   |v_i|·|dt| < α · |x_i - COM|     (paper §3.6)
     *
     * Adaptive dt-floor: if NO body has motion exceeding 10·α relative
     * to its COM distance, we're in the chaotic-pocket-escape regime
     * (sub-picosecond dt) where the fail-safe would over-exclude every
     * body. Bypass it then; rely on ias15.c's escape mechanism.
     *
     * α = 1e-8 is the paper-cited value. */
    const double K26_IAS15_ALPHA_FAILSAFE = 1.0e-8;
    {
        double total_m = 0.0;
        K26V3  com_offset = {0.0, 0.0, 0.0};  /* COM - bodies[0].pos in metres */
        for (int i = 0; i < n; i++) total_m += state->bodies[i].gm;
        if (total_m > 0.0) {
            for (int i = 0; i < n; i++) {
                double w = state->bodies[i].gm / total_m;
                com_offset.x += w * x0[i].x;
                com_offset.y += w * x0[i].y;
                com_offset.z += w * x0[i].z;
            }
        }
        double max_motion_ratio = 0.0;
        double r_com_arr[64];   /* small-n optimization, scratch overflow guarded */
        double v_mag_arr[64];
        double *r_com = (n <= 64) ? r_com_arr : malloc((size_t)n * sizeof(double));
        double *v_mag = (n <= 64) ? v_mag_arr : malloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; i++) {
            double rx = x0[i].x - com_offset.x;
            double ry = x0[i].y - com_offset.y;
            double rz = x0[i].z - com_offset.z;
            r_com[i] = sqrt(rx*rx + ry*ry + rz*rz);
            v_mag[i] = sqrt(v0[i].x*v0[i].x + v0[i].y*v0[i].y + v0[i].z*v0[i].z);
            double ratio = (r_com[i] > 0.0)
                ? (v_mag[i] * fabs(dt) / r_com[i])
                : 1.0;
            if (ratio > max_motion_ratio) max_motion_ratio = ratio;
        }
        int failsafe_active =
            (max_motion_ratio >= 10.0 * K26_IAS15_ALPHA_FAILSAFE);
        for (int i = 0; i < n; i++) {
            if (failsafe_active && r_com[i] > 0.0) {
                body_active[i] =
                    (v_mag[i] * fabs(dt)
                     >= K26_IAS15_ALPHA_FAILSAFE * r_com[i]) ? 1 : 0;
            } else {
                body_active[i] = 1;
            }
        }
        if (n > 64) { free(r_com); free(v_mag); }
    }

    /* Initial acceleration (a_0) at the start of the step. */
    k26astro_grav_accel_total(state, a0);
    memcpy(carry->at0, a0, (size_t)n * sizeof(K26V3));

    /* Predictor seed: existing b's from previous step (zero on
     * first call; caller must have initialised carry->b to zero). */

    int iter;
    double resid = 0.0;
    double max_b6 = 0.0;
    double max_a0 = 0.0;
    /* max |a0| diagnostic; computed once on the start-of-step accel.
     * Skip bodies excluded by the α-fail-safe; their accel
     * contribution to the b̃_6 ratio (eq. 9) is FP-noise dominated. */
    for (int i = 0; i < n; i++) {
        if (!body_active[i]) continue;
        double m = sqrt(a0[i].x*a0[i].x + a0[i].y*a0[i].y + a0[i].z*a0[i].z);
        if (m > max_a0) max_a0 = m;
    }
    for (iter = 0; iter < max_iter; iter++) {
        /* Stash b[6] for convergence check. */
        K26V3 *b6_prev = malloc((size_t)n * sizeof(K26V3));
        if (!b6_prev) { iter = -1; break; }
        memcpy(b6_prev, carry->b[6], (size_t)n * sizeof(K26V3));

        /* For each substep h_1..h_7: predict, evaluate, update g. */
        for (int sub = 1; sub <= 7; sub++) {
            double h = k26_ias15_h[sub];
            predict_positions_(state, dt, h,
                                x0, v0, a0, carry->b, r_pred);
            /* Substep positions are referenced relative to origin;
             * convert back to K26AstroPos and stash original. */
            for (int i = 0; i < n; i++) {
                K26V3 abs = {
                    r_pred[i].x + (i == 0 ? 0.0 : 0.0),
                    r_pred[i].y, r_pred[i].z
                };
                (void)abs;
            }
            /* The predictor's r_pred is already in metres relative
             * to origin (since x0 was). Build the absolute K26AstroPos
             * for each body. */
            for (int i = 0; i < n; i++) {
                pos_saved[i] = state->bodies[i].pos;
                state->bodies[i].pos = origin;
                k26astro_pos_add(&state->bodies[i].pos, r_pred[i]);
            }

            k26astro_grav_accel_total(state, a_n);
            update_g_(carry->g, a_n, a0, n, sub, rr);

            /* Restore body positions. */
            for (int i = 0; i < n; i++) state->bodies[i].pos = pos_saved[i];
        }

        /* Rebuild b's from g's. */
        update_b_from_g_(carry->b, carry->g, n, c);

        /* Convergence: max |Δb[6]| / max |b[6]| across all bodies.
         *
         * Per Rein & Spiegel 2015, the iteration must run enough
         * passes to refine the highest-order Newton-form coefficient.
         * Using |Δb_6|/|a_0| as the criterion is wrong for smooth
         * orbits where b_6 is naturally small; the initial b=0 seed
         * has |Δb_6| < |b_6| < tol·|a_0| trivially, terminating
         * before any refinement. The correct relative criterion is
         * |Δb_6|/|b_6|, which only converges once b_6 stops changing.
         *
         * Additionally, enforce a minimum of 2 iterations so the
         * "first pass" always gets refined at least once.
         *
         * max_db drives the iteration convergence test; all bodies
         * count (the iteration must converge on every body's b_6,
         * regardless of fail-safe status).
         *
         * max_b6 feeds the controller's eq. 11 step-size choice and
         * uses only active bodies (the fail-safe excludes bodies
         * whose b_6 is FP-noise dominated). */
        double max_db = 0.0;
        max_b6 = 0.0;
        for (int i = 0; i < n; i++) {
            double db = sqrt((carry->b[6][i].x - b6_prev[i].x)
                            *(carry->b[6][i].x - b6_prev[i].x)
                          + (carry->b[6][i].y - b6_prev[i].y)
                            *(carry->b[6][i].y - b6_prev[i].y)
                          + (carry->b[6][i].z - b6_prev[i].z)
                            *(carry->b[6][i].z - b6_prev[i].z));
            if (db > max_db) max_db = db;
            if (!body_active[i]) continue;
            double b6mag = sqrt(carry->b[6][i].x*carry->b[6][i].x
                              + carry->b[6][i].y*carry->b[6][i].y
                              + carry->b[6][i].z*carry->b[6][i].z);
            if (b6mag > max_b6) max_b6 = b6mag;
        }
        free(b6_prev);
        if (max_b6 > 0.0) {
            resid = max_db / max_b6;
        } else {
            resid = max_db;
        }
        if (iter >= 1 && resid < tol) break;
    }

    /* Diagnostic: store iter count + residual via a side channel so
     * ias15.c can include them in optional logging. */
    extern int    k26_ias15_last_iters;
    extern double k26_ias15_last_resid;
    k26_ias15_last_iters = iter + 1;
    k26_ias15_last_resid = resid;

    if (out_max_b6) *out_max_b6 = max_b6;
    if (out_max_a0) *out_max_a0 = max_a0;

    free(x0); free(v0); free(a0); free(pos_saved); free(a_n); free(r_pred);
    free(body_active);
    return (resid < tol) ? iter + 1 : -1;
}

int    k26_ias15_last_iters = 0;
double k26_ias15_last_resid = 0.0;
