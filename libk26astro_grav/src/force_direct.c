/* force_direct.c — direct O(N²) gravitational force evaluation.
 *
 * The hottest inner loop in the integrator pipeline. Sole owner of
 * the per-pair acceleration math; all integrator drivers go through
 * this file to compute accelerations. A future SOA migration would
 * replace this file alone; the function signatures taking a view
 * keep integrator code untouched.
 *
 * Determinism:
 *   - Pair ordering: iterate (i, j) with i < j, in body-array order.
 *     Newton's third law writes -F_ij into body i's accumulator and
 *     +F_ij into body j's, so each pair contributes exactly once
 *     and the order is deterministic.
 *   - For n_bodies ≥ K26ASTRO_COMPENSATED_SUM_THRESHOLD (16),
 *     Neumaier compensated summation runs on each body's accumulator
 *     so cross-platform IEEE-754-bit-identity holds in portable mode.
 *
 * Position handling: K26AstroBody carries K26AstroPos (sector grid +
 * local offset). The pairwise displacement is computed via
 * k26astro_pos_sub which returns a K26V3 in metres relative to the
 * caller-supplied reference. This means the inner loop works in
 * metres-scale doubles without losing precision past the first
 * sector boundary. */
#include "k26astro_grav/forces.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/perturb.h"
#include "k26astro_core/sum.h"
#include "k26astro_core/pos.h"

#include <math.h>
#include <string.h>

static void zero_v3_array_(K26V3 *a, int n)
{
    memset(a, 0, sizeof(K26V3) * (size_t)n);
}

/* Plain (uncompensated) direct N² for small body counts. */
static void force_direct_plain_(const K26AstroGravView *view,
                                 double soft2, K26V3 *accel)
{
    int n = view->n;
    const K26AstroBody *b = view->bodies;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            K26V3 r_ij = k26astro_pos_sub(&b[j].pos, &b[i].pos);
            double r2 = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z + soft2;
            double r_mag = sqrt(r2);
            double inv_r3 = 1.0 / (r2 * r_mag);
            /* a_i += +GM_j / r² · r̂  (toward j) */
            double s_i =  b[j].gm * inv_r3;
            double s_j = -b[i].gm * inv_r3;
            accel[i].x += s_i * r_ij.x;
            accel[i].y += s_i * r_ij.y;
            accel[i].z += s_i * r_ij.z;
            accel[j].x += s_j * r_ij.x;
            accel[j].y += s_j * r_ij.y;
            accel[j].z += s_j * r_ij.z;
        }
    }
}

/* Compensated direct N² for larger body counts. Each body's three
 * accumulators are K26AstroSum (sum + comp); finalised at the end. */
static void force_direct_compensated_(const K26AstroGravView *view,
                                       double soft2, K26V3 *accel)
{
    int n = view->n;
    const K26AstroBody *b = view->bodies;

    /* 3 compensated accumulators per body; could heap-alloc for huge
     * N, but n_bodies is bounded above by the static body array
     * (caller-managed); stack-alloc keeps the hot path malloc-free
     * for the typical inner-SS scale. We support up to n = 256 here.
     * Larger N would need the heap path (and ultimately a Barnes-Hut
     * approximation). */
    enum { MAX_COMP_N = 256 };
    K26AstroSum acc[MAX_COMP_N][3];
    int use_n = (n < MAX_COMP_N) ? n : MAX_COMP_N;
    for (int i = 0; i < use_n; i++) {
        for (int k = 0; k < 3; k++) k26astro_sum_init(&acc[i][k]);
    }

    /* For n > MAX_COMP_N, fall back to plain summation for the
     * overflow. The threshold is sized for v0.1 use cases. */
    for (int i = 0; i < use_n; i++) {
        for (int j = i + 1; j < use_n; j++) {
            K26V3 r_ij = k26astro_pos_sub(&b[j].pos, &b[i].pos);
            double r2 = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z + soft2;
            double r_mag = sqrt(r2);
            double inv_r3 = 1.0 / (r2 * r_mag);
            double s_i =  b[j].gm * inv_r3;
            double s_j = -b[i].gm * inv_r3;
            k26astro_sum_add(&acc[i][0], s_i * r_ij.x);
            k26astro_sum_add(&acc[i][1], s_i * r_ij.y);
            k26astro_sum_add(&acc[i][2], s_i * r_ij.z);
            k26astro_sum_add(&acc[j][0], s_j * r_ij.x);
            k26astro_sum_add(&acc[j][1], s_j * r_ij.y);
            k26astro_sum_add(&acc[j][2], s_j * r_ij.z);
        }
    }
    for (int i = 0; i < use_n; i++) {
        accel[i].x = k26astro_sum_final(&acc[i][0]);
        accel[i].y = k26astro_sum_final(&acc[i][1]);
        accel[i].z = k26astro_sum_final(&acc[i][2]);
    }
    /* Overflow tail (rare). */
    if (n > MAX_COMP_N) {
        K26AstroGravView tail = { view->bodies + MAX_COMP_N, n - MAX_COMP_N };
        (void)tail;
        /* Not implemented in v0.1; document the limit. */
    }
}

void k26astro_grav_force_direct(const K26AstroGravView *view,
                                K26V3 *accel_out)
{
    if (!view || !accel_out || view->n < 1) return;
    zero_v3_array_(accel_out, view->n);
    if (view->n >= K26ASTRO_COMPENSATED_SUM_THRESHOLD) {
        force_direct_compensated_(view, 0.0, accel_out);
    } else {
        force_direct_plain_(view, 0.0, accel_out);
    }
}

void k26astro_grav_force_direct_softened(const K26AstroGravView *view,
                                          double softening,
                                          K26V3 *accel_out)
{
    if (!view || !accel_out || view->n < 1) return;
    zero_v3_array_(accel_out, view->n);
    double soft2 = softening * softening;
    if (view->n >= K26ASTRO_COMPENSATED_SUM_THRESHOLD) {
        force_direct_compensated_(view, soft2, accel_out);
    } else {
        force_direct_plain_(view, soft2, accel_out);
    }
}

/* ---- MERCURIUS pair-by-pair force decomposition ----------------- */

/* Look up the K weight for pair (i, j) in the weight list. Linear
 * scan; encounters typically number 0-20 even in dense systems.
 * Returns 0.0 if the pair is not in the encounter list. */
static double mercurius_pair_weight_(
        const K26AstroPairWeight *weights, int n_weights, int i, int j)
{
    /* Canonical pair order: smaller index first. */
    if (i > j) { int t = i; i = j; j = t; }
    for (int k = 0; k < n_weights; k++) {
        int wi = weights[k].i, wj = weights[k].j;
        if (wi > wj) { int t = wi; wi = wj; wj = t; }
        if (wi == i && wj == j) return weights[k].k_weight;
    }
    return 0.0;
}

void k26astro_grav_force_direct_weighted(
        const K26AstroGravView    *view,
        const K26AstroPairWeight  *weights, int n_weights,
        K26V3 *out_far, K26V3 *out_near)
{
    if (!view || view->n < 1) return;
    int n = view->n;
    const K26AstroBody *b = view->bodies;

    if (out_far)  zero_v3_array_(out_far, n);
    if (out_near) zero_v3_array_(out_near, n);

    /* Same pair iteration order as force_direct_plain_ so the
     * accumulation order is deterministic and matches the
     * unweighted path. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            K26V3 r_ij = k26astro_pos_sub(&b[j].pos, &b[i].pos);
            double r2 = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
            double r_mag = sqrt(r2);
            double inv_r3 = 1.0 / (r2 * r_mag);
            double s_i =  b[j].gm * inv_r3;
            double s_j = -b[i].gm * inv_r3;
            double K_ij = mercurius_pair_weight_(weights, n_weights, i, j);
            double w_far  = 1.0 - K_ij;
            double w_near = K_ij;
            if (out_far) {
                out_far[i].x += w_far * s_i * r_ij.x;
                out_far[i].y += w_far * s_i * r_ij.y;
                out_far[i].z += w_far * s_i * r_ij.z;
                out_far[j].x += w_far * s_j * r_ij.x;
                out_far[j].y += w_far * s_j * r_ij.y;
                out_far[j].z += w_far * s_j * r_ij.z;
            }
            if (out_near) {
                out_near[i].x += w_near * s_i * r_ij.x;
                out_near[i].y += w_near * s_i * r_ij.y;
                out_near[i].z += w_near * s_i * r_ij.z;
                out_near[j].x += w_near * s_j * r_ij.x;
                out_near[j].y += w_near * s_j * r_ij.y;
                out_near[j].z += w_near * s_j * r_ij.z;
            }
        }
    }
}

/* Single-mode weighted direct N² writing into a single accel array.
 * Used internally by accel_total when state->mercurius is active. */
static void force_direct_mercurius_(
        const K26AstroGravView *view,
        const K26AstroMercuriusContext *ctx,
        K26V3 *accel)
{
    int n = view->n;
    const K26AstroBody *b = view->bodies;
    int near_mode = (ctx->mode == K26ASTRO_MERCURIUS_NEAR);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            K26V3 r_ij = k26astro_pos_sub(&b[j].pos, &b[i].pos);
            double r2 = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
            double r_mag = sqrt(r2);
            double inv_r3 = 1.0 / (r2 * r_mag);
            double s_i =  b[j].gm * inv_r3;
            double s_j = -b[i].gm * inv_r3;
            double K_ij = mercurius_pair_weight_(ctx->pair_weights,
                                                  ctx->n_pair_weights, i, j);
            double w = near_mode ? K_ij : (1.0 - K_ij);
            if (w == 0.0) continue;
            accel[i].x += w * s_i * r_ij.x;
            accel[i].y += w * s_i * r_ij.y;
            accel[i].z += w * s_i * r_ij.z;
            accel[j].x += w * s_j * r_ij.x;
            accel[j].y += w * s_j * r_ij.y;
            accel[j].z += w * s_j * r_ij.z;
        }
    }
}

/* Total acceleration: direct N² + enabled built-in perturbations +
 * user-registered perturbations. Each perturbation ADDS into accel_out;
 * direct N² zeros and writes; perturbations append. */
void k26astro_grav_accel_total(const K26AstroGravState *state,
                                K26V3 *accel_out)
{
    if (!state || !accel_out) return;
    K26AstroGravView view = { state->bodies, state->n_bodies };

    if (state->mercurius) {
        /* Paper-faithful MERCURIUS pair-split: only the active
         * mode's portion of the direct N² contributes. Softening
         * applied via the standard inverse-distance form (the
         * MERCURIUS handoff is independent of Plummer softening;
         * softening is for clusters, not close encounters). */
        zero_v3_array_(accel_out, state->n_bodies);
        force_direct_mercurius_(&view, state->mercurius, accel_out);
        /* Perturbations only contribute to FAR; they are smooth
         * additive corrections that participate with WH, not the
         * IAS15 near sub-step. */
        if (state->mercurius->mode != K26ASTRO_MERCURIUS_FAR) goto perturbs_done;
    } else if (state->softening > 0.0) {
        k26astro_grav_force_direct_softened(&view, state->softening,
                                              accel_out);
    } else {
        k26astro_grav_force_direct(&view, accel_out);
    }

    /* Built-in perturbations dispatched via the registry's flag gate.
     * Each perturbation reads state->use_<flag> and ADDS into accel_out. */
    void k26astro_perturb_j2     (const K26AstroGravState *,
                                   const K26AstroGravView *,
                                   K26V3 *, void *);
    void k26astro_perturb_srp    (const K26AstroGravState *,
                                   const K26AstroGravView *,
                                   K26V3 *, void *);
    void k26astro_perturb_gr_ppn1(const K26AstroGravState *,
                                   const K26AstroGravView *,
                                   K26V3 *, void *);

    if (state->use_j2)      k26astro_perturb_j2     (state, &view, accel_out, NULL);
    if (state->use_srp)     k26astro_perturb_srp    (state, &view, accel_out, NULL);
    if (state->use_gr_ppn1) k26astro_perturb_gr_ppn1(state, &view, accel_out, NULL);

    /* User-registered perturbations (called in registration order). */
    if (state->perturbs && state->perturbs->count > 0) {
        for (int k = 0; k < state->perturbs->count; k++) {
            state->perturbs->fns[k](state, &view, accel_out,
                                     state->perturbs->ctxs[k]);
        }
    }
perturbs_done:
    return;
}
