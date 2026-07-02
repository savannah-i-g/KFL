/* ias15_coeffs.c — Gauss-Radau nodes (hex-literal) + runtime-built
 * conversion matrices.
 *
 * The Radau-7 positive-endpoint node positions h[1..7] are stored as
 * hex-literal doubles so the integer bit pattern is bit-exact across
 * compilers (strtod-rounding-independent). The decimal comment is
 * the human-readable equivalent to 25 sig figs from Rein & Spiegel
 * 2015 (and the canonical REBOUND table at
 * src/integrator_ias15.c).
 *
 * Derived matrices c[], d[], r[] are computed once at first IAS15
 * step from the h[] table — this avoids hand-transcription errors
 * in the long lower-triangle packs that would otherwise be
 * transcribed from REBOUND source. Since h[] is fully determined,
 * the derived matrices are deterministic to within IEEE-754 rounding
 * (and bit-identical across compilers if the host arithmetic is
 * pinned, per the determinism contract — -ffp-contract=off, FPU
 * rounding mode pinned). */
#include "ias15_internal.h"
#include "k26astro_grav/ias15.h"

#include <stdlib.h>
#include <string.h>

/* Canonical Gauss-Radau positive-endpoint nodes (n=7).
 *
 * Originally transcribed as hex-literals for compiler-independence
 * but the hex transcription introduced a ~0.02% error in the
 * lower-order h values (h_1 came out 0.05625003 instead of 0.05626256)
 * which propagated through the c matrix and capped IAS15 at ~3rd
 * order accuracy. v0.1 uses decimal literals (15 sig figs); strtod
 * on glibc and musl both round-trip these to within 1 ULP so the
 * cross-libc determinism contract holds at 2e-16 precision. */
const double k26_ias15_h[8] = {
    0.0,
    0.0562625605369221464656522,
    0.180240691736892364987580,
    0.352624717113169637373908,
    0.547153626330555383001449,
    0.734210177215410531523211,
    0.885320946839095768090360,
    0.977520613561287501891174
};

/* Derived matrices, packed lower-triangle. */
static double c_lt[21];
static double d_lt[21];
static double r_lt[28];
static int    matrices_built = 0;

/* Lower-triangle index helpers.
 *
 * c, d are indexed as c[(k*(k-1))/2 + i] for k > i, k ∈ 1..6, i < k.
 *   Total: sum_{k=1..6} k = 21 entries.
 *
 * r is indexed as r[(i*(i-1))/2 + j] for i > j, i ∈ 1..7, j < i.
 *   Total: sum_{i=1..7} i = 28 entries. */
static inline int idx_lt2_(int k, int i) { return (k * (k - 1)) / 2 + i; }
static inline int idx_r_(int i, int j)   { return (i * (i - 1)) / 2 + j; }

void k26_ias15_init_matrices(void)
{
    if (matrices_built) return;

    /* r[i,j] = 1 / (h[i] - h[j]) for i > j */
    for (int i = 1; i <= 7; i++) {
        for (int j = 0; j < i; j++) {
            double denom = k26_ias15_h[i] - k26_ias15_h[j];
            r_lt[idx_r_(i, j)] = 1.0 / denom;
        }
    }

    /* c[k][i] are the Newton-form elementary symmetric polynomials of
     * h[1..k], signed alternately. For row k, c[k][0]..c[k][k-1] are
     * such that:
     *   product_{j=1..k} (x - h[j]) = sum_{i=0..k} c'[k][i] x^(k-i)
     * with c'[k][0] = 1 and c[k][i] = c'[k][i+1] (the coefficients of
     * lower-order terms).
     *
     * Generate via the standard polynomial-expansion recurrence:
     *   c[k][i] = c[k-1][i-1] - h[k]·c[k-1][i]
     * starting from c[0] = empty (the constant polynomial 1). */

    /* Working buffer for elementary-symm coefficients (1-indexed in
     * the polynomial sense). We keep them in a 2-D array of [row][col]
     * for clarity, then pack to the lower-triangle store. */
    double coef[7][7];
    memset(coef, 0, sizeof(coef));
    coef[0][0] = -k26_ias15_h[1];
    /* k = 1: polynomial is (x - h_1), so c'[1][0] = 1, c'[1][1] = -h_1.
     * Our c_lt[0] = c'[1][1] = -h_1.  */
    c_lt[0] = -k26_ias15_h[1];

    for (int k = 2; k <= 6; k++) {
        /* Build c'[k][i] from c'[k-1][i-1] - h[k]·c'[k-1][i].
         * Here coef[k-1][i] holds q_k(k-1-i), i.e. the coefficient
         * of τ^(k-1-i) in (τ-h_1)...(τ-h_k) — high degree first. */
        coef[k - 1][0] = coef[k - 2][0] - k26_ias15_h[k];      /* x^(k-1) coeff */
        for (int i = 1; i < k - 1; i++) {
            coef[k - 1][i] = coef[k - 2][i] - k26_ias15_h[k] * coef[k - 2][i - 1];
        }
        coef[k - 1][k - 1] = -k26_ias15_h[k] * coef[k - 2][k - 2];  /* constant term */

        /* Pack with index inversion: c_lt[idx_lt2_(k,i)] = q_k(i) =
         * coef[k-1][k-1-i]. */
        for (int i = 0; i < k; i++) {
            c_lt[idx_lt2_(k, i)] = coef[k - 1][k - 1 - i];
        }
    }

    /* d[k][i] are the inverse transform — products of h's (without
     * sign alternation). Generate via the dual recurrence:
     *   d[k][i] = d[k-1][i-1] + h[k]·d[k-1][i]
     * with d[1][0] = h_1. */
    double dcoef[7][7];
    memset(dcoef, 0, sizeof(dcoef));
    dcoef[0][0] = k26_ias15_h[1];
    d_lt[0] = k26_ias15_h[1];

    for (int k = 2; k <= 6; k++) {
        dcoef[k - 1][0] = dcoef[k - 2][0] + k26_ias15_h[k];
        for (int i = 1; i < k - 1; i++) {
            dcoef[k - 1][i] = dcoef[k - 2][i] + k26_ias15_h[k] * dcoef[k - 2][i - 1];
        }
        dcoef[k - 1][k - 1] = k26_ias15_h[k] * dcoef[k - 2][k - 2];
        for (int i = 0; i < k; i++) {
            d_lt[idx_lt2_(k, i)] = dcoef[k - 1][i];
        }
    }

    matrices_built = 1;
}

const double *k26_ias15_c(void)  { return c_lt; }
const double *k26_ias15_d(void)  { return d_lt; }
const double *k26_ias15_rr(void) { return r_lt; }

/* ---- Carry buffer alloc / release ------------------------------ */
int k26_ias15_carry_alloc(K26AstroIAS15Carry **out, int n)
{
    K26AstroIAS15Carry *c = calloc(1, sizeof(K26AstroIAS15Carry));
    if (!c) return K26ASTRO_E_ALLOC;
    for (int k = 0; k < 7; k++) {
        c->b[k] = calloc((size_t)n, sizeof(K26V3));
        c->e[k] = calloc((size_t)n, sizeof(K26V3));
        c->g[k] = calloc((size_t)n, sizeof(K26V3));
        if (!c->b[k] || !c->e[k] || !c->g[k]) {
            k26_ias15_carry_release(c); return K26ASTRO_E_ALLOC;
        }
    }
    c->at0   = calloc((size_t)n, sizeof(K26V3));
    c->r_sub = calloc((size_t)n, sizeof(K26V3));
    c->v_sub = calloc((size_t)n, sizeof(K26V3));
    c->a_sub = calloc((size_t)n, sizeof(K26V3));
    if (!c->at0 || !c->r_sub || !c->v_sub || !c->a_sub) {
        k26_ias15_carry_release(c); return K26ASTRO_E_ALLOC;
    }
    c->capacity = n;
    c->initialised = 0;
    c->dt_proposed = 0.0;
    *out = c;
    return K26ASTRO_E_OK;
}

void k26_ias15_carry_release(K26AstroIAS15Carry *c)
{
    if (!c) return;
    for (int k = 0; k < 7; k++) {
        free(c->b[k]); free(c->e[k]); free(c->g[k]);
    }
    free(c->at0); free(c->r_sub); free(c->v_sub); free(c->a_sub);
    free(c);
}
