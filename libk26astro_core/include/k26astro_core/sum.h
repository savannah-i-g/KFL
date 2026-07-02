/* k26astro_core/sum.h — Neumaier (Kahan-Babuška) compensated summation.
 *
 * The astro suite's force-accumulation loops are the hottest source of
 * cross-platform non-determinism: floating-point summation is order-
 * dependent, and adding tiny perturbations to large primary accelerations
 * loses bits to cancellation. Neumaier's variant of compensated summation
 * (Neumaier 1974, ZAMM 54:1-2) tracks the running rounding error in a
 * companion accumulator and folds it back at the end. With doubles, the
 * effective precision climbs from ~52 bits to ~104 — enough to make N²
 * force summation bit-identical across compilers and CPU vendors when the
 * surrounding FPU state is also pinned (FE_TONEAREST + DAZ/FTZ off, see
 * libk26astro_grav/src/fpu_pin.c).
 *
 * Usage:
 *
 *     K26AstroSum s;
 *     k26astro_sum_init(&s);
 *     for (int i = 0; i < n; ++i) k26astro_sum_add(&s, terms[i]);
 *     double total = k26astro_sum_final(&s);
 *
 * The threefold V3 form is provided as a convenience for the integrators'
 * accel accumulation:
 *
 *     K26AstroSum acc[3];
 *     for (int i = 0; i < 3; ++i) k26astro_sum_init(&acc[i]);
 *     for (int j = 0; j < n; ++j) k26astro_sum_add_v3(acc, &per_body_accel[j]);
 *     K26V3 total = k26astro_sum_final_v3(acc);
 *
 * The threshold K26ASTRO_COMPENSATED_SUM_THRESHOLD is the body-count
 * cut-off above which integrators switch into compensated mode; below it
 * plain summation is bit-identical anyway (no cancellation) and avoiding
 * the extra arithmetic keeps the inner loop fast.
 *
 * Reference: Neumaier, A. (1974) "Rundungsfehleranalyse einiger
 * Verfahren zur Summation endlicher Summen" ZAMM 54:39-51. */
#ifndef K26ASTRO_CORE_SUM_H
#define K26ASTRO_CORE_SUM_H

#include "k26m3d.h"   /* K26V3 */

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_COMPENSATED_SUM_THRESHOLD 16

typedef struct K26AstroSum {
    double sum;
    double comp;
} K26AstroSum;

static inline void k26astro_sum_init(K26AstroSum *s)
{
    s->sum  = 0.0;
    s->comp = 0.0;
}

void   k26astro_sum_add(K26AstroSum *s, double x);
double k26astro_sum_final(const K26AstroSum *s);

void  k26astro_sum_add_v3(K26AstroSum s[3], const K26V3 *v);
K26V3 k26astro_sum_final_v3(const K26AstroSum s[3]);

#ifdef __cplusplus
}
#endif

#endif
