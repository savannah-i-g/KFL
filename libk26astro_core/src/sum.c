/* sum.c — Neumaier compensated summation impl.
 *
 * FMA-off discipline arrives via Makefile CFLAGS (-ffp-contract=off);
 * the C11 `#pragma STDC FP_CONTRACT OFF` is incompletely implemented
 * in GCC (warns -Wunknown-pragmas) so we rely on the build flag. The
 * compiler must not contract the (t - sum) - x sequence below, or the
 * compensation collapses to a no-op. */

#include "k26astro_core/sum.h"

#include <math.h>

/* Single-scalar add. Neumaier's improvement over Kahan is that the
 * compensation term is selected by *magnitude* rather than always
 * the previous round-off: if |sum| >= |x| the low-order bits of x
 * fall into comp; if |x| > |sum| the low-order bits of sum do. This
 * fixes Kahan's failure mode when a tiny accumulator absorbs a large
 * value. */
void k26astro_sum_add(K26AstroSum *s, double x)
{
    double t = s->sum + x;
    if (fabs(s->sum) >= fabs(x)) {
        s->comp += (s->sum - t) + x;
    } else {
        s->comp += (x - t) + s->sum;
    }
    s->sum = t;
}

double k26astro_sum_final(const K26AstroSum *s)
{
    return s->sum + s->comp;
}

void k26astro_sum_add_v3(K26AstroSum s[3], const K26V3 *v)
{
    k26astro_sum_add(&s[0], v->x);
    k26astro_sum_add(&s[1], v->y);
    k26astro_sum_add(&s[2], v->z);
}

K26V3 k26astro_sum_final_v3(const K26AstroSum s[3])
{
    K26V3 r;
    r.x = k26astro_sum_final(&s[0]);
    r.y = k26astro_sum_final(&s[1]);
    r.z = k26astro_sum_final(&s[2]);
    return r;
}
