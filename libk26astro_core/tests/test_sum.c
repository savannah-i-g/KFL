/* test_sum.c — Neumaier compensated-summation precision tests.
 *
 * Acceptance:
 *   - 1e6 ones summed to a 1e8 accumulator: naive loses ~5 bits;
 *     compensated recovers exactly 1e6 in the residual.
 *   - Adversarial Kahan-failure case (large+tiny+(-large)): naive
 *     drops the tiny; compensated keeps it.
 *   - V3 form returns per-axis-correct totals.
 *   - Threshold constant equals the documented value (16).
 */
#include "k26astro_core/sum.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* ---- Threshold constant -------------------------------- */
    assert(K26ASTRO_COMPENSATED_SUM_THRESHOLD == 16);

    /* ---- Basic add + finalise ------------------------------ */
    K26AstroSum s;
    k26astro_sum_init(&s);
    assert(k26astro_sum_final(&s) == 0.0);

    k26astro_sum_add(&s, 1.0);
    k26astro_sum_add(&s, 2.0);
    k26astro_sum_add(&s, 3.0);
    assert(k26astro_sum_final(&s) == 6.0);

    /* ---- Cancellation-recovery test ------------------------ *
     * Add 1.0 a million times, starting from a 1e16-magnitude
     * accumulator. Naive sum loses the 1e6 entirely (1e16 + 1 == 1e16
     * in double); compensated should recover 1.0e6 exactly. */
    k26astro_sum_init(&s);
    k26astro_sum_add(&s, 1.0e16);
    for (int i = 0; i < 1000000; i++) {
        k26astro_sum_add(&s, 1.0);
    }
    k26astro_sum_add(&s, -1.0e16);
    double recovered = k26astro_sum_final(&s);
    /* Naive would give 0.0; compensated gives 1e6 (or within ULPs). */
    assert(fabs(recovered - 1.0e6) < 1.0);

    /* ---- Neumaier's specific advantage over Kahan ---------- *
     * The order (small, large, -large) traps Kahan but works for
     * Neumaier. Add 1.0, then 1e100, then -1e100: result should be 1.0. */
    k26astro_sum_init(&s);
    k26astro_sum_add(&s, 1.0);
    k26astro_sum_add(&s, 1.0e100);
    k26astro_sum_add(&s, -1.0e100);
    double r = k26astro_sum_final(&s);
    assert(fabs(r - 1.0) < 1e-10);

    /* ---- V3 accumulation ---------------------------------- */
    K26AstroSum acc[3];
    for (int i = 0; i < 3; i++) k26astro_sum_init(&acc[i]);

    K26V3 perturbs[5] = {
        { 1.0,  2.0,  3.0 },
        { 0.1, -0.2,  0.05 },
        { 1e-9, 1e-9, 1e-9 },
        { -1.0, -2.0, -3.0 },
        { 0.0,  0.0,  1e-15 }
    };
    for (int i = 0; i < 5; i++) k26astro_sum_add_v3(acc, &perturbs[i]);

    K26V3 total = k26astro_sum_final_v3(acc);
    /* Per-axis expected: x = 0.1 + 1e-9, y = -0.2 + 1e-9,
     * z = 0.05 + 1e-9 + 1e-15. */
    assert(fabs(total.x - (0.1  + 1e-9))         < 1e-12);
    assert(fabs(total.y - (-0.2 + 1e-9))         < 1e-12);
    assert(fabs(total.z - (0.05 + 1e-9 + 1e-15)) < 1e-12);

    /* ---- Determinism under reordering --------------------- *
     * Adding the same set in two different orders should produce
     * bit-identical results in compensated form (the whole point of
     * the compensation). */
    double terms[8] = { 1e10, 1.0, -1e10, 0.5, 1e-8, 1e-9, 1e-10, 1e-11 };
    K26AstroSum a, b;
    k26astro_sum_init(&a);
    k26astro_sum_init(&b);
    for (int i = 0; i < 8; i++)         k26astro_sum_add(&a, terms[i]);
    for (int i = 7; i >= 0; i--)        k26astro_sum_add(&b, terms[i]);
    /* Forward vs reverse order: should agree to many more bits than
     * naive sum would. (Bit-identical is too strong — orderings differ
     * in which path through the branch is taken; both are accurate.) */
    double diff = fabs(k26astro_sum_final(&a) - k26astro_sum_final(&b));
    assert(diff < 1e-12);

    printf("test_sum: OK (compensated sum recovers 1e6 over 1e16)\n");
    return 0;
}
