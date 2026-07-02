/* libk26compute — deterministic RNG.
 *
 * SplitMix64 for state advance: a simple, well-tested, splittable
 * PRNG (Steele, Lea, Flood 2014). Per-call output is independent of
 * call ordering as long as the seed advances deterministically. We
 * use it directly rather than going through xorshift / PCG because
 * it's the smallest implementation that passes BigCrush, and the
 * caller-supplied seed is the only state surface we need to advance
 * deterministically across machines.
 *
 * Box-Muller for normals: simple, deterministic, no rejection (so
 * given seed reproduces exactly), at the cost of being marginally
 * slower than the Marsaglia polar method on average. For research
 * reproducibility, exact-bit determinism wins. */

#include "k26compute.h"

#include <math.h>
#include <stdint.h>

void k26c_rng_init(K26CRng *r, uint64_t seed)
{
    if (!r) return;
    r->s = seed ? seed : 0xdeadbeefcafef00d;     /* nonzero default */
}

static uint64_t splitmix64_next(uint64_t *s)
{
    uint64_t z = (*s += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

double k26c_rng_uniform(K26CRng *r)
{
    if (!r) return 0.0;
    uint64_t u = splitmix64_next(&r->s);
    /* Convert top 53 bits into a [0, 1) double — same construction
     * std uniform_real_distribution uses, no bias. */
    return (double)(u >> 11) * (1.0 / (double)(1ull << 53));
}

double k26c_rng_normal(K26CRng *r, double mu, double sigma)
{
    if (!r) return mu;
    /* Box-Muller. We discard the second sample to keep state advance
     * 1:1 with caller invocations — essential for reproducibility
     * across implementations that may or may not cache the second
     * sample. The minor efficiency loss is acceptable. */
    double u1 = k26c_rng_uniform(r);
    double u2 = k26c_rng_uniform(r);
    if (u1 < 1e-300) u1 = 1e-300;             /* avoid log(0) */
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mu + sigma * z;
}
