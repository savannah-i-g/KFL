/* Ramer-Douglas-Peucker polyline simplification.
 *
 * Iterative stack-based implementation: never recurses, so the stack
 * depth is bounded by O(log n) for well-conditioned inputs and O(n)
 * for pathological ones. Operates on packed (x, z) double pairs in
 * place; the survivors are compacted to the front and the new count
 * is returned.
 *
 * Reference: mapConverter.ts:151-201. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

/* Perpendicular distance from point p to the line through a..b.
 * Returns the absolute distance in (x, z) coords. Degenerate segments
 * (a == b) collapse to |p - a|. */
static double perp_dist(double px, double pz,
                        double ax, double az,
                        double bx, double bz)
{
    const double dx = bx - ax;
    const double dz = bz - az;
    const double len2 = dx * dx + dz * dz;
    if (len2 == 0.0) {
        const double ex = px - ax, ez = pz - az;
        return sqrt(ex * ex + ez * ez);
    }
    /* |( (p-a) × (b-a) )| / |b-a| */
    const double cross = (px - ax) * dz - (pz - az) * dx;
    return fabs(cross) / sqrt(len2);
}

size_t k26geo_rdp_xz(double *coords, size_t n, double epsilon)
{
    if (!coords)        return 0;
    if (n  <  2)        return n;
    if (epsilon <= 0.0) return n;

    /* Mark survivors via a separate byte buffer; endpoints always 1.
     * Stack of (lo, hi) index pairs to inspect. */
    unsigned char *keep = (unsigned char *)calloc(n, 1);
    if (!keep) return n;   /* OOM degrades gracefully to no-op */

    typedef struct { size_t lo, hi; } Span;
    Span *stk = (Span *)malloc(sizeof(Span) * n);
    if (!stk) { free(keep); return n; }
    size_t sp = 0;

    keep[0]     = 1;
    keep[n - 1] = 1;
    stk[sp++].lo = 0; stk[sp - 1].hi = n - 1;

    while (sp > 0) {
        const Span s = stk[--sp];
        if (s.hi <= s.lo + 1) continue;

        const double ax = coords[2 * s.lo],     az = coords[2 * s.lo + 1];
        const double bx = coords[2 * s.hi],     bz = coords[2 * s.hi + 1];

        double best   = 0.0;
        size_t best_i = 0;
        for (size_t i = s.lo + 1; i < s.hi; i++) {
            const double d = perp_dist(coords[2 * i], coords[2 * i + 1],
                                       ax, az, bx, bz);
            if (d > best) { best = d; best_i = i; }
        }

        if (best > epsilon) {
            keep[best_i] = 1;
            stk[sp].lo = s.lo;   stk[sp].hi = best_i; sp++;
            stk[sp].lo = best_i; stk[sp].hi = s.hi;   sp++;
        }
    }

    /* Compact survivors to the front. */
    size_t out = 0;
    for (size_t i = 0; i < n; i++) {
        if (keep[i]) {
            if (out != i) {
                coords[2 * out]     = coords[2 * i];
                coords[2 * out + 1] = coords[2 * i + 1];
            }
            out++;
        }
    }

    free(stk);
    free(keep);
    return out;
}
