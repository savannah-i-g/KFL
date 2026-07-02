/* k26astro_core/pos.h — sector-grid position type for solar-system-
 * scale simulation.
 *
 * Why this is load-bearing
 * ------------------------
 * IEEE-754 binary64 from a single origin runs out of usable precision
 * past the asteroid belt within a decade of integrator time. At
 * Pluto (~6e12 m) the absolute precision is ~1.3 mm; integration
 * drift consumes the rest within a few hundred years. The sector
 * grid sidesteps the problem by keeping the working coordinate in a
 * small fold of binary64 (within-sector offset) plus an exact integer
 * sector index.
 *
 * Layout
 * ------
 * Each component is split:
 *   sx, sy, sz : int64 sector indices (always normalised so
 *                |lx|,|ly|,|lz| < SECTOR_EDGE/2 after a substep)
 *   lx, ly, lz : double local offset, in metres, inside the sector
 *
 * The sector edge is 2^36 m ≈ 6.872e10 m ≈ 0.459 AU. The size is
 * chosen as a compromise: small enough to keep local-coord
 * arithmetic at sub-mm precision, large enough that Earth-Mars
 * opposition only spans ~1-2 sectors.
 *
 * Substraction collapses to a precision-preserving relative vector:
 *   r = (a.s - b.s) * S + (a.l - b.l)
 * The integer subtraction is exact; the double residual is bounded
 * by 2*S so its precision is ~S * 2^-52 ≈ 23 μm. No precision is
 * lost regardless of where in the solar system the two bodies sit.
 *
 * Q64.64 fixed-point fallback
 * ---------------------------
 * `K26AstroPosFx` is opt-in for "I need bitwise-identical runs across
 * CPU vendors and even FPU modes". ±9.2e18 m range (≈970 ly),
 * precision ~5.4e-20 m. Add/sub via int128 arithmetic; no FPU mode,
 * no FMA, no reordering effects. ~2-3× slower than the sector path.
 */
#ifndef K26ASTRO_CORE_POS_H
#define K26ASTRO_CORE_POS_H

#include <stdint.h>

#include "k26m3d.h"   /* K26V3 — the binary64 vector */

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Sector grid ------------------------------------------------- */

#define K26ASTRO_SECTOR_EDGE_LOG2  36
#define K26ASTRO_SECTOR_EDGE_M     ((double)(1ULL << K26ASTRO_SECTOR_EDGE_LOG2))
/* 2^36 m ≈ 6.872 × 10^10 m ≈ 0.459 AU */

typedef struct {
    int64_t sx, sy, sz;   /* sector indices */
    double  lx, ly, lz;   /* local offset in metres, |.| < EDGE/2 after normalise */
} K26AstroPos;
/* sizeof == 48 (3*8 + 3*8). Distinct from K26V3 (24 bytes) — the type
 * system catches confusion between camera-space K26V3 and world-space
 * K26AstroPos at the call site. */

/* Zero-position constructor. */
K26AstroPos k26astro_pos_zero(void);

/* Construct from a raw metres-from-origin triple. Used for tests +
 * for small-scale (single-sector) scenes where the author is happy
 * to start with origin-relative inputs. Normalises the result. */
K26AstroPos k26astro_pos_from_m(double x, double y, double z);

/* Folds |lx|, |ly|, |lz| back into the [-EDGE/2, EDGE/2) range,
 * carrying the overflow into the sector index. Idempotent. Call
 * after every integrator substep on every body. */
void k26astro_pos_normalise(K26AstroPos *p);

/* Relative displacement vector: `a - b`. Returns a K26V3 in metres,
 * suitable for force-pair evaluation. Integer sector subtraction is
 * exact; the local residual is one binary64 subtract in [-EDGE, EDGE],
 * with precision ~EDGE × 2^-52 ≈ 23 μm. No precision is lost in the
 * relative vector even when both bodies are at Pluto distance from
 * the heliocentre.
 *
 * Pre-condition: both `a` and `b` are normalised. */
K26V3 k26astro_pos_sub(const K26AstroPos *a, const K26AstroPos *b);

/* In-place add of a K26V3 delta. Result is normalised. Useful for
 * the integrator's per-substep position update:
 *   k26astro_pos_add(&body->pos, k26m3d_v3_scale(body->vel, dt));
 */
void k26astro_pos_add(K26AstroPos *p, K26V3 delta);

/* In-place add of another K26AstroPos. Result is normalised.
 * Less common; useful for centre-of-mass accumulation. */
void k26astro_pos_add_pos(K26AstroPos *p, const K26AstroPos *q);

/* In-place scale of the position by a scalar. Result is normalised. */
void k26astro_pos_scale(K26AstroPos *p, double s);

/* Squared distance |a - b|² in metres². Cheaper than sqrt'ing
 * |sub| when only the magnitude-for-comparison is needed. */
double k26astro_pos_dist_sq(const K26AstroPos *a, const K26AstroPos *b);

/* Magnitude of |a - b| in metres. */
double k26astro_pos_dist(const K26AstroPos *a, const K26AstroPos *b);

/* Approximate metres-from-origin as a K26V3 — for diagnostic display
 * + render submission *only* when the camera is near origin. Loses
 * precision past the first sector boundary. Use `k26astro_pos_sub`
 * with a camera position for actual rendering. */
K26V3 k26astro_pos_to_m_approx(const K26AstroPos *p);

/* ---- Q64.64 fixed-point fallback -------------------------------- */
/* Opt-in determinism mode. Bitwise-identical add/sub across CPU
 * vendors via __int128_t intermediates. Range ±9.2e18 m (~970 ly);
 * precision ~5.4e-20 m. */

typedef struct {
    int64_t  hi;
    uint64_t lo;
} K26AstroQ6464;

typedef struct {
    K26AstroQ6464 x, y, z;
} K26AstroPosFx;

K26AstroPosFx k26astro_pos_fx_zero(void);
K26AstroPosFx k26astro_pos_fx_from_m(double x, double y, double z);
K26V3         k26astro_pos_fx_sub(const K26AstroPosFx *a,
                                  const K26AstroPosFx *b);
void          k26astro_pos_fx_add(K26AstroPosFx *p, K26V3 delta);

/* Lossy round-trip helpers for inter-conversion. The sector grid
 * has ~23 μm within-sector precision; the Q64.64 has ~5e-20 m. The
 * sector → Q64.64 direction loses nothing (well, loses ~23 μm
 * already inherent to the source). The other direction loses
 * everything below 23 μm. */
K26AstroPosFx k26astro_pos_to_fx(const K26AstroPos *p);
K26AstroPos   k26astro_pos_from_fx(const K26AstroPosFx *p);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CORE_POS_H */
