/* pos.c — sector grid + Q64.64 implementations. */
#include "k26astro_core/pos.h"

#include <math.h>

/* Sector edge in metres, repeated as a local for tighter codegen. */
static const double S = K26ASTRO_SECTOR_EDGE_M;

/* ---- Sector grid ------------------------------------------------- */

K26AstroPos k26astro_pos_zero(void)
{
    K26AstroPos p = { 0, 0, 0, 0.0, 0.0, 0.0 };
    return p;
}

static void normalise_axis_(int64_t *sec, double *loc)
{
    /* Fold |loc| into [-S/2, S/2). Loop instead of a single division
     * for two reasons: (a) the loop body is usually executed zero
     * times after a small substep, (b) the integer-divide path
     * (lroundf etc.) sometimes generates worse code on -O0 / -Og. */
    while (*loc >=  0.5 * S) { *loc -= S; (*sec)++; }
    while (*loc <  -0.5 * S) { *loc += S; (*sec)--; }
}

void k26astro_pos_normalise(K26AstroPos *p)
{
    if (!p) return;
    normalise_axis_(&p->sx, &p->lx);
    normalise_axis_(&p->sy, &p->ly);
    normalise_axis_(&p->sz, &p->lz);
}

K26AstroPos k26astro_pos_from_m(double x, double y, double z)
{
    K26AstroPos p = { 0, 0, 0, x, y, z };
    k26astro_pos_normalise(&p);
    return p;
}

K26V3 k26astro_pos_sub(const K26AstroPos *a, const K26AstroPos *b)
{
    K26V3 r;
    r.x = (double)(a->sx - b->sx) * S + (a->lx - b->lx);
    r.y = (double)(a->sy - b->sy) * S + (a->ly - b->ly);
    r.z = (double)(a->sz - b->sz) * S + (a->lz - b->lz);
    return r;
}

void k26astro_pos_add(K26AstroPos *p, K26V3 d)
{
    if (!p) return;
    p->lx += d.x;
    p->ly += d.y;
    p->lz += d.z;
    k26astro_pos_normalise(p);
}

void k26astro_pos_add_pos(K26AstroPos *p, const K26AstroPos *q)
{
    if (!p || !q) return;
    p->sx += q->sx; p->sy += q->sy; p->sz += q->sz;
    p->lx += q->lx; p->ly += q->ly; p->lz += q->lz;
    k26astro_pos_normalise(p);
}

void k26astro_pos_scale(K26AstroPos *p, double s)
{
    if (!p) return;
    /* Compose total-metres-along-axis, scale, re-split. The integer
     * part of the product can overflow int64 for sufficiently large
     * |s * sx|; this routine is intended for small-magnitude scales
     * (e.g. centroid weighting). Callers that need large-scale
     * multiplication should drop to Q64.64. */
    double tx = (double)p->sx * S + p->lx;
    double ty = (double)p->sy * S + p->ly;
    double tz = (double)p->sz * S + p->lz;
    tx *= s; ty *= s; tz *= s;
    p->sx = 0; p->sy = 0; p->sz = 0;
    p->lx = tx; p->ly = ty; p->lz = tz;
    k26astro_pos_normalise(p);
}

double k26astro_pos_dist_sq(const K26AstroPos *a, const K26AstroPos *b)
{
    K26V3 r = k26astro_pos_sub(a, b);
    return r.x * r.x + r.y * r.y + r.z * r.z;
}

double k26astro_pos_dist(const K26AstroPos *a, const K26AstroPos *b)
{
    return sqrt(k26astro_pos_dist_sq(a, b));
}

K26V3 k26astro_pos_to_m_approx(const K26AstroPos *p)
{
    K26V3 r;
    r.x = (double)p->sx * S + p->lx;
    r.y = (double)p->sy * S + p->ly;
    r.z = (double)p->sz * S + p->lz;
    return r;
}

/* ---- Q64.64 fixed-point ----------------------------------------- */

/* 2^64 as a double — exact (positive normal). */
#define K26A_TWO64 1.8446744073709551616e19

static K26AstroQ6464 q_from_double_(double x)
{
    K26AstroQ6464 q;
    /* floor() works correctly across negatives (e.g. floor(-1.5) == -2).
     * After splitting, lo is the [0,1) fractional residual. */
    double f_hi  = floor(x);
    double f_lo  = x - f_hi;
    q.hi = (int64_t)f_hi;
    q.lo = (uint64_t)(f_lo * K26A_TWO64);
    return q;
}

static double q_to_double_(K26AstroQ6464 q)
{
    return (double)q.hi + (double)q.lo * (1.0 / K26A_TWO64);
}

/* Signed subtraction. Result fits unless inputs differ by more than
 * 2^127; that's well outside the calibrated range. */
static K26AstroQ6464 q_sub_(K26AstroQ6464 a, K26AstroQ6464 b)
{
    __int128_t a128 = ((__int128_t)a.hi << 64) | (__int128_t)a.lo;
    __int128_t b128 = ((__int128_t)b.hi << 64) | (__int128_t)b.lo;
    __int128_t r    = a128 - b128;
    K26AstroQ6464 q;
    q.hi = (int64_t)(r >> 64);
    q.lo = (uint64_t)r;
    return q;
}

static K26AstroQ6464 q_add_(K26AstroQ6464 a, K26AstroQ6464 b)
{
    __int128_t a128 = ((__int128_t)a.hi << 64) | (__int128_t)a.lo;
    __int128_t b128 = ((__int128_t)b.hi << 64) | (__int128_t)b.lo;
    __int128_t r    = a128 + b128;
    K26AstroQ6464 q;
    q.hi = (int64_t)(r >> 64);
    q.lo = (uint64_t)r;
    return q;
}

K26AstroPosFx k26astro_pos_fx_zero(void)
{
    K26AstroPosFx p;
    p.x.hi = 0; p.x.lo = 0;
    p.y.hi = 0; p.y.lo = 0;
    p.z.hi = 0; p.z.lo = 0;
    return p;
}

K26AstroPosFx k26astro_pos_fx_from_m(double x, double y, double z)
{
    K26AstroPosFx p;
    p.x = q_from_double_(x);
    p.y = q_from_double_(y);
    p.z = q_from_double_(z);
    return p;
}

K26V3 k26astro_pos_fx_sub(const K26AstroPosFx *a, const K26AstroPosFx *b)
{
    K26V3 r;
    r.x = q_to_double_(q_sub_(a->x, b->x));
    r.y = q_to_double_(q_sub_(a->y, b->y));
    r.z = q_to_double_(q_sub_(a->z, b->z));
    return r;
}

void k26astro_pos_fx_add(K26AstroPosFx *p, K26V3 d)
{
    if (!p) return;
    p->x = q_add_(p->x, q_from_double_(d.x));
    p->y = q_add_(p->y, q_from_double_(d.y));
    p->z = q_add_(p->z, q_from_double_(d.z));
}

K26AstroPosFx k26astro_pos_to_fx(const K26AstroPos *p)
{
    /* total = sx * S + lx ; convert via two adds rather than one
     * multiply-then-add so the high-magnitude integer contribution
     * is exact (S is a power of two, so sx*S in Q64.64 is just a
     * left shift). */
    K26AstroPosFx out = k26astro_pos_fx_zero();
    if (!p) return out;
    /* sx * 2^36 as Q64.64: (sx << 36) in the int128 sense. Since
     * Q64.64's hi is int64, sx*2^36 must fit in int64 — i.e.
     * |sx| < 2^27. That's exactly the calibrated range
     * (±2^63 / 2^36 = ±2^27 sectors). Beyond it, the conversion
     * truncates and the caller is using Q64.64 beyond its design
     * range. */
    out.x.hi = p->sx << K26ASTRO_SECTOR_EDGE_LOG2;
    out.y.hi = p->sy << K26ASTRO_SECTOR_EDGE_LOG2;
    out.z.hi = p->sz << K26ASTRO_SECTOR_EDGE_LOG2;
    out.x = q_add_(out.x, q_from_double_(p->lx));
    out.y = q_add_(out.y, q_from_double_(p->ly));
    out.z = q_add_(out.z, q_from_double_(p->lz));
    return out;
}

K26AstroPos k26astro_pos_from_fx(const K26AstroPosFx *p)
{
    K26AstroPos out = k26astro_pos_zero();
    if (!p) return out;
    out.lx = q_to_double_(p->x);
    out.ly = q_to_double_(p->y);
    out.lz = q_to_double_(p->z);
    k26astro_pos_normalise(&out);
    return out;
}
