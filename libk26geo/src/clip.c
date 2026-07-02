/* Sutherland-Hodgman polygon clipping against an axis-aligned bbox.
 *
 * Standard SH: walk the polygon once per clip edge, accumulating a
 * new vertex list. After 4 passes (left/right/top/bottom of the bbox),
 * the polygon is clipped to the rect.
 *
 * The bbox is in (x, z) world metres; orientation matches our
 * Mercator projection (north is -Z, so min_z is the north edge). */

#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

/* "Inside" predicate for each of the 4 clip edges. */
typedef enum { EDGE_LEFT, EDGE_RIGHT, EDGE_TOP, EDGE_BOTTOM } ClipEdge;

static int inside_edge(double x, double z, ClipEdge e,
                       double min_x, double min_z,
                       double max_x, double max_z)
{
    switch (e) {
    case EDGE_LEFT:   return x >= min_x;
    case EDGE_RIGHT:  return x <= max_x;
    case EDGE_TOP:    return z >= min_z;
    case EDGE_BOTTOM: return z <= max_z;
    }
    return 1;
}

/* Intersect segment (p0, p1) with the clip line. Caller has already
 * confirmed the two endpoints straddle the line (one inside, one
 * outside). */
static void intersect_edge(double x0, double z0, double x1, double z1,
                           ClipEdge e,
                           double min_x, double min_z,
                           double max_x, double max_z,
                           double *out_x, double *out_z)
{
    double t = 0.0;
    switch (e) {
    case EDGE_LEFT:
        t = (min_x - x0) / (x1 - x0);
        *out_x = min_x;
        *out_z = z0 + t * (z1 - z0);
        break;
    case EDGE_RIGHT:
        t = (max_x - x0) / (x1 - x0);
        *out_x = max_x;
        *out_z = z0 + t * (z1 - z0);
        break;
    case EDGE_TOP:
        t = (min_z - z0) / (z1 - z0);
        *out_x = x0 + t * (x1 - x0);
        *out_z = min_z;
        break;
    case EDGE_BOTTOM:
        t = (max_z - z0) / (z1 - z0);
        *out_x = x0 + t * (x1 - x0);
        *out_z = max_z;
        break;
    }
}

size_t k26geo_clip_polygon_xz(const double *in_coords, size_t n_in,
                              double *out_coords, size_t cap_out,
                              double min_x, double min_z,
                              double max_x, double max_z)
{
    if (!in_coords || !out_coords || n_in < 3) return 0;
    if (cap_out < n_in * 2) return 0;     /* sanity floor */

    /* Worst-case SH growth: 2× per edge × 4 edges = 16× input verts.
     * In practice it grows much less. We bound the work buffer at
     * 32 × n_in doubles which covers any well-formed polygon. */
    const size_t cap = n_in * 32 + 64;
    double *bufA = (double *)malloc(cap * sizeof *bufA);
    double *bufB = (double *)malloc(cap * sizeof *bufB);
    if (!bufA || !bufB) { free(bufA); free(bufB); return 0; }

    memcpy(bufA, in_coords, n_in * 2 * sizeof(double));
    size_t cur_n = n_in;

    const ClipEdge edges[4] = { EDGE_LEFT, EDGE_RIGHT, EDGE_TOP, EDGE_BOTTOM };
    double *src = bufA;
    double *dst = bufB;

    for (int ei = 0; ei < 4; ei++) {
        const ClipEdge e = edges[ei];
        size_t out_n = 0;
        if (cur_n == 0) break;

        for (size_t i = 0; i < cur_n; i++) {
            const size_t j = (i + cur_n - 1) % cur_n;   /* previous vertex */
            const double xj = src[2 * j];
            const double zj = src[2 * j + 1];
            const double xi = src[2 * i];
            const double zi = src[2 * i + 1];
            const int curIn  = inside_edge(xi, zi, e, min_x, min_z, max_x, max_z);
            const int prevIn = inside_edge(xj, zj, e, min_x, min_z, max_x, max_z);

            if (curIn) {
                if (!prevIn) {
                    /* entering — emit intersection then current */
                    double ix, iz;
                    intersect_edge(xj, zj, xi, zi, e,
                                   min_x, min_z, max_x, max_z, &ix, &iz);
                    if (out_n * 2 + 4 > cap) goto bail;
                    dst[2 * out_n]     = ix; dst[2 * out_n + 1] = iz; out_n++;
                }
                if (out_n * 2 + 2 > cap) goto bail;
                dst[2 * out_n] = xi; dst[2 * out_n + 1] = zi; out_n++;
            } else if (prevIn) {
                /* leaving — emit intersection only */
                double ix, iz;
                intersect_edge(xj, zj, xi, zi, e,
                               min_x, min_z, max_x, max_z, &ix, &iz);
                if (out_n * 2 + 2 > cap) goto bail;
                dst[2 * out_n] = ix; dst[2 * out_n + 1] = iz; out_n++;
            }
            /* both outside — emit nothing */
        }

        cur_n = out_n;
        double *tmp = src; src = dst; dst = tmp;
    }

    /* `src` now holds the final clipped polygon. */
    if (cur_n < 3 || cur_n * 2 > cap_out) { free(bufA); free(bufB); return 0; }
    memcpy(out_coords, src, cur_n * 2 * sizeof(double));
    free(bufA); free(bufB);
    return cur_n;

bail:
    free(bufA); free(bufB);
    return 0;
}

/* ================================================================
 * Polyline clipping — Liang-Barsky per segment.
 *
 * For each segment of the input polyline, compute the parametric
 * (t_enter, t_exit) range that lies inside the bbox. If the segment
 * intersects the bbox, accumulate the clipped portion into the
 * current sub-polyline; emit + reset when we exit. Segments fully
 * inside extend the current sub-polyline by one vertex (the segment
 * end). Segments fully outside flush any open sub-polyline.
 * ================================================================ */

static int liang_barsky(double x0, double z0, double x1, double z1,
                        double min_x, double min_z,
                        double max_x, double max_z,
                        double *out_t_enter, double *out_t_exit)
{
    double dx = x1 - x0, dz = z1 - z0;
    double t_enter = 0.0, t_exit = 1.0;
    double p[4] = { -dx, +dx, -dz, +dz };
    double q[4] = { x0 - min_x, max_x - x0, z0 - min_z, max_z - z0 };

    for (int i = 0; i < 4; i++) {
        if (p[i] == 0.0) {
            if (q[i] < 0.0) return 0;  /* parallel and outside */
        } else {
            double t = q[i] / p[i];
            if (p[i] < 0.0) {
                if (t > t_enter) t_enter = t;
            } else {
                if (t < t_exit) t_exit = t;
            }
        }
    }
    if (t_enter > t_exit) return 0;
    *out_t_enter = t_enter;
    *out_t_exit  = t_exit;
    return 1;
}

static void poly_push(double *buf, size_t *buf_n, double x, double z)
{
    if (*buf_n > 0) {
        const double lx = buf[2 * (*buf_n - 1)];
        const double lz = buf[2 * (*buf_n - 1) + 1];
        if (lx == x && lz == z) return;
    }
    buf[2 * (*buf_n)]     = x;
    buf[2 * (*buf_n) + 1] = z;
    (*buf_n)++;
}

static void poly_flush(double *buf, size_t *buf_n,
                       K26GeoPolylineCb cb, void *user,
                       size_t *emitted)
{
    if (*buf_n >= 2) {
        cb(buf, *buf_n, user);
        (*emitted)++;
    }
    *buf_n = 0;
}

size_t k26geo_clip_polyline_xz(const double *in_coords, size_t n_in,
                               double min_x, double min_z,
                               double max_x, double max_z,
                               K26GeoPolylineCb cb, void *user)
{
    if (!in_coords || !cb || n_in < 2) return 0;

    /* Scratch buffer for the current open sub-polyline. Worst-case
     * length = 2·n_in (each in-bbox segment may add an entry-point
     * vertex on top of the kept endpoint). */
    double *buf = (double *)malloc(n_in * 4 * sizeof *buf);
    if (!buf) return 0;
    size_t buf_n = 0;
    size_t emitted = 0;

    for (size_t i = 0; i + 1 < n_in; i++) {
        const double x0 = in_coords[2 * i],         z0 = in_coords[2 * i + 1];
        const double x1 = in_coords[2 * (i + 1)],   z1 = in_coords[2 * (i + 1) + 1];
        double te = 0.0, tx = 1.0;
        if (!liang_barsky(x0, z0, x1, z1,
                          min_x, min_z, max_x, max_z, &te, &tx)) {
            poly_flush(buf, &buf_n, cb, user, &emitted);
            continue;
        }
        const double cx0 = x0 + te * (x1 - x0);
        const double cz0 = z0 + te * (z1 - z0);
        const double cx1 = x0 + tx * (x1 - x0);
        const double cz1 = z0 + tx * (z1 - z0);

        if (te > 0.0) {
            /* entering on this segment — start a fresh run at the
             * clip-line intersection. */
            poly_flush(buf, &buf_n, cb, user, &emitted);
            poly_push(buf, &buf_n, cx0, cz0);
        } else if (buf_n == 0) {
            poly_push(buf, &buf_n, cx0, cz0);
        }
        poly_push(buf, &buf_n, cx1, cz1);

        if (tx < 1.0) {
            /* exiting on this segment — emit and reset */
            poly_flush(buf, &buf_n, cb, user, &emitted);
        }
    }
    poly_flush(buf, &buf_n, cb, user, &emitted);
    free(buf);
    return emitted;
}
