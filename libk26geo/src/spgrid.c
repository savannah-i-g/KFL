/* Uniform G×G spatial grid.
 *
 * Backed by a per-cell dynamic int array. Adds are O(cells-touched);
 * a typical OSM polygon touches 1-2 cells at 250 m grid resolution.
 *
 * Memory ceiling: 256 cells × 32 ints starter = 32 KiB scratch on a
 * 1 km² grid; grows on demand with 1.5× expansion. */

#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

typedef struct {
    int    *ids;
    size_t  n;
    size_t  cap;
} K26GeoSpCell;

struct K26GeoSpGrid {
    double  min_x, min_z;
    double  cell_w, cell_h;
    int     nx,    nz;
    K26GeoSpCell *cells;   /* nx * nz */
};

K26GeoSpGrid *k26geo_spgrid_new(double min_x, double min_z,
                                double max_x, double max_z,
                                int n_cells_x, int n_cells_z)
{
    if (n_cells_x < 1 || n_cells_z < 1)             return NULL;
    if (!(max_x > min_x) || !(max_z > min_z))       return NULL;

    K26GeoSpGrid *g = (K26GeoSpGrid *)calloc(1, sizeof *g);
    if (!g) return NULL;
    g->min_x  = min_x;
    g->min_z  = min_z;
    g->nx     = n_cells_x;
    g->nz     = n_cells_z;
    g->cell_w = (max_x - min_x) / (double)n_cells_x;
    g->cell_h = (max_z - min_z) / (double)n_cells_z;
    g->cells  = (K26GeoSpCell *)calloc((size_t)n_cells_x * (size_t)n_cells_z,
                                       sizeof(K26GeoSpCell));
    if (!g->cells) { free(g); return NULL; }
    return g;
}

void k26geo_spgrid_free(K26GeoSpGrid *g)
{
    if (!g) return;
    if (g->cells) {
        const size_t n = (size_t)g->nx * (size_t)g->nz;
        for (size_t i = 0; i < n; i++) free(g->cells[i].ids);
        free(g->cells);
    }
    free(g);
}

int k26geo_spgrid_n_x(const K26GeoSpGrid *g) { return g ? g->nx : 0; }
int k26geo_spgrid_n_z(const K26GeoSpGrid *g) { return g ? g->nz : 0; }

static int clamp_cell(int v, int lo, int hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static K26GeoStatus cell_push(K26GeoSpCell *c, int id)
{
    if (c->n == c->cap) {
        size_t newcap = c->cap ? c->cap + (c->cap >> 1) + 1 : 8;
        int *p = (int *)realloc(c->ids, newcap * sizeof *p);
        if (!p) return K26GEO_ERR_OOM;
        c->ids = p;
        c->cap = newcap;
    }
    c->ids[c->n++] = id;
    return K26GEO_OK;
}

K26GeoStatus k26geo_spgrid_add(K26GeoSpGrid *g, int feature_id,
                               double minx, double minz,
                               double maxx, double maxz)
{
    if (!g)                  return K26GEO_ERR_INVAL;
    if (g->cell_w <= 0.0)    return K26GEO_ERR_INVAL;
    if (g->cell_h <= 0.0)    return K26GEO_ERR_INVAL;

    int cx0 = (int)((minx - g->min_x) / g->cell_w);
    int cx1 = (int)((maxx - g->min_x) / g->cell_w);
    int cz0 = (int)((minz - g->min_z) / g->cell_h);
    int cz1 = (int)((maxz - g->min_z) / g->cell_h);

    cx0 = clamp_cell(cx0, 0, g->nx - 1);
    cx1 = clamp_cell(cx1, 0, g->nx - 1);
    cz0 = clamp_cell(cz0, 0, g->nz - 1);
    cz1 = clamp_cell(cz1, 0, g->nz - 1);

    for (int cz = cz0; cz <= cz1; cz++) {
        for (int cx = cx0; cx <= cx1; cx++) {
            const size_t idx = (size_t)cz * (size_t)g->nx + (size_t)cx;
            K26GeoStatus s = cell_push(&g->cells[idx], feature_id);
            if (s != K26GEO_OK) return s;
        }
    }
    return K26GEO_OK;
}

void k26geo_spgrid_visit_rect(const K26GeoSpGrid *g,
                              double minx, double minz,
                              double maxx, double maxz,
                              K26GeoSpGridVisitFn fn, void *user)
{
    if (!g || !fn) return;

    int cx0 = (int)((minx - g->min_x) / g->cell_w);
    int cx1 = (int)((maxx - g->min_x) / g->cell_w);
    int cz0 = (int)((minz - g->min_z) / g->cell_h);
    int cz1 = (int)((maxz - g->min_z) / g->cell_h);

    cx0 = clamp_cell(cx0, 0, g->nx - 1);
    cx1 = clamp_cell(cx1, 0, g->nx - 1);
    cz0 = clamp_cell(cz0, 0, g->nz - 1);
    cz1 = clamp_cell(cz1, 0, g->nz - 1);

    for (int cz = cz0; cz <= cz1; cz++) {
        for (int cx = cx0; cx <= cx1; cx++) {
            const size_t idx = (size_t)cz * (size_t)g->nx + (size_t)cx;
            const K26GeoSpCell *c = &g->cells[idx];
            fn(cx, cz, c->ids, c->n, user);
        }
    }
}

void k26geo_spgrid_cell_bounds(const K26GeoSpGrid *g, int cx, int cz,
                               double *out_min_x, double *out_min_z,
                               double *out_max_x, double *out_max_z)
{
    if (!g) return;
    if (cx < 0 || cx >= g->nx || cz < 0 || cz >= g->nz) return;
    const double x0 = g->min_x + cx * g->cell_w;
    const double z0 = g->min_z + cz * g->cell_h;
    if (out_min_x) *out_min_x = x0;
    if (out_min_z) *out_min_z = z0;
    if (out_max_x) *out_max_x = x0 + g->cell_w;
    if (out_max_z) *out_max_z = z0 + g->cell_h;
}
