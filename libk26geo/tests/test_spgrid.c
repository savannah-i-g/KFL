#include <stdio.h>
#include <stdlib.h>

#include "k26geo.h"

static int fail = 0;

typedef struct {
    int n_visits;
    int n_total_ids;
    int saw_id[64];
} Visit;

static void visit(int cx, int cz, const int *ids, size_t n, void *user)
{
    (void)cx; (void)cz;
    Visit *v = (Visit *)user;
    v->n_visits++;
    v->n_total_ids += (int)n;
    for (size_t i = 0; i < n; i++) {
        if (ids[i] >= 0 && ids[i] < 64) v->saw_id[ids[i]] = 1;
    }
}

int main(void)
{
    /* Bad inputs */
    if (k26geo_spgrid_new(0, 0, 1, 1, 0, 1) != NULL) {
        fprintf(stderr, "FAIL zero cells should reject\n"); fail++;
    }
    if (k26geo_spgrid_new(0, 0, 0, 1, 4, 4) != NULL) {
        fprintf(stderr, "FAIL degenerate x extent should reject\n"); fail++;
    }

    /* Build a 10×10 grid over (0..1000, 0..1000). 100 m / cell. */
    K26GeoSpGrid *g = k26geo_spgrid_new(0, 0, 1000, 1000, 10, 10);
    if (!g) { fprintf(stderr, "FAIL alloc\n"); return 1; }

    if (k26geo_spgrid_n_x(g) != 10 || k26geo_spgrid_n_z(g) != 10) {
        fprintf(stderr, "FAIL dimensions\n"); fail++;
    }

    /* Add a point feature at (50, 50) — single cell. */
    k26geo_spgrid_add(g, 1, 50, 50, 50, 50);
    /* Add a 250 m square at (300..550, 300..550) — 3×3 cells (cells 3..5 inclusive). */
    k26geo_spgrid_add(g, 2, 300, 300, 550, 550);
    /* Add a sliver crossing a cell boundary (95..105, 0..1). */
    k26geo_spgrid_add(g, 3, 95, 0, 105, 1);

    /* Query [0..200, 0..200] — overlaps cells 0..2 inclusive on each
     * axis (200 lands exactly on the cell-1/cell-2 boundary; we include
     * cell 2 conservatively).  9 cells total. */
    {
        Visit v = {0};
        k26geo_spgrid_visit_rect(g, 0, 0, 200, 200, visit, &v);
        if (v.n_visits != 9) {
            fprintf(stderr, "FAIL visit count: got %d, want 9\n", v.n_visits); fail++;
        }
        if (!v.saw_id[1] || !v.saw_id[3]) {
            fprintf(stderr, "FAIL expected ids 1+3 in [0..200] window\n"); fail++;
        }
        if (v.saw_id[2]) {
            fprintf(stderr, "FAIL id 2 should be outside [0..200]\n"); fail++;
        }
    }

    /* Query around the 250m square. */
    {
        Visit v = {0};
        k26geo_spgrid_visit_rect(g, 250, 250, 600, 600, visit, &v);
        if (!v.saw_id[2]) {
            fprintf(stderr, "FAIL id 2 should be in 250..600 window\n"); fail++;
        }
    }

    /* Out-of-grid query clamps to bounds (should still visit valid cells). */
    {
        Visit v = {0};
        k26geo_spgrid_visit_rect(g, -1e6, -1e6, 1e6, 1e6, visit, &v);
        if (v.n_visits != 100) {
            fprintf(stderr, "FAIL clamped query visits: got %d, want 100\n", v.n_visits); fail++;
        }
    }

    /* Cell bounds. */
    {
        double a, b, c, d;
        k26geo_spgrid_cell_bounds(g, 3, 3, &a, &b, &c, &d);
        if (a != 300 || b != 300 || c != 400 || d != 400) {
            fprintf(stderr, "FAIL cell (3,3) bounds: (%g,%g,%g,%g)\n", a, b, c, d); fail++;
        }
    }

    k26geo_spgrid_free(g);

    if (fail) {
        fprintf(stderr, "%d failure(s)\n", fail);
        return 1;
    }
    puts("test_spgrid: OK");
    return 0;
}
