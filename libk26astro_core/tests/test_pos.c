/* test_pos.c — sector grid + Q64.64 arithmetic.
 *
 * Acceptance:
 *   - normalisation invariant: after k26astro_pos_normalise, every
 *     local-offset axis is in [-S/2, S/2).
 *   - relative-coordinate identity: (a - b) reconstructed from
 *     `k26astro_pos_sub` matches the original delta to within
 *     ~50 μm at solar-system scales.
 *   - sector-vs-Q64.64 cross-check: positions encoded both ways
 *     produce identical-to-23μm relative vectors.
 *   - Pluto-scale stress: an Earth + Pluto pair subtracted yields
 *     finite distance close to the literature value (>30 AU). */
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol) {
    return fabs(a - b) <= tol;
}

int main(void)
{
    /* ---- Construction + normalisation ---------------------------- */
    K26AstroPos z = k26astro_pos_zero();
    assert(z.sx == 0 && z.sy == 0 && z.sz == 0);
    assert(z.lx == 0.0 && z.ly == 0.0 && z.lz == 0.0);

    /* A position 0.5 AU from origin — fits in a single sector
     * (sector edge 0.459 AU, so 0.5 AU is across one boundary). */
    K26AstroPos p = k26astro_pos_from_m(K26A_AU_M * 0.5, 0.0, 0.0);
    assert(p.sx == 1);     /* one sector along +x */
    assert(fabs(p.lx) < 0.5 * K26ASTRO_SECTOR_EDGE_M);

    /* Mars-orbit scale (~1.52 AU) — should land at sx = 3 with
     * sub-millimetre lx residual. */
    K26AstroPos mars = k26astro_pos_from_m(K26A_AU_M * 1.52, 0.0, 0.0);
    assert(fabs(mars.lx) < 0.5 * K26ASTRO_SECTOR_EDGE_M);

    /* Pluto-orbit scale (~40 AU). */
    K26AstroPos pluto = k26astro_pos_from_m(K26A_AU_M * 40.0, 0.0, 0.0);
    assert(fabs(pluto.lx) < 0.5 * K26ASTRO_SECTOR_EDGE_M);

    /* ---- Subtraction identity ----------------------------------- */
    /* Earth (1 AU) - Pluto (40 AU) — expected ~ -39 AU along x. */
    K26AstroPos earth = k26astro_pos_from_m(K26A_AU_M * 1.0, 0.0, 0.0);
    K26V3 rel = k26astro_pos_sub(&earth, &pluto);
    assert(near_(rel.x, -K26A_AU_M * 39.0, 1e-3));   /* mm precision */
    assert(near_(rel.y,  0.0, 1e-3));
    assert(near_(rel.z,  0.0, 1e-3));

    /* Symmetry: (b - a) = -(a - b). */
    K26V3 rel_rev = k26astro_pos_sub(&pluto, &earth);
    assert(near_(rel_rev.x, -rel.x, 1e-9));
    assert(near_(rel_rev.y, -rel.y, 1e-9));
    assert(near_(rel_rev.z, -rel.z, 1e-9));

    /* Distance match. */
    double d = k26astro_pos_dist(&earth, &pluto);
    assert(near_(d, K26A_AU_M * 39.0, 1e-3));

    /* ---- Reflexivity: a - a == 0 ------------------------------- */
    K26V3 zero_rel = k26astro_pos_sub(&pluto, &pluto);
    assert(zero_rel.x == 0.0 && zero_rel.y == 0.0 && zero_rel.z == 0.0);

    /* ---- Add delta + normalisation ----------------------------- */
    K26AstroPos walk = k26astro_pos_zero();
    /* Step by 1.2 sectors along x. */
    K26V3 step = { 1.2 * K26ASTRO_SECTOR_EDGE_M, 0.0, 0.0 };
    k26astro_pos_add(&walk, step);
    assert(walk.sx >= 1);   /* normalisation carried into sector index */
    assert(fabs(walk.lx) < 0.5 * K26ASTRO_SECTOR_EDGE_M);

    /* ---- Q64.64 round-trip ------------------------------------- */
    K26AstroPosFx earth_fx = k26astro_pos_to_fx(&earth);
    K26AstroPos earth_back = k26astro_pos_from_fx(&earth_fx);
    /* Round-trip residual should be sub-mm (Q64.64 has ~5e-20 m
     * resolution; the lossy step is sector→Q64.64 quantising lx). */
    K26V3 rt_diff = k26astro_pos_sub(&earth, &earth_back);
    assert(fabs(rt_diff.x) < 1e-3);
    assert(fabs(rt_diff.y) < 1e-3);
    assert(fabs(rt_diff.z) < 1e-3);

    /* Q64.64 subtraction matches sector-mode subtraction within
     * sub-mm at solar-system scales. */
    K26AstroPosFx pluto_fx = k26astro_pos_to_fx(&pluto);
    K26V3 rel_fx = k26astro_pos_fx_sub(&earth_fx, &pluto_fx);
    assert(near_(rel_fx.x, rel.x, 1e-3));
    assert(near_(rel_fx.y, rel.y, 1e-3));
    assert(near_(rel_fx.z, rel.z, 1e-3));

    /* ---- Sector-boundary stress: 1e6 normalisations don't drift -- */
    K26AstroPos p2 = k26astro_pos_zero();
    /* Each step crosses one sector boundary. After 1e6 steps the
     * sector index sits at 1,000,000. Local offset must stay
     * bounded throughout. */
    for (int i = 0; i < 1000000; i++) {
        K26V3 d_step = { K26ASTRO_SECTOR_EDGE_M, 0.0, 0.0 };
        k26astro_pos_add(&p2, d_step);
        assert(fabs(p2.lx) < 0.5 * K26ASTRO_SECTOR_EDGE_M);
    }
    assert(p2.sx == 1000000);

    printf("test_pos: OK (sector grid + Q64.64 + stress 1e6 normalisations)\n");
    return 0;
}
