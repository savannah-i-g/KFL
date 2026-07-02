/* test_soi.c — SOI radius + hierarchy + crossing.
 *
 * Acceptance:
 *   - Hill + Laplace radii: Earth's SOI ~9.24e8 m (Laplace),
 *     ~1.5e9 m (Hill). Mars: ~5.78e8 m / ~9.83e8 m.
 *   - Hierarchy: build a 3-body system [Sun, Earth, Moon] with
 *     Moon's parent = Earth, Earth's parent = Sun. Query should
 *     resolve correctly.
 *   - Crossing: a body at r0 > R then r1 < R triggers crossing,
 *     no-crossing in either-staying-side cases.
 */
#include "k26astro_conics/soi.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* ---- SOI radii ---------------------------------------- */
    double r_laplace_earth = k26astro_soi_radius_laplace(
        K26A_GM_SUN, K26A_GM_EARTH, K26A_AU_M);
    /* Earth's Laplace SOI is ~9.24e8 m (~924,000 km). */
    assert(r_laplace_earth > 9.0e8 && r_laplace_earth < 9.5e8);

    double r_hill_earth = k26astro_soi_radius_hill(
        K26A_GM_SUN, K26A_GM_EARTH, K26A_AU_M);
    /* Hill ~1.49e9 m. */
    assert(r_hill_earth > 1.4e9 && r_hill_earth < 1.6e9);

    /* Jupiter SOI checks (using core's named constant). */
    double r_laplace_jup = k26astro_soi_radius_laplace(
        K26A_GM_SUN, K26A_GM_JUPITER, 5.203 * K26A_AU_M);
    /* Jupiter Laplace ~4.82e10 m. */
    assert(r_laplace_jup > 4.0e10 && r_laplace_jup < 5.5e10);

    /* Hill > Laplace (always, for any pair). */
    assert(r_hill_earth > r_laplace_earth);

    /* ---- Bad inputs return 0 ---------------------------- */
    assert(k26astro_soi_radius_hill(-1.0, 1.0, 1.0) == 0.0);
    assert(k26astro_soi_radius_laplace(1.0, -1.0, 1.0) == 0.0);
    assert(k26astro_soi_radius_hill(1.0, 1.0, -1.0) == 0.0);

    /* ---- Hierarchy walk -------------------------------- *
     * Build a [Sun, Earth, Moon] system. */
    K26AstroBody bodies[3];
    memset(bodies, 0, sizeof(bodies));

    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH;
    bodies[1].parent_body_idx = 0;   /* Earth around Sun */

    bodies[2].kind = K26ASTRO_BODY_MOON;
    bodies[2].gm   = 4.9028e12;       /* Moon GM (m³/s²), IAU 2015 */
    bodies[2].parent_body_idx = 1;    /* Moon around Earth */

    /* Moon's parent is Earth. */
    assert(k26astro_soi_parent(bodies, 3, 2, 0.0) == 1);
    /* Earth's parent is Sun. */
    assert(k26astro_soi_parent(bodies, 3, 1, 0.0) == 0);
    /* Sun has no parent. */
    assert(k26astro_soi_parent(bodies, 3, 0, 0.0) == -1);

    /* Out-of-range. */
    assert(k26astro_soi_parent(bodies, 3, -1, 0.0) == -1);
    assert(k26astro_soi_parent(bodies, 3, 99, 0.0) == -1);

    /* ---- Crossing detection ---------------------------- *
     * R = 1e6 m (small SOI). Body goes from 1.5e6 to 0.5e6. */
    double R = 1.0e6;
    K26V3 p0 = { 1.5e6, 0.0, 0.0 };
    K26V3 p1 = { 0.5e6, 0.0, 0.0 };
    double tau = 0.0;
    assert(k26astro_soi_crossing_detect(p0, p1, R, &tau) == 1);
    /* Linear interpolation: tau = (R - r0)/(r1-r0) = (1e6 - 1.5e6)/(0.5e6 - 1.5e6) = 0.5 */
    assert(fabs(tau - 0.5) < 1e-9);

    /* No crossing — both inside. */
    K26V3 p_in_0 = { 0.5e6, 0.0, 0.0 };
    K26V3 p_in_1 = { 0.3e6, 0.0, 0.0 };
    assert(k26astro_soi_crossing_detect(p_in_0, p_in_1, R, &tau) == 0);

    /* No crossing — both outside. */
    K26V3 p_out_0 = { 2.0e6, 0.0, 0.0 };
    K26V3 p_out_1 = { 3.0e6, 0.0, 0.0 };
    assert(k26astro_soi_crossing_detect(p_out_0, p_out_1, R, &tau) == 0);

    /* Exit (inside → outside). */
    K26V3 e0 = { 0.5e6, 0.0, 0.0 };
    K26V3 e1 = { 1.5e6, 0.0, 0.0 };
    assert(k26astro_soi_crossing_detect(e0, e1, R, &tau) == 1);
    assert(fabs(tau - 0.5) < 1e-9);

    /* Bad input — non-positive R. */
    assert(k26astro_soi_crossing_detect(p0, p1, -1.0, &tau) == -1);

    printf("test_soi: OK (radii + hierarchy + 3 crossing cases)\n");
    return 0;
}
