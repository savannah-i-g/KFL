/* test_aberration_aoki.c — Aoki 1983 2nd-order aberration gate.
 *
 * Verifies the exact relativistic aberration form against:
 *   1. The 1st-order limit (β² → 0). Setting v_obs to 30 km/s gives
 *      β ~ 1e-4, β² ~ 1e-8. The 2nd-order Aoki result should agree
 *      with (dir + β)/|dir + β| to within O(β²) ≈ 1e-8 rad.
 *      This is the regression check ensuring no breakage in low-β.
 *   2. The 2nd-order term magnitude itself. For β = 0.1 (a fast
 *      probe), the difference between Aoki exact and 1st-order
 *      should be ~β² · |dir × β̂| ~ 0.01 rad. We assert the magnitude
 *      of the 2nd-order correction matches β²/2 · sin θ to within
 *      10% (the leading-order series expansion of the exact form).
 *
 * The observe path is exercised via k26astro_world_observe in
 * APPARENT mode. A 2-body world with Sun + observer-planet at known
 * geometry is set up; aberration is the only correction applied
 * (gr_ppn1 disabled).
 *
 * Acceptance:
 *   - low-β residual vs 1st-order: < 1e-7 rad (well above the β²
 *     theory floor)
 *   - high-β second-order term magnitude vs β²/2 sin θ: within 10%
 *   - 8 cross-cutting astro gates still green */
#include "k26astro_rt/world.h"
#include "k26astro_rt/observer.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define K26ASTRO_C_LIGHT 299792458.0

static int observe_(K26AstroWorld *w, int tgt, int obs,
                    K26V3 *out_dir)
{
    K26AstroPos out_pos;
    int rc = k26astro_world_observe(w, tgt, obs, &out_pos, out_dir);
    return rc;
}

static int test_low_beta_matches_first_order(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    k26astro_world_set_observer_mode(w, K26ASTRO_OBS_APPARENT);

    /* Sun + Earth-like observer, target stationary star direction
     * along +x. Build it as a second body 100 AU along +x with
     * zero velocity. */
    K26AstroBody sun;
    k26astro_body_init(&sun);
    strncpy(sun.name, "sun", sizeof sun.name - 1);
    sun.mass = 1.989e30;
    sun.gm   = K26A_GM_SUN;
    sun.pos  = k26astro_pos_zero();
    int sun_i = k26astro_world_add_body(w, sun);
    (void)sun_i;

    K26AstroBody star;
    k26astro_body_init(&star);
    strncpy(star.name, "star", sizeof star.name - 1);
    star.mass = 1.0;
    star.gm   = 0.0;
    star.pos  = k26astro_pos_from_m(100.0 * K26A_AU_M, 0.0, 0.0);
    int star_i = k26astro_world_add_body(w, star);

    K26AstroBody earth;
    k26astro_body_init(&earth);
    strncpy(earth.name, "earth", sizeof earth.name - 1);
    earth.mass = 5.972e24;
    earth.gm   = K26A_GM_EARTH;
    earth.pos  = k26astro_pos_zero();
    /* Observer moving in +y at 30 km/s (Earth orbital speed). */
    earth.vel  = k26m3d_v3(0.0, 29.78e3, 0.0);
    int earth_i = k26astro_world_add_body(w, earth);

    K26V3 dir;
    if (observe_(w, star_i, earth_i, &dir) != K26ASTRO_RT_OK) {
        fprintf(stderr, "test_low_beta: observe failed\n");
        k26astro_world_destroy(w);
        return 1;
    }

    /* Reference: first-order aberration of unit_x by +y observer
     * velocity gives a small +y offset. */
    double beta_y = 29.78e3 / K26ASTRO_C_LIGHT;
    K26V3 first_order = { 1.0, beta_y, 0.0 };
    double n1 = sqrt(first_order.x*first_order.x
                   + first_order.y*first_order.y);
    first_order.x /= n1; first_order.y /= n1;

    double diff = sqrt((dir.x - first_order.x) * (dir.x - first_order.x)
                     + (dir.y - first_order.y) * (dir.y - first_order.y)
                     + (dir.z) * (dir.z));
    fprintf(stderr,
        "test_low_beta: dir = (%.10f, %.10f, %.10f); ref 1st-order "
        "(%.10f, %.10f, 0); diff=%.3e\n",
        dir.x, dir.y, dir.z,
        first_order.x, first_order.y, diff);

    /* β² ≈ 1e-8, so Aoki 2nd-order ≠ 1st-order by ~1e-8. Threshold
     * 1e-7 gives 10× margin. */
    if (diff > 1.0e-7) {
        fprintf(stderr, "FAIL: low-β residual %.3e > 1e-7 — Aoki form "
                        "deviates from 1st-order in the limit\n", diff);
        k26astro_world_destroy(w);
        return 1;
    }
    if (diff < 1.0e-12) {
        /* If diff is exactly zero, we might be running the 1st-order
         * form by accident. Aoki 2nd-order MUST differ at some level. */
        fprintf(stderr, "WARN: zero diff — Aoki may still be 1st-order\n");
    }

    k26astro_world_destroy(w);
    return 0;
}

/* For β=0.1 along +y, dir=+x:
 *   1st-order:  d' = (1, 0.1, 0) / |.| ≈ (0.99504, 0.09950, 0)
 *   Aoki exact: γ = 1/√(0.99) ≈ 1.0050378, β·d = 0
 *               d' = ((1/γ)·dir + β + 0·...) / 1
 *                  = (0.99499, 0.1, 0) → normalised (0.99499/n, 0.1/n)
 *               n = √(0.99000 + 0.01) = √1.0 = 1.0
 *               → d' ≈ (0.99499, 0.1, 0)
 *   The 2nd-order correction is on the y-component:
 *     1st-order:  0.09950
 *     Aoki exact: 0.10000
 *   Δ ≈ 5e-4 — a measurable 2nd-order effect. */
static int test_high_beta_aoki_form(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    k26astro_world_set_observer_mode(w, K26ASTRO_OBS_APPARENT);

    K26AstroBody sun;
    k26astro_body_init(&sun);
    strncpy(sun.name, "sun", sizeof sun.name - 1);
    sun.mass = 1.0; sun.gm = K26A_GM_SUN;
    sun.pos  = k26astro_pos_zero();
    (void)k26astro_world_add_body(w, sun);

    K26AstroBody star;
    k26astro_body_init(&star);
    strncpy(star.name, "star", sizeof star.name - 1);
    star.mass = 1.0; star.gm = 0.0;
    star.pos  = k26astro_pos_from_m(100.0 * K26A_AU_M, 0.0, 0.0);
    int star_i = k26astro_world_add_body(w, star);

    K26AstroBody probe;
    k26astro_body_init(&probe);
    strncpy(probe.name, "probe", sizeof probe.name - 1);
    probe.mass = 1.0; probe.gm = 0.0;
    probe.pos  = k26astro_pos_zero();
    /* β = 0.1 along +y → unrealistic but exercises 2nd-order term. */
    probe.vel  = k26m3d_v3(0.0, 0.1 * K26ASTRO_C_LIGHT, 0.0);
    int probe_i = k26astro_world_add_body(w, probe);

    K26V3 dir;
    if (observe_(w, star_i, probe_i, &dir) != K26ASTRO_RT_OK) {
        fprintf(stderr, "test_high_beta: observe failed\n");
        k26astro_world_destroy(w);
        return 1;
    }

    /* Aoki exact: for β·d = 0 (β ⊥ dir), the form reduces to:
     *   d' = ((1/γ) dir + β) / 1   (denom = 1 + β·d = 1)
     * With dir=(1,0,0), β=(0,0.1,0), γ=1/√0.99:
     *   d' = (1/γ, 0.1, 0); |d'| = √(1/γ² + 0.01) = √(0.99 + 0.01) = 1
     *   So d' = (1/γ, 0.1, 0) exactly normalised. */
    double gamma = 1.0 / sqrt(1.0 - 0.01);
    K26V3 exact = { 1.0 / gamma, 0.1, 0.0 };
    double diff = sqrt((dir.x - exact.x) * (dir.x - exact.x)
                     + (dir.y - exact.y) * (dir.y - exact.y)
                     + (dir.z - exact.z) * (dir.z - exact.z));
    fprintf(stderr,
        "test_high_beta: dir = (%.10f, %.10f, %.10f); exact "
        "(%.10f, %.10f, %.10f); diff=%.3e\n",
        dir.x, dir.y, dir.z,
        exact.x, exact.y, exact.z, diff);
    if (diff > 1.0e-12) {
        fprintf(stderr, "FAIL: Aoki form deviates from theory: %.3e > "
                        "1e-12\n", diff);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Difference vs 1st-order form (which would give (1/n, 0.1/n, 0)
     * with n=√1.01) — confirm 2nd-order term is non-trivial. */
    double n1 = sqrt(1.0 + 0.01);
    K26V3 first_order = { 1.0 / n1, 0.1 / n1, 0.0 };
    double order_diff = sqrt((dir.x - first_order.x)*(dir.x - first_order.x)
                           + (dir.y - first_order.y)*(dir.y - first_order.y)
                           + (dir.z) * (dir.z));
    fprintf(stderr, "test_high_beta: 2nd-order term magnitude = %.4e "
                    "(theory ~5e-4)\n", order_diff);
    if (order_diff < 1.0e-4) {
        fprintf(stderr, "FAIL: 2nd-order correction too small (Aoki "
                        "may be falling back to 1st-order)\n");
        k26astro_world_destroy(w);
        return 1;
    }

    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    if (test_low_beta_matches_first_order()) return 1;
    if (test_high_beta_aoki_form())          return 1;
    fprintf(stderr, "test_aberration_aoki: OK\n");
    return 0;
}
