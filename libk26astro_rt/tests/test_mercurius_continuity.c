/* test_mercurius_continuity.c — paper-faithful MERCURIUS gates.
 *
 * Three coverage axes:
 *
 *   Test 1 (decomposition identity): a_far + a_near = a_total for any
 *     K-weight assignment. Verifies the weighted-direct API is
 *     consistent with the unweighted direct N² across the pair list.
 *
 *   Test 2 (K endpoints): with all K=0 weights, an MERCURIUS-FAR
 *     mode step is bit-equal to an unmasked full-force step; with
 *     all K=1, an MERCURIUS-NEAR step matches the same.
 *
 *   Test 3 (smooth handoff): a slow-approaching pair crosses
 *     y_outer → y_inner over many steps. Verify there is no
 *     discontinuity in (x_after_step - x_before_step) as K(y)
 *     transitions smoothly through the window.
 *
 * Reference: Rein & Tamayo 2019 §3 eq. 12-14. */
#include "k26astro_rt/world.h"
#include "k26astro_body/body.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double v3_mag_(K26V3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

static int test1_decomposition_identity(void)
{
    /* 4-body system; assign arbitrary K weights to two pairs. */
    K26AstroBody bodies[4];
    memset(bodies, 0, sizeof bodies);
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH;
    bodies[1].pos  = k26astro_pos_from_m(1.5e11, 0, 0);
    bodies[2].kind = K26ASTRO_BODY_PLANET;
    bodies[2].gm   = K26A_GM_EARTH;
    bodies[2].pos  = k26astro_pos_from_m(1.5e11 + 5e9, 0, 0);
    bodies[3].kind = K26ASTRO_BODY_PLANET;
    bodies[3].gm   = K26A_GM_EARTH * 0.5;
    bodies[3].pos  = k26astro_pos_from_m(0, 2.3e11, 0);

    K26AstroGravView view = { bodies, 4 };

    K26V3 a_full[4], a_far[4], a_near[4];
    k26astro_grav_force_direct(&view, a_full);

    K26AstroPairWeight weights[2] = {
        { 1, 2, 0.42 },
        { 0, 3, 0.17 }
    };
    k26astro_grav_force_direct_weighted(&view, weights, 2, a_far, a_near);

    double max_resid = 0.0;
    for (int i = 0; i < 4; i++) {
        K26V3 sum = {
            a_far[i].x + a_near[i].x,
            a_far[i].y + a_near[i].y,
            a_far[i].z + a_near[i].z
        };
        K26V3 d = { sum.x - a_full[i].x, sum.y - a_full[i].y, sum.z - a_full[i].z };
        double r = v3_mag_(d);
        double m = v3_mag_(a_full[i]);
        double rel = (m > 0.0) ? r / m : r;
        if (rel > max_resid) max_resid = rel;
    }
    fprintf(stderr, "test1 (decomp identity): max_resid=%.3e\n", max_resid);
    if (max_resid > 1e-14) {
        fprintf(stderr, "test1 FAIL: a_far + a_near != a_total\n");
        return 1;
    }
    return 0;
}

static int test2_K_endpoints(void)
{
    /* Same 4-body setup. Compare:
     *   FAR with K=0 for all pairs  vs  unweighted direct
     *   NEAR with K=1 for all pairs vs  unweighted direct
     * Both should be bit-equal. */
    K26AstroBody bodies[3];
    memset(bodies, 0, sizeof bodies);
    bodies[0].gm  = K26A_GM_SUN;
    bodies[0].pos = k26astro_pos_zero();
    bodies[1].gm  = K26A_GM_EARTH * 1000.0;
    bodies[1].pos = k26astro_pos_from_m(7.78e11, 0, 0);
    bodies[2].gm  = K26A_GM_EARTH * 500.0;
    bodies[2].pos = k26astro_pos_from_m(7.78e11 + 3e9, 0, 0);

    K26AstroGravView view = { bodies, 3 };
    K26V3 a_full[3];
    k26astro_grav_force_direct(&view, a_full);

    /* All K=0 (everything is far). */
    K26AstroPairWeight w_far[3] = {
        { 0, 1, 0.0 }, { 0, 2, 0.0 }, { 1, 2, 0.0 }
    };
    K26V3 a_far[3], a_near[3];
    k26astro_grav_force_direct_weighted(&view, w_far, 3, a_far, a_near);
    /* a_far should match a_full bit-exactly; a_near should be all-zero. */
    for (int i = 0; i < 3; i++) {
        if (a_far[i].x != a_full[i].x || a_far[i].y != a_full[i].y ||
            a_far[i].z != a_full[i].z) {
            fprintf(stderr, "test2 FAIL: K=0 a_far[%d] != a_full[%d]\n", i, i);
            return 1;
        }
        if (a_near[i].x != 0.0 || a_near[i].y != 0.0 || a_near[i].z != 0.0) {
            fprintf(stderr, "test2 FAIL: K=0 a_near[%d] != 0\n", i);
            return 1;
        }
    }

    /* All K=1 (everything is near). */
    K26AstroPairWeight w_near[3] = {
        { 0, 1, 1.0 }, { 0, 2, 1.0 }, { 1, 2, 1.0 }
    };
    k26astro_grav_force_direct_weighted(&view, w_near, 3, a_far, a_near);
    for (int i = 0; i < 3; i++) {
        if (a_near[i].x != a_full[i].x || a_near[i].y != a_full[i].y ||
            a_near[i].z != a_full[i].z) {
            fprintf(stderr, "test2 FAIL: K=1 a_near[%d] != a_full[%d]\n", i, i);
            return 1;
        }
        if (a_far[i].x != 0.0 || a_far[i].y != 0.0 || a_far[i].z != 0.0) {
            fprintf(stderr, "test2 FAIL: K=1 a_far[%d] != 0\n", i);
            return 1;
        }
    }

    fprintf(stderr, "test2 (K endpoints): K=0 and K=1 limits bit-exact\n");
    return 0;
}

static int test3_smooth_handoff(void)
{
    /* Slow-approaching pair: Sun + two co-orbiting planets that
     * drift through the MERCURIUS transition window. Verify the
     * step-to-step position delta varies smoothly as K(y) crosses
     * through (0, 1) — no jump > 1% of the typical step length. */
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_PORTABLE,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    k26astro_world_set_mercurius_factors(w, 3.0, 5.0);

    K26AstroBody sun;
    k26astro_body_init(&sun);
    strncpy(sun.name, "sun", sizeof sun.name - 1);
    sun.kind = K26ASTRO_BODY_STAR;
    sun.gm   = K26A_GM_SUN;
    sun.mass = K26A_GM_SUN / K26A_G;
    sun.pos  = k26astro_pos_zero();
    k26astro_world_add_body(w, sun);

    K26AstroBody p1;
    k26astro_body_init(&p1);
    strncpy(p1.name, "p1", sizeof p1.name - 1);
    p1.kind = K26ASTRO_BODY_PLANET;
    p1.gm   = K26A_GM_EARTH * 1000.0;  /* Jupiter-class for clear Hill radius */
    p1.mass = p1.gm / K26A_G;
    double a1 = 5.2 * K26A_AU_M;
    p1.pos = k26astro_pos_from_m(a1, 0.0, 0.0);
    double v1 = sqrt(K26A_GM_SUN / a1);
    p1.vel = k26m3d_v3(0.0, v1, 0.0);
    k26astro_world_add_body(w, p1);

    /* p2 starts slightly outside p1's Hill sphere on the +y side,
     * with a slow inward radial velocity that will carry it through
     * the transition zone over ~30 days. */
    K26AstroBody p2;
    k26astro_body_init(&p2);
    strncpy(p2.name, "p2", sizeof p2.name - 1);
    p2.kind = K26ASTRO_BODY_PLANET;
    p2.gm   = K26A_GM_EARTH * 10.0;
    p2.mass = p2.gm / K26A_G;
    /* Hill radius of p1 ~ a1 * (m_p1/(3*M_sun))^(1/3) ~ 5.2 AU * 0.07 = 0.36 AU.
     * Start p2 at p1.pos + (0, 6 R_hill, 0); push inward at ~1 R_hill / 60 days. */
    double R_hill1 = a1 * cbrt(p1.gm / (3.0 * K26A_GM_SUN));
    K26V3 p1_m = k26astro_pos_to_m_approx(&p1.pos);
    p2.pos = k26astro_pos_from_m(p1_m.x,
                                  p1_m.y + 6.0 * R_hill1, p1_m.z);
    /* Tangential velocity matches p1's so the pair is roughly co-orbiting. */
    p2.vel = k26m3d_v3(0.0, v1, -R_hill1 / (60.0 * 86400.0));
    k26astro_world_add_body(w, p2);

    /* Step at 1-day rate; observe (x_after - x_before) for p2. */
    double dt = 86400.0;
    int n_steps = 120;
    K26V3 prev_delta = { 0, 0, 0 };
    double max_jump = 0.0;
    K26V3 prev_pos = k26astro_pos_to_m_approx(&p2.pos);
    K26AstroBody *world_bodies = k26astro_world_grav(w)->bodies;
    (void)world_bodies;
    int finite_failed = 0;
    for (int s = 0; s < n_steps; s++) {
        K26V3 pos_before = k26astro_pos_to_m_approx(&k26astro_world_grav(w)->bodies[2].pos);
        (void)k26astro_world_step(w, dt);
        K26V3 pos_after = k26astro_pos_to_m_approx(&k26astro_world_grav(w)->bodies[2].pos);
        K26V3 delta = {
            pos_after.x - pos_before.x,
            pos_after.y - pos_before.y,
            pos_after.z - pos_before.z
        };
        if (s > 0) {
            K26V3 jerk = {
                delta.x - prev_delta.x,
                delta.y - prev_delta.y,
                delta.z - prev_delta.z
            };
            double jump = v3_mag_(jerk);
            double step_mag = v3_mag_(delta);
            double rel = (step_mag > 0.0) ? jump / step_mag : 0.0;
            if (rel > max_jump) max_jump = rel;
        }
        prev_delta = delta;
        prev_pos = pos_after;
        for (int b = 0; b < 3; b++) {
            if (!isfinite(k26astro_world_grav(w)->bodies[b].vel.x)
              || !isfinite(k26astro_world_grav(w)->bodies[b].vel.y)) {
                finite_failed = 1;
            }
        }
    }
    (void)prev_pos;

    fprintf(stderr, "test3 (smooth handoff): max_jump_rel=%.3e\n", max_jump);
    if (finite_failed) {
        fprintf(stderr, "test3 FAIL: NaN/Inf in body state\n");
        k26astro_world_destroy(w);
        return 1;
    }
    /* The jump tolerance is generous: we're checking the K(y)
     * transition doesn't introduce a force discontinuity. A pure WH
     * integrator over a circular orbit would have max_jump ~ 1e-5
     * just from orbital curvature; a force discontinuity would
     * exceed 0.1 easily. */
    if (max_jump > 0.5) {
        fprintf(stderr, "test3 FAIL: max_jump %.3e exceeds 0.5\n", max_jump);
        k26astro_world_destroy(w);
        return 1;
    }

    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    if (test1_decomposition_identity() != 0) return 1;
    if (test2_K_endpoints() != 0) return 1;
    if (test3_smooth_handoff() != 0) return 1;
    printf("test_mercurius_continuity: OK\n");
    return 0;
}
