/* test_body_step.c — per-body Kepler advance gate.
 *
 * Validates that k26astro_world_body_step() advances a body along
 * its current Kepler orbit around its SOI parent without invoking
 * the full integrator and without touching other bodies. */
#include "k26astro_rt/world.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* IAU-2015 GM_SUN, m³/s². */
#define K26_GM_SUN 1.32712440018e20
#define AU         1.495978707e11

static double v3_norm_(K26V3 v)
{
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

int main(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    /* Sun at origin. */
    K26AstroBody sun;
    k26astro_body_init(&sun);
    strncpy(sun.name, "sun", sizeof sun.name - 1);
    sun.mass = 1.989e30; sun.gm = K26_GM_SUN;
    sun.pos = k26astro_pos_from_m(0.0, 0.0, 0.0);
    sun.vel = k26m3d_v3(0.0, 0.0, 0.0);
    int sun_i = k26astro_world_add_body(w, sun);

    /* Earth in a circular 1 AU orbit on +x → +y. v_circ = sqrt(GM/r). */
    double v_circ = sqrt(K26_GM_SUN / AU);
    K26AstroBody earth;
    k26astro_body_init(&earth);
    strncpy(earth.name, "earth", sizeof earth.name - 1);
    earth.mass = 5.972e24; earth.gm = 3.986004418e14;
    earth.pos = k26astro_pos_from_m(AU, 0.0, 0.0);
    earth.vel = k26m3d_v3(0.0, v_circ, 0.0);
    earth.parent_body_idx = sun_i;
    int earth_i = k26astro_world_add_body(w, earth);

    /* A second body that should NOT change: a "marker" at (2 AU, 0). */
    K26AstroBody mark;
    k26astro_body_init(&mark);
    strncpy(mark.name, "mark", sizeof mark.name - 1);
    mark.mass = 1.0; mark.gm = 0.0;
    mark.pos = k26astro_pos_from_m(2.0 * AU, 0.0, 0.0);
    mark.vel = k26m3d_v3(0.0, 0.0, 0.0);
    int mark_i = k26astro_world_add_body(w, mark);

    /* Snapshot mark's state. */
    K26AstroPos mark_pos_before = k26astro_world_body_at(w, mark_i)->pos;

    /* Quarter-period advance: T = 2π√(r³/μ). Quarter ≈ 90 days for
     * Earth's orbit. After this step Earth should be near (0, AU, 0). */
    double T = 2.0 * M_PI * sqrt(AU * AU * AU / K26_GM_SUN);
    int rc = k26astro_world_body_step(w, earth_i, T / 4.0);
    if (rc != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: body_step rc=%d\n", rc);
        k26astro_world_destroy(w);
        return 1;
    }

    K26AstroBody *e = k26astro_world_body_at(w, earth_i);
    K26V3 r_e = k26astro_pos_to_m_approx(&e->pos);
    K26V3 expected = { 0.0, AU, 0.0 };
    K26V3 d = { r_e.x - expected.x, r_e.y - expected.y, r_e.z - expected.z };
    double err = v3_norm_(d);
    fprintf(stderr,
        "body_step quarter-period: |r|=%.6e m, target=(0,1AU,0), miss=%.3e m\n",
        v3_norm_(r_e), err);
    /* Numerical universal-anomaly Newton-Raphson should give sub-mm
     * agreement on circular orbits. Allow generous 1 km tolerance. */
    if (err > 1.0e3) {
        fprintf(stderr, "FAIL: Kepler quarter-period miss too large\n");
        k26astro_world_destroy(w);
        return 1;
    }

    /* Verify marker body is unchanged. */
    K26AstroBody *m = k26astro_world_body_at(w, mark_i);
    K26V3 dmark = k26astro_pos_sub(&m->pos, &mark_pos_before);
    if (v3_norm_(dmark) != 0.0) {
        fprintf(stderr, "FAIL: marker body moved %.3e m (body_step "
                        "should only mutate target)\n", v3_norm_(dmark));
        k26astro_world_destroy(w);
        return 1;
    }

    /* Error cases. */
    if (k26astro_world_body_step(NULL, 0, 1.0) != -K26ASTRO_RT_E_NULL) {
        fprintf(stderr, "FAIL: NULL world should return E_NULL\n");
        return 1;
    }
    if (k26astro_world_body_step(w, -1, 1.0) != -K26ASTRO_RT_E_BAD_ARG) {
        fprintf(stderr, "FAIL: bad body_idx should return E_BAD_ARG\n");
        return 1;
    }
    if (k26astro_world_body_step(w, sun_i, 1.0) != -K26ASTRO_RT_E_BAD_ARG) {
        fprintf(stderr, "FAIL: Sun (parent==self) should return E_BAD_ARG\n");
        return 1;
    }

    k26astro_world_destroy(w);
    fprintf(stderr, "test_body_step: OK\n");
    return 0;
}
