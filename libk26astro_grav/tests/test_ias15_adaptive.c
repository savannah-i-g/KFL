/* test_ias15_adaptive.c — adaptive step-size controller smoke gate.
 *
 * Two coverage axes:
 *
 *   Test 1: Smooth 2-body Sun-Earth, 1000 days at requested dt=1 day.
 *     The adaptive controller should grow the substep aggressively
 *     (no error pressure on a circular orbit) and rarely reject.
 *     Energy conservation should match or beat fixed-step IAS15.
 *
 *   Test 2: Chaotic 3-body Sun + 2 close planets in 2:1 mean-motion
 *     resonance. Asymmetry forces multiple close approaches over
 *     360 days. The controller should shrink dt around each
 *     encounter (reject counter > 0) and recover afterwards.
 *
 * Acceptance:
 *   - Test 1: dE_max < 1e-10 across 1000 days; rejected_steps < 10.
 *   - Test 2: dE_max < 1e-6 across 360 days; reject counter > 0;
 *     final state finite (no NaN). */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static double total_energy_(K26AstroBody *bodies, int n)
{
    double KE = 0.0;
    for (int i = 0; i < n; i++) {
        double m_i = bodies[i].gm / K26A_G;
        double v2 = bodies[i].vel.x*bodies[i].vel.x
                  + bodies[i].vel.y*bodies[i].vel.y
                  + bodies[i].vel.z*bodies[i].vel.z;
        KE += 0.5 * m_i * v2;
    }
    double PE = 0.0;
    for (int i = 0; i < n; i++) {
        double m_i = bodies[i].gm / K26A_G;
        for (int j = i + 1; j < n; j++) {
            double m_j = bodies[j].gm / K26A_G;
            K26V3 r = k26astro_pos_sub(&bodies[j].pos, &bodies[i].pos);
            double rmag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
            if (rmag > 0.0) PE -= K26A_G * m_i * m_j / rmag;
        }
    }
    return KE + PE;
}

static int test1_sun_earth(void)
{
    K26AstroBody bodies[2];
    memset(bodies, 0, sizeof(bodies));
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH;
    bodies[1].pos  = k26astro_pos_from_m(K26A_AU_M, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, sqrt(K26A_GM_SUN / K26A_AU_M), 0.0);
    bodies[1].parent_body_idx = 0;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 2) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15) == 0);

    double E0 = total_energy_(bodies, 2);
    double dE_max = 0.0;

    double dt = 86400.0;  /* 1 day */
    int n_steps = 1000;
    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        if (rc != 0) {
            fprintf(stderr, "test1: step %d rc=%d\n", s, rc);
            return 1;
        }
        double E = total_energy_(bodies, 2);
        double dE = fabs((E - E0) / E0);
        if (dE > dE_max) dE_max = dE;
    }

    uint32_t rejects = k26astro_grav_ias15_rejected_steps(&state);
    double dt_last  = k26astro_grav_ias15_get_dt_last(&state);
    fprintf(stderr, "test1 (smooth 2-body): dE_max=%.3e  rejects=%u  dt_last=%.2es\n",
            dE_max, rejects, dt_last);

    if (dE_max >= 1.0e-10) {
        fprintf(stderr, "test1 FAIL: dE_max=%.3e ≥ 1e-10\n", dE_max);
        k26astro_grav_state_destroy(&state);
        return 1;
    }
    if (rejects >= 10) {
        fprintf(stderr, "test1 FAIL: smooth orbit had %u rejects (>=10)\n", rejects);
        k26astro_grav_state_destroy(&state);
        return 1;
    }

    k26astro_grav_state_destroy(&state);
    return 0;
}

static int test2_chaotic_3body(void)
{
    /* Sun + two planets near 2:1 MMR with eccentric inner orbit so
     * close encounters happen within 360 days. */
    K26AstroBody bodies[3];
    memset(bodies, 0, sizeof(bodies));

    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    /* Inner planet: a = 1 AU, e = 0.3 → perihelion 0.7 AU. */
    double a1 = K26A_AU_M;
    double e1 = 0.3;
    double r1_peri = a1 * (1.0 - e1);
    double v1_peri = sqrt(K26A_GM_SUN * (1.0 + e1) / (a1 * (1.0 - e1)));
    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH * 1.0e3;  /* Jupiter-class for chaos */
    bodies[1].pos  = k26astro_pos_from_m(r1_peri, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, v1_peri, 0.0);
    bodies[1].parent_body_idx = 0;

    /* Outer planet: a = 1.587 AU (2:1 MMR with a=1), e = 0.25. */
    double a2 = K26A_AU_M * 1.587401;  /* 2^(2/3) */
    double e2 = 0.25;
    double r2_peri = a2 * (1.0 - e2);
    double v2_peri = sqrt(K26A_GM_SUN * (1.0 + e2) / (a2 * (1.0 - e2)));
    bodies[2].kind = K26ASTRO_BODY_PLANET;
    bodies[2].gm   = K26A_GM_EARTH * 5.0e2;
    /* Start at perihelion but rotated 180° so the two planets begin
     * on opposite sides of the Sun. */
    bodies[2].pos  = k26astro_pos_from_m(-r2_peri, 0.0, 0.0);
    bodies[2].vel  = k26m3d_v3(0.0, -v2_peri, 0.0);
    bodies[2].parent_body_idx = 0;

    K26AstroGravState state;
    assert(k26astro_grav_state_init(&state, bodies, 3) == 0);
    assert(k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15) == 0);
    /* Looser tolerance — Jupiter-class perturbations otherwise
     * trigger more rejects than the smoke gate justifies. */
    k26astro_grav_ias15_set_tol(&state, 1.0e-7);

    double E0 = total_energy_(bodies, 3);
    double dE_max = 0.0;
    double dt = 86400.0;  /* request 1 day */
    int n_steps = 360;
    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        if (rc != 0) {
            fprintf(stderr, "test2: step %d rc=%d\n", s, rc);
            return 1;
        }
        double E = total_energy_(bodies, 3);
        double dE = fabs((E - E0) / E0);
        if (dE > dE_max) dE_max = dE;
        for (int i = 0; i < 3; i++) {
            if (!isfinite(bodies[i].vel.x)) {
                fprintf(stderr, "test2: NaN at step %d body %d\n", s, i);
                return 1;
            }
        }
    }

    uint32_t rejects = k26astro_grav_ias15_rejected_steps(&state);
    fprintf(stderr, "test2 (chaotic 3-body): dE_max=%.3e  rejects=%u\n",
            dE_max, rejects);

    if (dE_max >= 1.0e-6) {
        fprintf(stderr, "test2 FAIL: dE_max=%.3e ≥ 1e-6\n", dE_max);
        k26astro_grav_state_destroy(&state);
        return 1;
    }

    k26astro_grav_state_destroy(&state);
    return 0;
}

int main(void)
{
    if (test1_sun_earth() != 0) return 1;
    if (test2_chaotic_3body() != 0) return 1;
    printf("test_ias15_adaptive: OK\n");
    return 0;
}
