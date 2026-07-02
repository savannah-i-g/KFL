/* test_ias15_rollback.c — IAS15 explicit snapshot/rollback gate.
 *
 * Two coverage axes for the byte-snapshot rollback machinery:
 *
 *   Test 1: Determinism under rejects. Run a chaotic Jupiter-class
 *     3-body scenario across 360 days at a tight tolerance that
 *     forces rejects during close approaches. Run the same scenario
 *     twice in fully independent K26AstroGravState instances.
 *     Reject counters must match; final body[] arrays must be
 *     bit-identical (memcmp == 0). Proves the rollback mechanism
 *     doesn't introduce reject-path non-determinism.
 *
 *   Test 2: Reject paths actually fire. The same scenario must
 *     produce strictly positive `ias15_rejected_steps` — otherwise
 *     test 1 trivially passes without exercising rollback. Tight
 *     tol (1e-12) on a chaotic 3-body forces the controller to
 *     reject around close approaches.
 *
 * Acceptance:
 *   - rejects > 0 (confirms snapshot/restore path actually fired)
 *   - memcmp(A, B) == 0 over the bodies arrays after 360 steps
 *   - reject counts match between runs A and B
 *   - dE_max stays below 1e-6 (sanity check that the integrator
 *     is still doing useful work, not running into NaN territory) */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void seed_chaotic_3body_(K26AstroBody *bodies)
{
    memset(bodies, 0, 3 * sizeof(K26AstroBody));

    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    /* Inner planet: a=1 AU, e=0.4 → perihelion 0.6 AU. */
    double a1 = K26A_AU_M;
    double e1 = 0.4;
    double r1_peri = a1 * (1.0 - e1);
    double v1_peri = sqrt(K26A_GM_SUN * (1.0 + e1) / (a1 * (1.0 - e1)));
    bodies[1].kind = K26ASTRO_BODY_PLANET;
    bodies[1].gm   = K26A_GM_EARTH * 1.0e3;
    bodies[1].pos  = k26astro_pos_from_m(r1_peri, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, v1_peri, 0.0);
    bodies[1].parent_body_idx = 0;

    double a2 = K26A_AU_M * 1.587401;
    double e2 = 0.3;
    double r2_peri = a2 * (1.0 - e2);
    double v2_peri = sqrt(K26A_GM_SUN * (1.0 + e2) / (a2 * (1.0 - e2)));
    bodies[2].kind = K26ASTRO_BODY_PLANET;
    bodies[2].gm   = K26A_GM_EARTH * 5.0e2;
    bodies[2].pos  = k26astro_pos_from_m(-r2_peri, 0.0, 0.0);
    bodies[2].vel  = k26m3d_v3(0.0, -v2_peri, 0.0);
    bodies[2].parent_body_idx = 0;
}

static double total_energy_(const K26AstroBody *bodies, int n)
{
    double KE = 0.0, PE = 0.0;
    for (int i = 0; i < n; i++) {
        double m_i = bodies[i].gm / K26A_G;
        double v2  = bodies[i].vel.x*bodies[i].vel.x
                   + bodies[i].vel.y*bodies[i].vel.y
                   + bodies[i].vel.z*bodies[i].vel.z;
        KE += 0.5 * m_i * v2;
        for (int j = i + 1; j < n; j++) {
            double m_j = bodies[j].gm / K26A_G;
            K26V3 r = k26astro_pos_sub(&bodies[j].pos, &bodies[i].pos);
            double rmag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
            if (rmag > 0.0) PE -= K26A_G * m_i * m_j / rmag;
        }
    }
    return KE + PE;
}

static int run_chaotic_(K26AstroBody *bodies,
                        uint32_t *out_rejects,
                        double *out_dE)
{
    K26AstroGravState state;
    if (k26astro_grav_state_init(&state, bodies, 3) != 0) return 1;
    if (k26astro_grav_set_integrator(&state,
                                     K26ASTRO_INTEGRATOR_IAS15) != 0) return 1;
    /* Default tol (1e-9). Rejects are forced not by tightening tol
     * but by requesting an enormous step size that the controller
     * must shrink — every shrink-reject iteration exercises the
     * snapshot save/restore path. dt=100 days requested in one
     * call → for the e=0.4 inner planet (period ~365 days) the
     * controller will reject and re-subdivide many times. */
    double E0 = total_energy_(bodies, 3);
    double dE_max = 0.0;

    double dt = 86400.0 * 100.0;  /* 100 days requested in one shot */
    int n_steps = 4;             /* 4 × 100 = 400 days total */
    for (int s = 0; s < n_steps; s++) {
        int rc = k26astro_grav_step(&state, dt);
        if (rc != 0) {
            fprintf(stderr, "run: step %d rc=%d\n", s, rc);
            k26astro_grav_state_destroy(&state);
            return 1;
        }
        double E = total_energy_(bodies, 3);
        double dE = fabs((E - E0) / E0);
        if (dE > dE_max) dE_max = dE;
        for (int i = 0; i < 3; i++) {
            if (!isfinite(bodies[i].vel.x)) {
                fprintf(stderr, "run: NaN at step %d body %d\n", s, i);
                k26astro_grav_state_destroy(&state);
                return 1;
            }
        }
    }

    *out_rejects = k26astro_grav_ias15_rejected_steps(&state);
    *out_dE      = dE_max;
    k26astro_grav_state_destroy(&state);
    return 0;
}

static int test_deterministic_rejects(void)
{
    K26AstroBody A[3], B[3];
    seed_chaotic_3body_(A);
    seed_chaotic_3body_(B);
    assert(memcmp(A, B, sizeof A) == 0);

    uint32_t rj_A = 0, rj_B = 0;
    double dE_A = 0.0, dE_B = 0.0;
    if (run_chaotic_(A, &rj_A, &dE_A)) return 1;
    if (run_chaotic_(B, &rj_B, &dE_B)) return 1;

    fprintf(stderr,
        "test (determinism under rejects): rejects A=%u B=%u  "
        "dE_A=%.3e dE_B=%.3e\n", rj_A, rj_B, dE_A, dE_B);

    if (rj_A == 0) {
        fprintf(stderr,
            "FAIL: zero rejects — scenario didn't exercise rollback "
            "(tighten tol or extend horizon)\n");
        return 1;
    }
    if (rj_A != rj_B) {
        fprintf(stderr,
            "FAIL: reject counts differ A=%u B=%u — non-deterministic "
            "rollback path\n", rj_A, rj_B);
        return 1;
    }
    if (memcmp(A, B, sizeof A) != 0) {
        fprintf(stderr,
            "FAIL: bodies[] not bit-identical between runs\n");
        for (int i = 0; i < 3; i++) {
            const unsigned char *a = (const unsigned char *)&A[i];
            const unsigned char *b = (const unsigned char *)&B[i];
            for (size_t k = 0; k < sizeof(K26AstroBody); k++) {
                if (a[k] != b[k]) {
                    fprintf(stderr,
                      "  body[%d] byte %zu: A=0x%02x B=0x%02x\n",
                      i, k, a[k], b[k]);
                    break;
                }
            }
        }
        return 1;
    }
    if (dE_A > 1.0e-6) {
        fprintf(stderr, "FAIL: dE_A=%.3e exceeds 1e-6\n", dE_A);
        return 1;
    }

    return 0;
}

int main(void)
{
    if (test_deterministic_rejects()) return 1;
    fprintf(stderr, "test_ias15_rollback: OK\n");
    return 0;
}
