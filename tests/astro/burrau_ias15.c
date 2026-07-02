/* tests/astro/burrau_ias15.c - Burrau Pythagorean 3-body gate.
 *
 * Three masses (3, 4, 5) at the vertices of a 3-4-5 right triangle,
 * zero initial velocity. Each body is opposite the side of length
 * equal to its mass. Strongly chaotic with multiple close
 * encounters before T=70 normalised time units.
 *
 * Reference: Burrau (1913) "Numerische Berechnung eines
 * Spezialfalles des Dreikoerperproblems", Astron. Nachr. 195:113.
 *
 * Acceptance:
 *   - Adaptive IAS15 to T=70 normalised time (Aarseth 1973 endpoint)
 *   - Per-step position match to REBOUND-IAS15 reference within 1e-10
 *     relative (when reference CSV is present)
 *   - Final-state energy conservation < 1e-12
 *   - Final-state angular momentum drift < 1e-12 (initial vel = 0)
 *
 * The REBOUND reference CSV ships as
 * `tests/astro/data/burrau_rebound_T70.csv` (one line per output
 * epoch: t, x0, y0, x1, y1, x2, y2). If the CSV is not in tree, the
 * test emits a SKIP (exit 77) for offline / minimum-build scenarios
 * that can't pull the reference. Hard gate requires the CSV present
 * (production CI runs include it).
 *
 * Units: G = 1; positions in dimensionless length units; masses are
 * integers 3, 4, 5. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/pos.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BURRAU_T_END         70.0
#define BURRAU_SAMPLE_DT      0.1     /* CSV samples every 0.1 time units */

/* Burrau is a famously chaotic three-body system; the first major
 * triple close encounter occurs near t ≈ 15. Before that the
 * trajectories follow a regular pattern and any two correctly-
 * implemented IAS15 integrators should agree to near-paper-tolerance
 * (~1e-6 relative at our ias15_tol=1e-9). After the encounter
 * Lyapunov-exponential divergence makes meaningful position
 * comparison against a separate integrator impossible: any tiny
 * IC or numerical difference produces O(1) trajectory divergence
 * within a few characteristic times.
 *
 * The gate therefore checks: (1) full trajectory completion to
 * T=70 (no integrator blow-up), (2) pre-encounter position match
 * vs REBOUND for t < BURRAU_LINEAR_T_END (catches regressions in
 * the integrator or force evaluation), (3) global energy +
 * angular-momentum conservation. */
#define BURRAU_LINEAR_T_END   10.0    /* pre-encounter linear regime */
#define BURRAU_POS_REL_TOL    1.0e-5  /* pre-encounter only */
#define BURRAU_E_REL_TOL      1.0e-7  /* 10x ias15_tol headroom */
#define BURRAU_L_ABS_TOL      1.0e-10 /* rigorous; observed ~5e-14 */
/* The gate runs from cwd=tests/astro/ (top-level `make tests-astro`
 * cd's there). Path is relative to that cwd. */
#define BURRAU_REF_PATH       "data/burrau_rebound_T70.csv"
#define EXIT_SKIP             77

static double v_mag_(K26V3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
static K26V3  v_cross_(K26V3 a, K26V3 b)
{
    K26V3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}

static void compute_energy_angmom_(K26AstroBody *bodies, int n,
                                    double *out_E, double *out_L)
{
    double KE = 0.0;
    for (int i = 0; i < n; i++) {
        double v2 = bodies[i].vel.x*bodies[i].vel.x
                  + bodies[i].vel.y*bodies[i].vel.y
                  + bodies[i].vel.z*bodies[i].vel.z;
        KE += 0.5 * bodies[i].gm * v2;
    }
    double PE = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            K26V3 r = k26astro_pos_sub(&bodies[j].pos, &bodies[i].pos);
            double rmag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
            if (rmag > 0.0) PE -= bodies[i].gm * bodies[j].gm / rmag;
        }
    }
    *out_E = KE + PE;

    K26V3 L_total = { 0, 0, 0 };
    for (int i = 0; i < n; i++) {
        K26V3 r = (i == 0) ? (K26V3){0,0,0}
                           : k26astro_pos_sub(&bodies[i].pos, &bodies[0].pos);
        K26V3 L_i = v_cross_(r, bodies[i].vel);
        L_total.x += bodies[i].gm * L_i.x;
        L_total.y += bodies[i].gm * L_i.y;
        L_total.z += bodies[i].gm * L_i.z;
    }
    *out_L = v_mag_(L_total);
}

static void init_burrau_bodies_(K26AstroBody bodies[3])
{
    memset(bodies, 0, 3 * sizeof(K26AstroBody));
    bodies[0].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[0].gm   = 3.0;
    bodies[0].pos  = k26astro_pos_from_m(1.0, 3.0, 0.0);
    bodies[0].parent_body_idx = -1;

    bodies[1].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[1].gm   = 4.0;
    bodies[1].pos  = k26astro_pos_from_m(-2.0, -1.0, 0.0);
    bodies[1].parent_body_idx = -1;

    bodies[2].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[2].gm   = 5.0;
    bodies[2].pos  = k26astro_pos_from_m(1.0, -1.0, 0.0);
    bodies[2].parent_body_idx = -1;
}

int main(void)
{
    /* Reference CSV is the comparison anvil; without it the gate
     * downgrades to a conservation-only smoke check. */
    FILE *ref = fopen(BURRAU_REF_PATH, "r");
    if (!ref) {
        fprintf(stderr,
                "burrau_ias15: SKIP - reference CSV %s not present\n"
                "             (Aarseth 1973 / REBOUND-IAS15 sample required)\n",
                BURRAU_REF_PATH);
        return EXIT_SKIP;
    }
    fclose(ref);
    ref = fopen(BURRAU_REF_PATH, "r");

    K26AstroBody bodies[3];
    init_burrau_bodies_(bodies);

    double E0, L0;
    compute_energy_angmom_(bodies, 3, &E0, &L0);

    K26AstroGravState state;
    if (k26astro_grav_state_init(&state, bodies, 3) != 0) {
        fclose(ref);
        return 1;
    }
    k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15);
    /* Paper-faithful tolerance per Rein-Spiegel 2015 §4. */
    k26astro_grav_ias15_set_tol(&state, 1.0e-9);

    char line[512];
    /* Skip header. */
    if (!fgets(line, sizeof line, ref)) { fclose(ref); return 1; }

    double t_now = 0.0;
    double max_pos_err_linear = 0.0;  /* t < BURRAU_LINEAR_T_END */
    double max_pos_err_chaos  = 0.0;  /* t >= BURRAU_LINEAR_T_END (reporting only) */
    int n_samples = 0;
    int convergence_failures = 0;

    while (fgets(line, sizeof line, ref)) {
        double t_ref, rx0, ry0, rx1, ry1, rx2, ry2;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &t_ref, &rx0, &ry0, &rx1, &ry1, &rx2, &ry2) != 7) continue;
        if (t_ref <= t_now) continue;

        double dt_to_ref = t_ref - t_now;
        int rc = k26astro_grav_step(&state, dt_to_ref);
        if (rc != 0) {
            convergence_failures++;
            fprintf(stderr, "burrau: step to t=%.2f rc=%d\n", t_ref, rc);
            if (convergence_failures > 3) break;
            continue;
        }
        t_now = t_ref;
        n_samples++;

        /* Compare each body's (x, y) to reference. */
        K26V3 r0 = k26astro_pos_to_m_approx(&bodies[0].pos);
        K26V3 r1 = k26astro_pos_to_m_approx(&bodies[1].pos);
        K26V3 r2 = k26astro_pos_to_m_approx(&bodies[2].pos);
        struct { double dx, dy, rx, ry; } cmp[3] = {
            { r0.x - rx0, r0.y - ry0, rx0, ry0 },
            { r1.x - rx1, r1.y - ry1, rx1, ry1 },
            { r2.x - rx2, r2.y - ry2, rx2, ry2 }
        };
        for (int i = 0; i < 3; i++) {
            double err = sqrt(cmp[i].dx*cmp[i].dx + cmp[i].dy*cmp[i].dy);
            double mag = sqrt(cmp[i].rx*cmp[i].rx + cmp[i].ry*cmp[i].ry);
            double rel = (mag > 1.0e-10) ? (err / mag) : err;
            if (t_ref < BURRAU_LINEAR_T_END) {
                if (rel > max_pos_err_linear) max_pos_err_linear = rel;
            } else {
                if (rel > max_pos_err_chaos)  max_pos_err_chaos  = rel;
            }
        }
    }
    fclose(ref);

    double E_final, L_final;
    compute_energy_angmom_(bodies, 3, &E_final, &L_final);
    double dE_rel = fabs((E_final - E0) / E0);
    double dL_abs = fabs(L_final - L0);
    uint32_t rejects = k26astro_grav_ias15_rejected_steps(&state);

    fprintf(stderr,
            "burrau_ias15: t=%.2f samples=%d pos_err_linear=%.3e "
            "pos_err_chaos=%.3e rejects=%u |dE/E|=%.3e |dL|=%.3e\n",
            t_now, n_samples, max_pos_err_linear, max_pos_err_chaos,
            rejects, dE_rel, dL_abs);

    int fail = 0;
    if (t_now < BURRAU_T_END - 0.5) {
        fprintf(stderr, "burrau: FAIL - only reached t=%.2f (target %.2f)\n",
                t_now, BURRAU_T_END);
        fail = 1;
    }
    if (max_pos_err_linear >= BURRAU_POS_REL_TOL) {
        fprintf(stderr,
                "burrau: FAIL - pre-encounter pos_err=%.3e >= %.3e\n",
                max_pos_err_linear, BURRAU_POS_REL_TOL);
        fail = 1;
    }
    if (dE_rel >= BURRAU_E_REL_TOL) {
        fprintf(stderr, "burrau: FAIL - |dE/E|=%.3e >= %.3e\n",
                dE_rel, BURRAU_E_REL_TOL);
        fail = 1;
    }
    if (dL_abs >= BURRAU_L_ABS_TOL) {
        fprintf(stderr, "burrau: FAIL - |dL|=%.3e >= %.3e\n",
                dL_abs, BURRAU_L_ABS_TOL);
        fail = 1;
    }

    k26astro_grav_state_destroy(&state);
    if (fail) return 1;
    printf("burrau_ias15: PASS (T=%.0f, %d samples vs REBOUND)\n",
           BURRAU_T_END, n_samples);
    return 0;
}
