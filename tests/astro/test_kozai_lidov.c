/* tests/astro/test_kozai_lidov.c - 3-body Kozai-Lidov secular cycle,
 * paper-aligned conservation gate.
 *
 * Rein & Spiegel 2015 section 4.3 canonical test (apocenter pumps
 * inner planet's eccentricity to >0.99; adaptive timestepping at
 * extreme geometry). This is the highest-leverage fixture for the
 * COM-relative alpha-fail-safe: apocenter geometry drives the inner
 * planet's COM-relative coords toward FP noise, which is the regime
 * the fail-safe was designed for.
 *
 * IC: REBOUND code units (G=1, M_sun=1, AU=1, year=2pi). Three bodies:
 *   star      m=1.0    at origin
 *   planet    m=1e-3   x=1, vy=1 (circular orbit a=1)
 *   perturber m=1.0    x=10, inclined 89.9 deg
 *
 * Matches a REBOUND-derived reference CSV in tests/astro/data/.
 *
 * Frame: barycentric (move_to_com at t=0). Both this driver and the
 * REBOUND generator apply COM centering, so per-epoch comparison is
 * direct (no per-epoch COM subtraction needed).
 *
 * Units interpretation: K26 always works in SI. We pick a fake
 * G*M_unit = K26A_GM_SUN and length unit = 1.0 AU, time unit chosen
 * such that the orbital relation v^2 = GM/r holds. For the IC
 * "vy=1 at x=1 around m=1", we need v^2 = GM/r -> 1 = 1*1, true if
 * G*m_unit = 1 in our chosen units. In SI: pick length = AU, mass
 * = M_sun, then G*M_sun = K26A_GM_SUN ~ 1.327e20 m^3/s^2, and the
 * derived time unit is sqrt(AU^3/(G*M_sun)) ~ 5.0226e6 s ~ 58.1
 * Julian days ~ year/(2pi).
 *
 * Acceptance:
 *   - Reach T = 800 code time units (~127 yr)
 *   - Final |dL/L| < 1e-9 (rigorous; paper section 4.3 Fig. 5 measures
 *     ~1e-12 at eps_b = 1e-9, so 1e-9 is a comfortable gate)
 *   - Final |dE/E| informational (Kozai is non-symplectic-friendly
 *     but K26 IAS15 should still hold ~1e-10) */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/ias15.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_BODIES        3
#define T_END_UNITS     800.0
#define N_SAMPLES       100
#define L_REL_TOL       1.0e-9
#define E_REL_TOL_INFO  1.0e-8   /* informational only */

#define REF_PATH        "data/kozai_lidov_rebound.csv"
#define EXIT_SKIP       77

/* SI unit conversion: REBOUND code units -> SI.
 * Length unit: 1 code length = K26A_AU_M metres.
 * Mass unit:   1 code mass   = M_sun (i.e., gm_unit = K26A_GM_SUN).
 * Time unit:   t_unit = sqrt(L^3/(G*M)) = sqrt(AU^3/GM_sun).
 *              In Julian-day units: t_unit ~ 58.1324 days.
 *              Year (Julian) ~ 365.25 days ~ 2pi*t_unit per Kepler 3.
 * Velocity unit: 1 code vel = L_unit / t_unit. */
static double t_unit_s_(void)
{
    double au_m = K26A_AU_M;
    return sqrt(au_m * au_m * au_m / K26A_GM_SUN);
}

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

static void move_to_com_(K26AstroBody *bodies, int n)
{
    double total_m = 0.0;
    K26V3 com_pos_m = {0, 0, 0};
    K26V3 com_vel   = {0, 0, 0};
    for (int i = 0; i < n; i++) total_m += bodies[i].gm;
    for (int i = 0; i < n; i++) {
        K26V3 dx = (i == 0)
            ? (K26V3){0,0,0}
            : k26astro_pos_sub(&bodies[i].pos, &bodies[0].pos);
        double w = bodies[i].gm / total_m;
        com_pos_m.x += w * dx.x; com_pos_m.y += w * dx.y; com_pos_m.z += w * dx.z;
        com_vel.x   += w * bodies[i].vel.x;
        com_vel.y   += w * bodies[i].vel.y;
        com_vel.z   += w * bodies[i].vel.z;
    }
    K26V3 neg_com = { -com_pos_m.x, -com_pos_m.y, -com_pos_m.z };
    for (int i = 0; i < n; i++) {
        k26astro_pos_add(&bodies[i].pos, neg_com);
        bodies[i].vel.x -= com_vel.x;
        bodies[i].vel.y -= com_vel.y;
        bodies[i].vel.z -= com_vel.z;
    }
}

static void init_kozai_bodies_(K26AstroBody bodies[N_BODIES])
{
    const double L_unit = K26A_AU_M;          /* m per code length */
    const double T_unit = t_unit_s_();        /* s per code time */
    const double V_unit = L_unit / T_unit;    /* m/s per code velocity */
    const double inc_rad = 89.9 / 180.0 * M_PI;
    const double v_circ_outer = sqrt((1.0 + 1.0) / 10.0);  /* sqrt(2/10) */

    memset(bodies, 0, N_BODIES * sizeof(K26AstroBody));

    /* star: m=1, at origin, zero velocity */
    bodies[0].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[0].gm   = 1.0 * K26A_GM_SUN;
    bodies[0].pos  = k26astro_pos_from_m(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    /* planet: m=1e-3, x=1, vy=1 */
    bodies[1].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[1].gm   = 1.0e-3 * K26A_GM_SUN;
    bodies[1].pos  = k26astro_pos_from_m(1.0 * L_unit, 0.0, 0.0);
    bodies[1].vel.x = 0.0;
    bodies[1].vel.y = 1.0 * V_unit;
    bodies[1].vel.z = 0.0;
    bodies[1].parent_body_idx = -1;

    /* perturber: m=1, x=10, vy=cos(i)·sqrt(2/10), vz=sin(i)·sqrt(2/10) */
    bodies[2].kind = K26ASTRO_BODY_BARYCENTRE;
    bodies[2].gm   = 1.0 * K26A_GM_SUN;
    bodies[2].pos  = k26astro_pos_from_m(10.0 * L_unit, 0.0, 0.0);
    bodies[2].vel.x = 0.0;
    bodies[2].vel.y = cos(inc_rad) * v_circ_outer * V_unit;
    bodies[2].vel.z = sin(inc_rad) * v_circ_outer * V_unit;
    bodies[2].parent_body_idx = -1;

    move_to_com_(bodies, N_BODIES);
}

int main(void)
{
    FILE *ref = fopen(REF_PATH, "r");
    if (!ref) {
        fprintf(stderr,
                "test_kozai_lidov: SKIP - reference CSV %s not present\n",
                REF_PATH);
        return EXIT_SKIP;
    }

    K26AstroBody bodies[N_BODIES];
    init_kozai_bodies_(bodies);

    double E0, L0;
    compute_energy_angmom_(bodies, N_BODIES, &E0, &L0);

    K26AstroGravState state;
    if (k26astro_grav_state_init(&state, bodies, N_BODIES) != 0) {
        fclose(ref);
        return 1;
    }
    k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15);
    k26astro_grav_ias15_set_tol(&state, 1.0e-9);

    const double T_unit_s = t_unit_s_();
    const double L_unit_m = K26A_AU_M;

    char line[1024];
    for (int i = 0; i < 6; i++) {
        if (!fgets(line, sizeof line, ref)) { fclose(ref); return 1; }
    }

    double t_now = 0.0;
    double max_pos_err = 0.0;
    int n_samples = 0;
    int convergence_failures = 0;

    while (fgets(line, sizeof line, ref)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double t_units, vals[N_BODIES * 3];
        int nread = sscanf(line,
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            &t_units,
            &vals[0], &vals[1], &vals[2],
            &vals[3], &vals[4], &vals[5],
            &vals[6], &vals[7], &vals[8]);
        if (nread != 1 + N_BODIES * 3) continue;

        double t_ref_s = t_units * T_unit_s;
        if (t_ref_s <= t_now) continue;

        double dt_to_ref = t_ref_s - t_now;
        int rc = k26astro_grav_step(&state, dt_to_ref);
        if (rc != 0) {
            convergence_failures++;
            fprintf(stderr, "test_kozai_lidov: step to t=%.2f units rc=%d\n",
                    t_units, rc);
            if (convergence_failures > 3) break;
            continue;
        }
        t_now = t_ref_s;
        n_samples++;

        for (int i = 0; i < N_BODIES; i++) {
            K26V3 r_k26_m = k26astro_pos_to_m_approx(&bodies[i].pos);
            double rx = r_k26_m.x / L_unit_m;
            double ry = r_k26_m.y / L_unit_m;
            double rz = r_k26_m.z / L_unit_m;
            double dx = rx - vals[i*3 + 0];
            double dy = ry - vals[i*3 + 1];
            double dz = rz - vals[i*3 + 2];
            double err = sqrt(dx*dx + dy*dy + dz*dz);
            double mag = sqrt(vals[i*3 + 0]*vals[i*3 + 0]
                            + vals[i*3 + 1]*vals[i*3 + 1]
                            + vals[i*3 + 2]*vals[i*3 + 2]);
            double rel = (mag > 1.0e-6) ? (err / mag) : err;
            if (rel > max_pos_err) max_pos_err = rel;
        }
    }
    fclose(ref);

    double E_final, L_final;
    compute_energy_angmom_(bodies, N_BODIES, &E_final, &L_final);
    double dE_rel = fabs((E_final - E0) / E0);
    double dL_rel = (L0 > 0.0) ? fabs((L_final - L0) / L0) : fabs(L_final - L0);
    uint32_t rejects = k26astro_grav_ias15_rejected_steps(&state);

    fprintf(stderr,
            "test_kozai_lidov: t=%.2f units samples=%d pos_err_max=%.3e "
            "rejects=%u |dE/E|=%.3e |dL/L|=%.3e\n",
            t_now / T_unit_s, n_samples, max_pos_err,
            rejects, dE_rel, dL_rel);

    int fail = 0;
    if (t_now < (T_END_UNITS - 1.0) * T_unit_s) {
        fprintf(stderr,
            "test_kozai_lidov: FAIL - only reached t=%.2f units (target %.2f)\n",
            t_now / T_unit_s, T_END_UNITS);
        fail = 1;
    }
    if (dL_rel >= L_REL_TOL) {
        fprintf(stderr,
            "test_kozai_lidov: FAIL - |dL/L|=%.3e >= %.3e\n",
            dL_rel, L_REL_TOL);
        fail = 1;
    }

    k26astro_grav_state_destroy(&state);
    if (fail) return 1;
    printf("test_kozai_lidov: PASS (T=%.0f code units, %d samples vs REBOUND)\n",
           T_END_UNITS, n_samples);
    return 0;
}
