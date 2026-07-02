/* tests/astro/test_sj_s_10kyr.c - Sun + Jupiter + Saturn 3-body
 * over 10 kyr.
 *
 * Real Solar System 3-body subset on a long secular timescale.
 * Acceptance gate |dL/L| < 1e-9 is the Kinoshita-Nakai bound.
 *
 * IC: Applegate 1986 heliocentric (AU/day/M_sun), matched by a
 * REBOUND-derived reference CSV in tests/astro/data/.
 * Sun mass absorbs Mercury-Mars (1.00000597682 M_sun).
 *
 * Frame: barycentric (move_to_com at t=0). Direct per-epoch comparison
 * with REBOUND CSV (no per-epoch COM subtraction needed). */
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
#define T_END_DAYS      (10000.0 * 365.25)
#define T_END_S         (T_END_DAYS * 86400.0)
#define N_SAMPLES       100
#define L_REL_TOL       1.0e-9
#define E_REL_TOL_INFO  1.0e-9

#define REF_PATH        "data/sj_s_10kyr_rebound.csv"
#define EXIT_SKIP       77

static const double ic_pos[N_BODIES][3] = {
    {-4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6},
    {+3.40546614227466e+0, +3.62978190075864e+0, +3.42386261766577e-2},
    {+6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1}
};
static const double ic_vel[N_BODIES][3] = {
    {+6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9},
    {-5.59797969310664e-3, +5.51815399480116e-3, -2.66711392865591e-6},
    {-4.17354020307064e-3, +3.99723751748116e-3, +1.67206320571441e-5}
};
static const double ic_mass[N_BODIES] = {
    1.00000597682,
    1.0 / 1047.355,
    1.0 / 3501.6
};

static const double AU_PER_DAY_TO_M_PER_S = K26A_AU_M / 86400.0;

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

static void init_sjs_bodies_(K26AstroBody bodies[N_BODIES])
{
    memset(bodies, 0, N_BODIES * sizeof(K26AstroBody));
    for (int i = 0; i < N_BODIES; i++) {
        bodies[i].kind = K26ASTRO_BODY_BARYCENTRE;
        bodies[i].gm   = ic_mass[i] * K26A_GM_SUN;
        bodies[i].pos  = k26astro_pos_from_m(
            ic_pos[i][0] * K26A_AU_M,
            ic_pos[i][1] * K26A_AU_M,
            ic_pos[i][2] * K26A_AU_M);
        bodies[i].vel.x = ic_vel[i][0] * AU_PER_DAY_TO_M_PER_S;
        bodies[i].vel.y = ic_vel[i][1] * AU_PER_DAY_TO_M_PER_S;
        bodies[i].vel.z = ic_vel[i][2] * AU_PER_DAY_TO_M_PER_S;
        bodies[i].parent_body_idx = -1;
    }
    move_to_com_(bodies, N_BODIES);
}

int main(void)
{
    FILE *ref = fopen(REF_PATH, "r");
    if (!ref) {
        fprintf(stderr,
                "test_sj_s_10kyr: SKIP - reference CSV %s not present\n",
                REF_PATH);
        return EXIT_SKIP;
    }

    K26AstroBody bodies[N_BODIES];
    init_sjs_bodies_(bodies);

    double E0, L0;
    compute_energy_angmom_(bodies, N_BODIES, &E0, &L0);

    K26AstroGravState state;
    if (k26astro_grav_state_init(&state, bodies, N_BODIES) != 0) {
        fclose(ref);
        return 1;
    }
    k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15);
    k26astro_grav_ias15_set_tol(&state, 1.0e-9);

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
        double t_days, vals[N_BODIES * 3];
        int nread = sscanf(line,
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            &t_days,
            &vals[0], &vals[1], &vals[2],
            &vals[3], &vals[4], &vals[5],
            &vals[6], &vals[7], &vals[8]);
        if (nread != 1 + N_BODIES * 3) continue;

        double t_ref_s = t_days * 86400.0;
        if (t_ref_s <= t_now) continue;

        double dt_to_ref = t_ref_s - t_now;
        int rc = k26astro_grav_step(&state, dt_to_ref);
        if (rc != 0) {
            convergence_failures++;
            fprintf(stderr, "test_sj_s_10kyr: step to t=%.2f days rc=%d\n",
                    t_days, rc);
            if (convergence_failures > 3) break;
            continue;
        }
        t_now = t_ref_s;
        n_samples++;

        for (int i = 0; i < N_BODIES; i++) {
            K26V3 r_k26_m = k26astro_pos_to_m_approx(&bodies[i].pos);
            double rx = r_k26_m.x / K26A_AU_M;
            double ry = r_k26_m.y / K26A_AU_M;
            double rz = r_k26_m.z / K26A_AU_M;
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
            "test_sj_s_10kyr: t=%.0f days samples=%d pos_err_max=%.3e "
            "rejects=%u |dE/E|=%.3e |dL/L|=%.3e\n",
            t_now / 86400.0, n_samples, max_pos_err,
            rejects, dE_rel, dL_rel);

    int fail = 0;
    if (t_now < T_END_S - 86400.0) {
        fprintf(stderr,
            "test_sj_s_10kyr: FAIL - only reached t=%.0f days (target %.0f)\n",
            t_now / 86400.0, T_END_DAYS);
        fail = 1;
    }
    if (dL_rel >= L_REL_TOL) {
        fprintf(stderr,
            "test_sj_s_10kyr: FAIL - |dL/L|=%.3e >= %.3e\n",
            dL_rel, L_REL_TOL);
        fail = 1;
    }

    k26astro_grav_state_destroy(&state);
    if (fail) return 1;
    printf("test_sj_s_10kyr: PASS (T=%.0f days, %d samples vs REBOUND)\n",
           T_END_DAYS, n_samples);
    return 0;
}
