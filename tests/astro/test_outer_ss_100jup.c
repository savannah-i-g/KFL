/* tests/astro/test_outer_ss_100jup.c - outer Solar System
 * (Sun + Jup + Sat + Ura + Nep + Pluto) over 100 Jupiter orbits,
 * paper-aligned conservation gate.
 *
 * Rein & Spiegel 2015 section 4.1 canonical test. Paper Fig. 1
 * shows IAS15 reaching |dE/E| ~ 1e-13 across this run at
 * eps_b = 1e-9.
 *
 * IC: Applegate 1986 heliocentric (AU/day/M_sun), matched by a
 * REBOUND-derived reference CSV in tests/astro/data/.
 * Sun mass absorbs Mercury-Mars (1.00000597682 M_sun). Pluto is a
 * test particle (m=0).
 *
 * Frame: barycentric. Both this driver and the REBOUND generator
 * apply move-to-COM at t=0 -> COM at origin, comparison is direct
 * (no per-epoch COM subtraction needed; K26's COM stays at origin
 * modulo numerical drift).
 *
 * Acceptance:
 *   - Reach T = 100 * 4332.59 days = 433259 days exactly
 *   - Final |dE/E| < 1e-12 (paper section 4.1 bound)
 *   - Final |dL/L| < 1e-12 (informational)
 *   - Pre-chaos position match vs REBOUND CSV: rel < 1e-7 AU for first
 *     50 Jupiter orbits (informational; outer SS is non-chaotic on
 *     these timescales but small phase drift accumulates) */
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

#define T_JUPITER_DAYS      4332.59
#define T_END_DAYS          (100.0 * T_JUPITER_DAYS)
#define SAMPLE_DT_DAYS      (10.0 * 365.25)
#define T_END_S             (T_END_DAYS * 86400.0)
#define SAMPLE_DT_S         (SAMPLE_DT_DAYS * 86400.0)

#define N_BODIES            6
#define E_REL_TOL           1.0e-12
#define POS_REL_TOL_LINEAR  1.0e-7   /* per-epoch, first 50 Jup orbits */

#define REF_PATH            "data/outer_ss_100jup_rebound.csv"
#define EXIT_SKIP           77

/* Applegate 1986 IC (AU, AU/day, M_sun). Identical to the REBOUND
 * generator and the source it was derived from. */
static const double ic_pos[N_BODIES][3] = {
    {-4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6},
    {+3.40546614227466e+0, +3.62978190075864e+0, +3.42386261766577e-2},
    {+6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1},
    {+1.11636331405597e+1, +1.60373479057256e+1, +3.61783279369958e-1},
    {-3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1},
    {-2.13858977531573e+1, +3.20719104739886e+1, +2.49245689556096e+0}
};
static const double ic_vel[N_BODIES][3] = {
    {+6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9},
    {-5.59797969310664e-3, +5.51815399480116e-3, -2.66711392865591e-6},
    {-4.17354020307064e-3, +3.99723751748116e-3, +1.67206320571441e-5},
    {-3.25884806151064e-3, +2.06438412905916e-3, -2.17699042180559e-5},
    {-2.17471785045538e-4, -3.11361111025884e-3, +3.58344705491441e-5},
    {-1.76936577252484e-3, -2.06720938381724e-3, +6.58091931493844e-4}
};
static const double ic_mass[N_BODIES] = {
    1.00000597682,
    1.0 / 1047.355,
    1.0 / 3501.6,
    1.0 / 22869.0,
    1.0 / 19314.0,
    0.0
};

/* AU/day -> m/s */
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

/* Apply barycentric centering: subtract COM pos+vel from each body
 * so the system COM sits at origin with zero net momentum at t=0.
 * Matches REBOUND's `reb_simulation_move_to_com` semantics. */
static void move_to_com_(K26AstroBody *bodies, int n)
{
    double total_m = 0.0;
    K26V3 com_pos_m = {0, 0, 0};  /* COM offset relative to body[0] in metres */
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
    /* Shift positions: body[i].pos -= (body[0].pos + com_pos_m).
     * Easiest: subtract com_pos_m from each body's offset, then
     * shift body[0] explicitly. */
    K26V3 neg_com = { -com_pos_m.x, -com_pos_m.y, -com_pos_m.z };
    for (int i = 0; i < n; i++) {
        k26astro_pos_add(&bodies[i].pos, neg_com);
        bodies[i].vel.x -= com_vel.x;
        bodies[i].vel.y -= com_vel.y;
        bodies[i].vel.z -= com_vel.z;
    }
}

static void init_outer_ss_bodies_(K26AstroBody bodies[N_BODIES])
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
                "test_outer_ss_100jup: SKIP - reference CSV %s not present\n",
                REF_PATH);
        return EXIT_SKIP;
    }

    K26AstroBody bodies[N_BODIES];
    init_outer_ss_bodies_(bodies);

    double E0, L0;
    compute_energy_angmom_(bodies, N_BODIES, &E0, &L0);

    K26AstroGravState state;
    if (k26astro_grav_state_init(&state, bodies, N_BODIES) != 0) {
        fclose(ref);
        return 1;
    }
    k26astro_grav_set_integrator(&state, K26ASTRO_INTEGRATOR_IAS15);
    k26astro_grav_ias15_set_tol(&state, 1.0e-9);

    /* Skip 5 header lines + 1 column header. */
    char line[1024];
    for (int i = 0; i < 5; i++) {
        if (!fgets(line, sizeof line, ref)) { fclose(ref); return 1; }
    }

    double t_now = 0.0;       /* seconds */
    double max_pos_err_linear = 0.0;
    double max_pos_err_late   = 0.0;
    int n_samples = 0;
    int convergence_failures = 0;
    const double linear_t_end_s = 50.0 * T_JUPITER_DAYS * 86400.0;

    while (fgets(line, sizeof line, ref)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double t_days, vals[N_BODIES * 3];
        int nread = sscanf(line,
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,"
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            &t_days,
            &vals[0], &vals[1], &vals[2],
            &vals[3], &vals[4], &vals[5],
            &vals[6], &vals[7], &vals[8],
            &vals[9], &vals[10], &vals[11],
            &vals[12], &vals[13], &vals[14],
            &vals[15], &vals[16], &vals[17]);
        if (nread != 1 + N_BODIES * 3) continue;

        double t_ref_s = t_days * 86400.0;
        if (t_ref_s <= t_now) continue;

        double dt_to_ref = t_ref_s - t_now;
        int rc = k26astro_grav_step(&state, dt_to_ref);
        if (rc != 0) {
            convergence_failures++;
            fprintf(stderr, "test_outer_ss_100jup: step to t=%.2f days rc=%d\n",
                    t_days, rc);
            if (convergence_failures > 3) break;
            continue;
        }
        t_now = t_ref_s;
        n_samples++;

        for (int i = 0; i < N_BODIES; i++) {
            K26V3 r_k26_m = k26astro_pos_to_m_approx(&bodies[i].pos);
            double rx_AU = r_k26_m.x / K26A_AU_M;
            double ry_AU = r_k26_m.y / K26A_AU_M;
            double rz_AU = r_k26_m.z / K26A_AU_M;
            double dx = rx_AU - vals[i*3 + 0];
            double dy = ry_AU - vals[i*3 + 1];
            double dz = rz_AU - vals[i*3 + 2];
            double err = sqrt(dx*dx + dy*dy + dz*dz);
            double mag = sqrt(vals[i*3 + 0]*vals[i*3 + 0]
                            + vals[i*3 + 1]*vals[i*3 + 1]
                            + vals[i*3 + 2]*vals[i*3 + 2]);
            double rel = (mag > 1.0e-6) ? (err / mag) : err;
            if (t_now < linear_t_end_s) {
                if (rel > max_pos_err_linear) max_pos_err_linear = rel;
            } else {
                if (rel > max_pos_err_late)   max_pos_err_late   = rel;
            }
        }
    }
    fclose(ref);

    double E_final, L_final;
    compute_energy_angmom_(bodies, N_BODIES, &E_final, &L_final);
    double dE_rel = fabs((E_final - E0) / E0);
    double dL_rel = (L0 > 0.0) ? fabs((L_final - L0) / L0) : fabs(L_final - L0);
    uint32_t rejects = k26astro_grav_ias15_rejected_steps(&state);

    fprintf(stderr,
            "test_outer_ss_100jup: t=%.0f days samples=%d "
            "pos_err_linear=%.3e pos_err_late=%.3e rejects=%u "
            "|dE/E|=%.3e |dL/L|=%.3e\n",
            t_now / 86400.0, n_samples,
            max_pos_err_linear, max_pos_err_late,
            rejects, dE_rel, dL_rel);

    int fail = 0;
    if (t_now < T_END_S - 86400.0) {
        fprintf(stderr,
            "test_outer_ss_100jup: FAIL - only reached t=%.0f days (target %.0f)\n",
            t_now / 86400.0, T_END_DAYS);
        fail = 1;
    }
    if (dE_rel >= E_REL_TOL) {
        fprintf(stderr,
            "test_outer_ss_100jup: FAIL - |dE/E|=%.3e >= %.3e\n",
            dE_rel, E_REL_TOL);
        fail = 1;
    }

    k26astro_grav_state_destroy(&state);
    if (fail) return 1;
    printf("test_outer_ss_100jup: PASS (T=%.0f days, %d samples vs REBOUND)\n",
           T_END_DAYS, n_samples);
    return 0;
}
