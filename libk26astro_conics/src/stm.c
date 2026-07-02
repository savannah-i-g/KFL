/* stm.c — 6 x 6 Cartesian state-transition matrix via centered
 * finite differences through the universal-variable Kepler
 * propagator.
 *
 * Battin (1999) paragraph 10.5; Vallado (2013) paragraph 10.7.2. */

#include "k26astro_conics/stm.h"
#include "k26astro_conics/kepler.h"

#include <math.h>

#define KEPLER_NEWTON_MAX_ITER 32

/* Centered FD step-size constants. The position step is taken
 * relative to |r0| and the velocity step relative to |v0|; the
 * geometric mean of FLT_EPSILON and DBL_EPSILON sets the relative
 * floor at roughly 1e-6, which balances roundoff vs truncation for
 * the centered O(h^2) stencil. */
#define K26_STM_REL_STEP        1.0e-6
#define K26_STM_R_ABS_FLOOR_M   1.0e-3
#define K26_STM_V_ABS_FLOOR_MS  1.0e-6

int k26astro_kepler_stm(K26V3 r0, K26V3 v0, double dt, double mu,
                        double M_out[6][6])
{
    if (!M_out) return 1;
    if (!(mu > 0.0)) return 1;

    double r0_norm = sqrt(r0.x * r0.x + r0.y * r0.y + r0.z * r0.z);
    double v0_norm = sqrt(v0.x * v0.x + v0.y * v0.y + v0.z * v0.z);
    if (!(r0_norm > 0.0)) return 1;

    double h_r = K26_STM_REL_STEP * r0_norm;
    double h_v = K26_STM_REL_STEP * ((v0_norm > 0.0) ? v0_norm : 1.0);
    if (h_r < K26_STM_R_ABS_FLOOR_M)  h_r = K26_STM_R_ABS_FLOOR_M;
    if (h_v < K26_STM_V_ABS_FLOOR_MS) h_v = K26_STM_V_ABS_FLOOR_MS;

    double s0[6] = { r0.x, r0.y, r0.z, v0.x, v0.y, v0.z };
    double hs[6] = { h_r, h_r, h_r, h_v, h_v, h_v };

    for (int j = 0; j < 6; j++) {
        double sp[6], sm[6];
        for (int k = 0; k < 6; k++) {
            sp[k] = s0[k];
            sm[k] = s0[k];
        }
        sp[j] += hs[j];
        sm[j] -= hs[j];

        K26V3 r_p  = { sp[0], sp[1], sp[2] };
        K26V3 v_p  = { sp[3], sp[4], sp[5] };
        K26V3 r_m  = { sm[0], sm[1], sm[2] };
        K26V3 v_m  = { sm[3], sm[4], sm[5] };

        K26V3 r_pos_p, v_pos_p;
        K26V3 r_pos_m, v_pos_m;
        if (k26astro_kepler_propagate(&r_pos_p, &v_pos_p,
                                       r_p, v_p, mu, dt,
                                       KEPLER_NEWTON_MAX_ITER) != 0) {
            return 2;
        }
        if (k26astro_kepler_propagate(&r_pos_m, &v_pos_m,
                                       r_m, v_m, mu, dt,
                                       KEPLER_NEWTON_MAX_ITER) != 0) {
            return 2;
        }

        double inv2h = 0.5 / hs[j];
        M_out[0][j] = (r_pos_p.x - r_pos_m.x) * inv2h;
        M_out[1][j] = (r_pos_p.y - r_pos_m.y) * inv2h;
        M_out[2][j] = (r_pos_p.z - r_pos_m.z) * inv2h;
        M_out[3][j] = (v_pos_p.x - v_pos_m.x) * inv2h;
        M_out[4][j] = (v_pos_p.y - v_pos_m.y) * inv2h;
        M_out[5][j] = (v_pos_p.z - v_pos_m.z) * inv2h;
    }

    return 0;
}
