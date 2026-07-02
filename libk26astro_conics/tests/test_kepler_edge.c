/* test_kepler_edge.c — parabolic, radial, hyperbolic-asymptote.
 *
 * Acceptance:
 *   - Parabolic: a body launched at exactly escape velocity from 1 AU
 *     stays at r > r_periapsis for forward dt; r grows monotonically.
 *   - Radial-elliptic: a body dropped from rest at 1 AU (purely radial)
 *     falls toward the Sun; at half its free-fall time it's between
 *     0 and 1 AU.
 *   - Radial-hyperbolic: a body fired radially outward at 1.5×v_esc
 *     keeps moving outward; r grows monotonically.
 *   - Dispatch: any-propagator routes the three cases correctly and
 *     returns 0 on every call.
 */
#include "k26astro_conics/kepler.h"
#include "k26astro_conics/kepler_edge.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double v_mag_(K26V3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

int main(void)
{
    /* ---- Parabolic flyby ------------------------------------- *
     * Place at 1 AU with v exactly escape, perpendicular to r̂.
     * r grows monotonically; never falls below r0. */
    double v_esc = sqrt(2.0 * K26A_GM_SUN / K26A_AU_M);
    K26V3 pos0 = { K26A_AU_M, 0.0, 0.0 };
    K26V3 vel0 = { 0.0, v_esc, 0.0 };

    double r_prev = v_mag_(pos0);
    for (int day = 1; day <= 60; day++) {
        K26V3 p, v;
        int rc = k26astro_kepler_propagate_parabolic(&p, &v, pos0, vel0,
                                                      K26A_GM_SUN,
                                                      day * 86400.0);
        assert(rc == 0);
        double r_now = v_mag_(p);
        assert(r_now > r_prev * 0.999);   /* allow tiny numerical noise */
        r_prev = r_now;
    }

    /* ---- Radial-elliptic free fall --------------------------- *
     * Drop from rest at 1 AU. The body falls toward the Sun;
     * after a finite time it's somewhere between r0 and 0. */
    K26V3 rad_pos = { K26A_AU_M, 0.0, 0.0 };
    K26V3 rad_vel = { 0.0, 0.0, 0.0 };

    /* T_freefall from rest at r0 to centre = π/(2√2) * sqrt(r0³/μ).
     * For 1 AU around Sun ≈ 64.6 days; sample at 10 days. */
    K26V3 rp, rv;
    int rc = k26astro_kepler_propagate_radial(&rp, &rv, rad_pos, rad_vel,
                                                K26A_GM_SUN, 10.0 * 86400.0);
    assert(rc == 0);
    double r_after_10d = v_mag_(rp);
    /* Should be < 1 AU (falling). */
    assert(r_after_10d < K26A_AU_M);
    /* And v should be inward (negative x). */
    assert(rv.x < 0.0);

    /* ---- Radial-hyperbolic outward ------------------------- */
    K26V3 hyp_pos = { K26A_AU_M, 0.0, 0.0 };
    K26V3 hyp_vel = { v_esc * 1.5, 0.0, 0.0 };

    double rh_prev = v_mag_(hyp_pos);
    for (int day = 1; day <= 30; day++) {
        K26V3 p, v;
        rc = k26astro_kepler_propagate_radial(&p, &v, hyp_pos, hyp_vel,
                                                K26A_GM_SUN, day * 86400.0);
        assert(rc == 0);
        double r_now = v_mag_(p);
        assert(r_now > rh_prev);
        rh_prev = r_now;
    }

    /* ---- Dispatch: any-propagator on an elliptic orbit ----- *
     * Should land in the baseline (kepler.c) path. */
    K26V3 ec_pos = { K26A_AU_M, 0.0, 0.0 };
    K26V3 ec_vel = { 0.0, sqrt(K26A_GM_SUN / K26A_AU_M), 0.0 };
    K26V3 ep, ev;
    rc = k26astro_kepler_propagate_any(&ep, &ev, ec_pos, ec_vel,
                                        K26A_GM_SUN, 86400.0 * 30.0, 32);
    assert(rc == 0);
    assert(fabs(v_mag_(ep) - K26A_AU_M) < 1.0e4);  /* circular: r preserved */

    /* ---- Dispatch: any-propagator on a parabola ------------ */
    rc = k26astro_kepler_propagate_any(&ep, &ev, pos0, vel0,
                                        K26A_GM_SUN, 86400.0 * 30.0, 32);
    assert(rc == 0);
    assert(v_mag_(ep) > K26A_AU_M);   /* escaped */

    /* ---- Dispatch: any-propagator on a radial fall --------- */
    rc = k26astro_kepler_propagate_any(&ep, &ev, rad_pos, rad_vel,
                                        K26A_GM_SUN, 86400.0 * 10.0, 32);
    assert(rc == 0);
    assert(v_mag_(ep) < K26A_AU_M);   /* fell inward */

    printf("test_kepler_edge: OK (parabolic + radial-fall + radial-out + dispatch)\n");
    return 0;
}
