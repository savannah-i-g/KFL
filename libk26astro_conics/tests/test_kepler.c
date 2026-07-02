/* test_kepler.c — universal-variable Kepler propagator.
 *
 * Migrated from libk26astro_body/tests/test_elements.c lines 86-112,
 * extended with:
 *   - non-zero eccentricity round-trip
 *   - hyperbolic-orbit round-trip
 *   - elements-API convenience wrapper round-trip
 *   - mid-step monotonicity (radial distance evolves smoothly)
 */
#include "k26astro_conics/kepler.h"
#include "k26astro_body/elements.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double v_mag_(K26V3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

int main(void)
{
    /* ---- Circular orbit: Earth at 1 AU around Sun --------------- */
    K26V3 pos0 = { K26A_AU_M, 0.0, 0.0 };
    double v_circ = sqrt(K26A_GM_SUN / K26A_AU_M);
    K26V3 vel0 = { 0.0, v_circ, 0.0 };
    double period = K26A_TWO_PI * sqrt(K26A_AU_M * K26A_AU_M * K26A_AU_M
                                       / K26A_GM_SUN);

    /* Full-period round-trip: <1 km residual at 1 AU. */
    K26V3 pos1, vel1;
    int rc = k26astro_kepler_propagate(&pos1, &vel1, pos0, vel0,
                                        K26A_GM_SUN, period, 32);
    assert(rc == 0);
    double dx = pos1.x - pos0.x;
    double dy = pos1.y - pos0.y;
    double dz = pos1.z - pos0.z;
    double res = sqrt(dx * dx + dy * dy + dz * dz);
    assert(res < 1000.0);

    /* Quarter-period: pos should be ~ (0, R, 0). */
    rc = k26astro_kepler_propagate(&pos1, &vel1, pos0, vel0,
                                    K26A_GM_SUN, period * 0.25, 32);
    assert(rc == 0);
    assert(fabs(pos1.x)           < 1.0e5);
    assert(fabs(pos1.y - K26A_AU_M) < 1.0e5);
    assert(fabs(pos1.z)           < 1.0e5);

    /* ---- Elliptic e=0.5 orbit: aphelion → perihelion ----------- */
    /* a = 1 AU, e = 0.5. Perihelion at a(1-e) = 0.5 AU, vp =
     * sqrt(mu(2/rp - 1/a)). Place body at perihelion, propagate half
     * a period, expect to be at aphelion (-1.5 AU on x). */
    double a_e = K26A_AU_M;
    double e_e = 0.5;
    double rp  = a_e * (1.0 - e_e);
    double vp  = sqrt(K26A_GM_SUN * (2.0 / rp - 1.0 / a_e));
    K26V3 pos_e0 = { rp, 0.0, 0.0 };
    K26V3 vel_e0 = { 0.0, vp, 0.0 };
    double period_e = K26A_TWO_PI * sqrt(a_e * a_e * a_e / K26A_GM_SUN);

    K26V3 pos_eh, vel_eh;
    rc = k26astro_kepler_propagate(&pos_eh, &vel_eh, pos_e0, vel_e0,
                                    K26A_GM_SUN, period_e * 0.5, 64);
    assert(rc == 0);
    double r_aphelion = a_e * (1.0 + e_e);
    /* x should be negative aphelion (within ~10 km). */
    assert(fabs(pos_eh.x + r_aphelion) < 1.0e4);
    assert(fabs(pos_eh.y)              < 1.0e4);

    /* Full-period round-trip on elliptic orbit. */
    K26V3 pos_ef, vel_ef;
    rc = k26astro_kepler_propagate(&pos_ef, &vel_ef, pos_e0, vel_e0,
                                    K26A_GM_SUN, period_e, 64);
    assert(rc == 0);
    double res_e = sqrt((pos_ef.x - pos_e0.x) * (pos_ef.x - pos_e0.x)
                      + (pos_ef.y - pos_e0.y) * (pos_ef.y - pos_e0.y)
                      + (pos_ef.z - pos_e0.z) * (pos_ef.z - pos_e0.z));
    assert(res_e < 1000.0);

    /* ---- Hyperbolic flyby ------------------------------------- *
     * Start at 1 AU heading away from Sun at 1.5x escape velocity.
     * Propagate 30 days; r should monotonically grow. */
    double v_esc = sqrt(2.0 * K26A_GM_SUN / K26A_AU_M);
    K26V3 pos_h = { K26A_AU_M, 0.0, 0.0 };
    K26V3 vel_h = { 0.0, v_esc * 1.5, 0.0 };
    double r_prev = v_mag_(pos_h);
    for (int step = 1; step <= 30; step++) {
        K26V3 ph, vh;
        rc = k26astro_kepler_propagate(&ph, &vh, pos_h, vel_h,
                                        K26A_GM_SUN, step * 86400.0, 64);
        assert(rc == 0);
        double r_now = v_mag_(ph);
        assert(r_now > r_prev);
        r_prev = r_now;
    }

    /* ---- Element-API convenience wrapper round-trip ----------- */
    K26AstroKeplerian k0;
    k0.a    = K26A_AU_M * 1.2;
    k0.e    = 0.05;
    k0.i    = 0.1;
    k0.raan = 0.5;
    k0.argp = 1.0;
    k0.M    = 0.3;
    k0.t0   = k26astro_epoch_j2000_tt();
    k0.mu   = K26A_GM_SUN;

    K26AstroKeplerian k_after;
    double dt_quarter = K26A_TWO_PI * sqrt(k0.a*k0.a*k0.a / k0.mu) * 0.25;
    rc = k26astro_kepler_propagate_elements(&k_after, &k0, dt_quarter);
    assert(rc == 0);
    /* a, e, i should be preserved; M should advance by π/2. */
    assert(fabs(k_after.a - k0.a) / k0.a < 1e-9);
    assert(fabs(k_after.e - k0.e)        < 1e-9);
    assert(fabs(k_after.i - k0.i)        < 1e-9);
    /* M_after - M_before mod 2π should be ~π/2. */
    double dM = fmod(k_after.M - k0.M, K26A_TWO_PI);
    if (dM < 0) dM += K26A_TWO_PI;
    assert(fabs(dM - K26A_HALF_PI) < 1e-6);

    /* ---- Zero-dt no-op --------------------------------------- */
    K26V3 pzero, vzero;
    rc = k26astro_kepler_propagate(&pzero, &vzero, pos0, vel0,
                                    K26A_GM_SUN, 0.0, 32);
    assert(rc == 0);
    assert(pzero.x == pos0.x && pzero.y == pos0.y && pzero.z == pos0.z);
    assert(vzero.x == vel0.x && vzero.y == vel0.y && vzero.z == vel0.z);

    /* ---- Reverse propagation ---------------------------------- *
     * Propagate forward by dt then back by -dt; expect round-trip
     * to within ~mm. */
    K26V3 pa, va, pb, vb;
    rc = k26astro_kepler_propagate(&pa, &va, pos0, vel0,
                                    K26A_GM_SUN, period * 0.3, 64);
    assert(rc == 0);
    rc = k26astro_kepler_propagate(&pb, &vb, pa, va,
                                    K26A_GM_SUN, -period * 0.3, 64);
    assert(rc == 0);
    double rev_res = sqrt((pb.x - pos0.x) * (pb.x - pos0.x)
                        + (pb.y - pos0.y) * (pb.y - pos0.y)
                        + (pb.z - pos0.z) * (pb.z - pos0.z));
    assert(rev_res < 100.0);  /* <100 m round-trip over Earth orbit */

    printf("test_kepler: OK (circ + ellip + hyperbolic + elements + reverse)\n");
    return 0;
}
