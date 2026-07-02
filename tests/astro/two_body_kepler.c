/* tests/astro/two_body_kepler.c - Two-body Kepler propagator gate.
 *
 * Setup: Sun + Earth, Keplerian initial conditions. Propagate 1000
 * orbits via libk26astro_conics' universal-variable propagator.
 *
 * Acceptance:
 *   - Energy conservation:           |dE/E|        < 1e-12
 *   - Angular momentum conservation: |dL/L|        < 1e-12
 *   - Orbital element drift:                         < 1e-10 per orbit
 *
 * The Kepler propagator is exact for two-body; non-conservation here
 * indicates a bug in the universal-variable solver. */
#include "k26astro_conics/kepler.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static double v_dot_(K26V3 a, K26V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double v_norm_(K26V3 a)         { return sqrt(v_dot_(a, a)); }
static K26V3  v_cross_(K26V3 a, K26V3 b)
{
    K26V3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}

int main(void)
{
    /* Initial state: Earth-like at periapsis. */
    double a0 = K26A_AU_M;
    double e0 = 0.0167;                     /* Earth's eccentricity */
    double rp = a0 * (1.0 - e0);
    double vp = sqrt(K26A_GM_SUN * (2.0 / rp - 1.0 / a0));
    K26V3 pos = { rp, 0.0, 0.0 };
    K26V3 vel = { 0.0, vp, 0.0 };

    double mu = K26A_GM_SUN;
    double period = K26A_TWO_PI * sqrt(a0*a0*a0 / mu);

    /* Conserved quantities at t=0. */
    double r0_mag = v_norm_(pos);
    double v0_mag = v_norm_(vel);
    double E0 = 0.5 * v0_mag * v0_mag - mu / r0_mag;
    K26V3  L0 = v_cross_(pos, vel);
    double L0_mag = v_norm_(L0);

    /* 1000-orbit run. Propagate one orbit at a time to keep χ
     * convergence well-conditioned. */
    int n_orbits = 1000;
    double dE_max = 0.0;
    double dL_max = 0.0;
    K26V3  pos_n = pos, vel_n = vel;
    for (int k = 0; k < n_orbits; k++) {
        K26V3 pos_next, vel_next;
        int rc = k26astro_kepler_propagate(&pos_next, &vel_next,
                                            pos_n, vel_n,
                                            mu, period, 64);
        assert(rc == 0);
        pos_n = pos_next;
        vel_n = vel_next;

        double rn = v_norm_(pos_n);
        double vn = v_norm_(vel_n);
        double En = 0.5 * vn * vn - mu / rn;
        K26V3  Ln = v_cross_(pos_n, vel_n);
        double Ln_mag = v_norm_(Ln);

        double dE = fabs((En - E0) / E0);
        double dL = fabs((Ln_mag - L0_mag) / L0_mag);
        if (dE > dE_max) dE_max = dE;
        if (dL > dL_max) dL_max = dL;
    }

    /* Acceptance: 1e-12 ceiling. Universal-variable Kepler is exact
     * in IEEE-754 doubles for non-degenerate orbits; observed drift
     * comes from Newton-Raphson convergence threshold + accumulated
     * round-off in the Lagrange f/g formula (~10 ULPs per step). */
    printf("two_body_kepler: %d orbits, dE_max=%.3e, dL_max=%.3e\n",
           n_orbits, dE_max, dL_max);
    assert(dE_max < 1.0e-12);
    assert(dL_max < 1.0e-12);

    /* Element drift: after N orbits the periapsis distance should
     * have drifted by < 1e-10 · N · a0. */
    double rp_final = v_norm_(pos_n);
    double rp_drift = fabs(rp_final - rp);
    double rp_drift_per_orbit = rp_drift / n_orbits;
    printf("                         rp_drift/orbit=%.3e m\n",
           rp_drift_per_orbit);
    assert(rp_drift_per_orbit < 1.0e-10 * a0);

    printf("two_body_kepler: PASS\n");
    return 0;
}
