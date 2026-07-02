/* test_elements.c — Keplerian + equinoctial round-trip tests.
 *
 * The universal-variable propagator's round-trip assertions migrated
 * to libk26astro_conics/tests/test_kepler.c on 2026-05-22 along with
 * the propagator source. */
#include "k26astro_body/elements.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void)
{
    /* ---- Major-body lookup ------------------------------------ */
    const K26AstroMajorBody *earth = k26astro_major_body_find("earth");
    assert(earth != NULL);
    assert(earth->naif_id == 399);
    assert(near_(earth->gm, K26A_GM_EARTH, 1e-12));

    /* ---- Keplerian → state → Keplerian round-trip ----------- *
     * Earth's orbit: a = 1 AU, e = 0.0167, i = ~0 (ecliptic). */
    K26AstroKeplerian k0;
    k0.a    = K26A_AU_M;
    k0.e    = 0.0167;
    k0.i    = 0.0;            /* placed in equatorial for simplicity */
    k0.raan = 0.0;
    k0.argp = 0.0;
    k0.M    = K26A_HALF_PI;   /* quarter-orbit past periapsis */
    k0.t0   = k26astro_epoch_j2000_tt();
    k0.mu   = K26A_GM_SUN;

    K26AstroPos sun_pos = k26astro_pos_zero();
    K26AstroStateVector s;
    int rc = k26astro_state_from_elements(&s, &k0, &sun_pos);
    assert(rc == 0);

    K26AstroKeplerian k1;
    rc = k26astro_elements_from_state(&k1, &s, &sun_pos);
    assert(rc == 0);
    /* a, e should round-trip cleanly. */
    assert(near_(k1.a, k0.a, 1e-9));
    assert(near_(k1.e, k0.e, 1e-9));
    /* M can wrap differently when i = 0 (raan / argp undefined); we
     * check M + raan + argp = M_total which is well-defined. */
    double sum0 = k0.M + k0.raan + k0.argp;
    double sum1 = k1.M + k1.raan + k1.argp;
    double diff = fmod(sum1 - sum0, K26A_TWO_PI);
    if (diff > K26A_PI)  diff -= K26A_TWO_PI;
    if (diff < -K26A_PI) diff += K26A_TWO_PI;
    assert(fabs(diff) < 1e-8);

    /* ---- Keplerian ↔ equinoctial round-trip ----------------- *
     * Use a non-zero eccentricity + non-zero inclination to avoid
     * the degenerate cases. */
    K26AstroKeplerian kk;
    kk.a = K26A_AU_M * 1.5;
    kk.e = 0.1;
    kk.i = 0.05;
    kk.raan = 1.0;
    kk.argp = 2.0;
    kk.M    = 0.5;
    kk.t0   = k26astro_epoch_j2000_tt();
    kk.mu   = K26A_GM_SUN;

    K26AstroEquinoctial eq;
    k26astro_equinoctial_from_keplerian(&eq, &kk);
    K26AstroKeplerian kk2;
    k26astro_keplerian_from_equinoctial(&kk2, &eq);
    assert(near_(kk2.a, kk.a, 1e-12));
    assert(near_(kk2.e, kk.e, 1e-9));
    assert(near_(kk2.i, kk.i, 1e-9));

    /* ---- Anomaly conversions ------------------------------ */
    /* M=0 → ν=0, M=π → ν=π (for any e < 1). */
    for (double e = 0.0; e < 0.9; e += 0.1) {
        double nu0 = k26astro_anomaly_mean_to_true(0.0, e);
        assert(fabs(nu0) < 1e-9);
        double nu_pi = k26astro_anomaly_mean_to_true(K26A_PI, e);
        assert(fabs(nu_pi - K26A_PI) < 1e-9);
    }

    printf("test_elements: OK\n");
    return 0;
}
