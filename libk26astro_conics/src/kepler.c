/* kepler.c — universal-variable Kepler propagator.
 *
 * The Newton-Raphson iterates on χ until
 *   F(χ) = r0·σ0/√μ · χ²·C(αχ²) +
 *          (1 - r0·α) · χ³·S(αχ²) + r0·χ - √μ·dt ≈ 0,
 * where α = 1/a, σ0 = r⃗·v⃗/√μ. The conic type is implicit in α (>0
 * elliptic, <0 hyperbolic, 0 parabolic; see kepler_edge.c).
 *
 * Reference: Vallado §2.4, Algorithm 8. */
#include "k26astro_conics/kepler.h"
#include "stumpff_internal.h"

#include <math.h>

/* ---- Local V3 helpers (deliberately duplicated from body to keep the
 *      conics inner loop self-contained; promoting them to libk26m3d's
 *      public surface is a future housekeeping task). */
static double v_dot_(K26V3 a, K26V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double v_norm_(K26V3 a)         { return sqrt(v_dot_(a, a)); }
static K26V3  v_scale_(K26V3 a, double s)
{
    K26V3 r = { a.x * s, a.y * s, a.z * s };
    return r;
}
static K26V3 v_add_(K26V3 a, K26V3 b)
{
    K26V3 r = { a.x + b.x, a.y + b.y, a.z + b.z };
    return r;
}

int k26astro_kepler_propagate(K26V3 *out_pos, K26V3 *out_vel,
                              K26V3 pos0, K26V3 vel0,
                              double mu, double dt,
                              int max_iter)
{
    if (!out_pos || !out_vel) return 1;
    if (mu <= 0.0) return 2;
    if (dt == 0.0) { *out_pos = pos0; *out_vel = vel0; return 0; }

    double r0 = v_norm_(pos0);
    if (r0 <= 0.0) return 3;
    double v0 = v_norm_(vel0);
    double sqrt_mu = sqrt(mu);
    double alpha = 2.0 / r0 - v0 * v0 / mu;
    double sigma0 = v_dot_(pos0, vel0) / sqrt_mu;

    double chi = sqrt_mu * fabs(alpha) * dt;
    if (alpha > 0.0) {
        chi = sqrt_mu * alpha * dt;
    } else if (alpha < 0.0) {
        double sign_dt = (dt > 0.0) ? 1.0 : -1.0;
        double mag = sqrt(-1.0 / alpha)
                   * log((-2.0 * mu * alpha * dt)
                         / (sigma0 + sign_dt * sqrt(-mu / alpha) * (1.0 - r0 * alpha)));
        chi = sign_dt * mag;
        if (!isfinite(chi)) chi = sqrt_mu * dt / r0;
    } else {
        chi = sqrt_mu * dt / r0;
    }

    if (max_iter <= 0) max_iter = 32;
    double C = 0.0, S = 0.0, r = r0, psi = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        psi = chi * chi * alpha;
        C = k26astro_conics_stumpff_C(psi);
        S = k26astro_conics_stumpff_S(psi);
        double F = sigma0 * chi * chi * C
                 + (1.0 - r0 * alpha) * chi * chi * chi * S
                 + r0 * chi
                 - sqrt_mu * dt;
        double dF = sigma0 * chi * (1.0 - psi * S)
                  + (1.0 - r0 * alpha) * chi * chi * C
                  + r0;
        double d = F / dF;
        chi -= d;
        if (fabs(d) < 1.0e-12) break;
    }

    double f = 1.0 - chi * chi / r0 * C;
    double g = dt - chi * chi * chi / sqrt_mu * S;
    r = sigma0 * chi * (1.0 - psi * S)
      + (1.0 - r0 * alpha) * chi * chi * C
      + r0;
    if (r <= 0.0) return 4;
    double fdot = sqrt_mu / (r * r0) * (psi * S - 1.0) * chi;
    double gdot = 1.0 - chi * chi / r * C;

    *out_pos = v_add_(v_scale_(pos0, f),    v_scale_(vel0, g));
    *out_vel = v_add_(v_scale_(pos0, fdot), v_scale_(vel0, gdot));
    return 0;
}

/* ---- Convenience: propagate Keplerian elements forward by dt ------ */
int k26astro_kepler_propagate_elements(K26AstroKeplerian *out,
                                       const K26AstroKeplerian *in,
                                       double dt)
{
    if (!out || !in) return 1;

    /* For elliptic orbits the cleanest path is M += n·dt where
     * n = sqrt(mu / a³). For hyperbolic we round-trip through state
     * to avoid mean-anomaly sign confusion. The state-vector path
     * works for both. */
    K26AstroPos central = k26astro_pos_zero();
    K26AstroStateVector sv;
    int rc = k26astro_state_from_elements(&sv, in, &central);
    if (rc != 0) return rc;

    K26V3 pos0 = k26astro_pos_sub(&sv.pos, &central);
    K26V3 pos1, vel1;
    rc = k26astro_kepler_propagate(&pos1, &vel1, pos0, sv.vel,
                                    in->mu, dt, 32);
    if (rc != 0) return rc;

    sv.pos = central;
    k26astro_pos_add(&sv.pos, pos1);
    sv.vel = vel1;
    k26astro_epoch_add_seconds(&sv.t0, dt);

    return k26astro_elements_from_state(out, &sv, &central);
}
