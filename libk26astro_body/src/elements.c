/* elements.c — Keplerian + equinoctial elements, state conversions.
 *
 * The universal-variable Kepler propagator + Stumpff helpers live in
 * libk26astro_conics/src/kepler.c (and src/stumpff.c). body keeps the
 * elements-shaped state-and-conversion machinery; conics owns
 * propagation.
 *
 * References:
 *   - Bate, Mueller & White, Fundamentals of Astrodynamics (1971),
 *     §2.5 (state → elements)
 *   - Broucke & Cefola (1972), "On the Equinoctial Orbit Elements",
 *     CeMec 5:303 — Keplerian↔equinoctial conversion
 *
 * Conventions:
 *   - Right-handed inertial frame (ICRF-style)
 *   - Inclination measured from +z
 *   - Argument of periapsis measured from ascending node, in
 *     direction of motion
 *   - Mean anomaly: M = n (t - t_peri), n = sqrt(mu / |a|^3)
 *   - Negative `a` indicates hyperbolic orbit (e > 1)
 */
#include "k26astro_body/elements.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

/* ---- Vector helpers ------------------------------------------- */

static double v_dot_(K26V3 a, K26V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static K26V3  v_cross_(K26V3 a, K26V3 b) {
    K26V3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}
static double v_norm_(K26V3 a) { return sqrt(v_dot_(a, a)); }
static K26V3  v_scale_(K26V3 a, double s) {
    K26V3 r = { a.x * s, a.y * s, a.z * s };
    return r;
}
static K26V3  v_sub_(K26V3 a, K26V3 b) {
    K26V3 r = { a.x - b.x, a.y - b.y, a.z - b.z };
    return r;
}

/* Wrap angle into [0, 2π). */
static double wrap_2pi_(double a)
{
    double two_pi = K26A_TWO_PI;
    a = fmod(a, two_pi);
    if (a < 0.0) a += two_pi;
    return a;
}

/* ---- State → Keplerian ---------------------------------------- */

int k26astro_elements_from_state(K26AstroKeplerian *out,
                                  const K26AstroStateVector *s,
                                  const K26AstroPos *central_pos)
{
    if (!out || !s) return 1;

    /* Position relative to central body. */
    K26V3 r;
    if (central_pos) {
        r = k26astro_pos_sub(&s->pos, central_pos);
    } else {
        r = k26astro_pos_to_m_approx(&s->pos);
    }
    K26V3 v   = s->vel;
    double mu = s->mu;

    double r_mag = v_norm_(r);
    if (r_mag <= 0.0) return 2;
    double v2    = v_dot_(v, v);

    /* Specific orbital energy & angular momentum. */
    double energy = 0.5 * v2 - mu / r_mag;
    K26V3  h      = v_cross_(r, v);
    double h_mag  = v_norm_(h);

    /* Semi-major axis. */
    double a;
    if (fabs(energy) < 1.0e-12) {
        /* Parabolic — caller should use the equinoctial form. */
        return 3;
    }
    a = -mu / (2.0 * energy);

    /* Eccentricity vector. */
    K26V3 e_vec;
    {
        K26V3 t1 = v_cross_(v, h);
        K26V3 t1s = v_scale_(t1, 1.0 / mu);
        K26V3 t2  = v_scale_(r, 1.0 / r_mag);
        e_vec = v_sub_(t1s, t2);
    }
    double e = v_norm_(e_vec);

    /* Inclination. */
    double i = acos(h.z / h_mag);
    if (i < 0.0) i = 0.0;
    if (i > K26A_PI) i = K26A_PI;

    /* Node vector. */
    K26V3 n = { -h.y, h.x, 0.0 };
    double n_mag = v_norm_(n);

    /* Right ascension of ascending node. */
    double raan;
    if (n_mag > 1.0e-12) {
        raan = acos(n.x / n_mag);
        if (n.y < 0.0) raan = K26A_TWO_PI - raan;
    } else {
        raan = 0.0;   /* equatorial orbit — convention: node at +x */
    }

    /* Argument of periapsis. */
    double argp;
    if (n_mag > 1.0e-12 && e > 1.0e-12) {
        argp = acos(v_dot_(n, e_vec) / (n_mag * e));
        if (e_vec.z < 0.0) argp = K26A_TWO_PI - argp;
    } else if (e > 1.0e-12) {
        /* Equatorial, eccentric — measure argp from +x. */
        argp = atan2(e_vec.y, e_vec.x);
        argp = wrap_2pi_(argp);
    } else {
        argp = 0.0;
    }

    /* True anomaly + mean anomaly. */
    double nu;
    if (e > 1.0e-12) {
        double cnv = v_dot_(e_vec, r) / (e * r_mag);
        if (cnv >  1.0) cnv =  1.0;
        if (cnv < -1.0) cnv = -1.0;
        nu = acos(cnv);
        if (v_dot_(r, v) < 0.0) nu = K26A_TWO_PI - nu;
    } else {
        /* Circular — measure ν from node (or +x for equatorial). */
        if (n_mag > 1.0e-12) {
            double cnu = v_dot_(n, r) / (n_mag * r_mag);
            if (cnu >  1.0) cnu =  1.0;
            if (cnu < -1.0) cnu = -1.0;
            nu = acos(cnu);
            if (r.z < 0.0) nu = K26A_TWO_PI - nu;
        } else {
            nu = atan2(r.y, r.x);
            nu = wrap_2pi_(nu);
        }
    }

    double M = k26astro_anomaly_true_to_mean(nu, e);

    out->a    = a;
    out->e    = e;
    out->i    = i;
    out->raan = wrap_2pi_(raan);
    out->argp = wrap_2pi_(argp);
    out->M    = wrap_2pi_(M);
    out->t0   = s->t0;
    out->mu   = mu;
    return 0;
}

/* ---- Keplerian → state ---------------------------------------- */

int k26astro_state_from_elements(K26AstroStateVector *out,
                                  const K26AstroKeplerian *k,
                                  const K26AstroPos *central_pos)
{
    if (!out || !k) return 1;
    double mu = k->mu;
    if (mu <= 0.0) return 2;

    /* Solve Kepler for E (elliptic) or hyperbolic anomaly. */
    double nu = k26astro_anomaly_mean_to_true(k->M, k->e);

    /* Distance r from true anomaly. */
    double p = k->a * (1.0 - k->e * k->e);
    double r_mag = p / (1.0 + k->e * cos(nu));

    /* Position + velocity in perifocal frame (PQW). */
    double cos_nu = cos(nu), sin_nu = sin(nu);
    K26V3 r_pqw = { r_mag * cos_nu, r_mag * sin_nu, 0.0 };
    double h_pqw_scale = sqrt(mu * p);
    K26V3 v_pqw = {
        -h_pqw_scale / p * sin_nu,
         h_pqw_scale / p * (k->e + cos_nu),
         0.0
    };

    /* Rotate PQW → inertial via R3(-raan) * R1(-i) * R3(-argp). */
    double cR = cos(k->raan), sR = sin(k->raan);
    double cI = cos(k->i),    sI = sin(k->i);
    double cA = cos(k->argp), sA = sin(k->argp);

    double R[3][3];
    R[0][0] =  cR * cA - sR * sA * cI;
    R[0][1] = -cR * sA - sR * cA * cI;
    R[0][2] =  sR * sI;
    R[1][0] =  sR * cA + cR * sA * cI;
    R[1][1] = -sR * sA + cR * cA * cI;
    R[1][2] = -cR * sI;
    R[2][0] =  sA * sI;
    R[2][1] =  cA * sI;
    R[2][2] =  cI;

    K26V3 r_in = {
        R[0][0] * r_pqw.x + R[0][1] * r_pqw.y + R[0][2] * r_pqw.z,
        R[1][0] * r_pqw.x + R[1][1] * r_pqw.y + R[1][2] * r_pqw.z,
        R[2][0] * r_pqw.x + R[2][1] * r_pqw.y + R[2][2] * r_pqw.z
    };
    K26V3 v_in = {
        R[0][0] * v_pqw.x + R[0][1] * v_pqw.y + R[0][2] * v_pqw.z,
        R[1][0] * v_pqw.x + R[1][1] * v_pqw.y + R[1][2] * v_pqw.z,
        R[2][0] * v_pqw.x + R[2][1] * v_pqw.y + R[2][2] * v_pqw.z
    };

    if (central_pos) {
        out->pos = *central_pos;
        k26astro_pos_add(&out->pos, r_in);
    } else {
        out->pos = k26astro_pos_from_m(r_in.x, r_in.y, r_in.z);
    }
    out->vel = v_in;
    out->t0  = k->t0;
    out->mu  = mu;
    return 0;
}

/* ---- Keplerian ↔ Equinoctial ---------------------------------- */

void k26astro_equinoctial_from_keplerian(K26AstroEquinoctial *out,
                                          const K26AstroKeplerian *k)
{
    if (!out || !k) return;
    out->a      = k->a;
    out->h      = k->e * sin(k->argp + k->raan);
    out->k      = k->e * cos(k->argp + k->raan);
    out->p      = tan(k->i * 0.5) * sin(k->raan);
    out->q      = tan(k->i * 0.5) * cos(k->raan);
    out->lambda = wrap_2pi_(k->M + k->argp + k->raan);
    out->t0     = k->t0;
    out->mu     = k->mu;
}

void k26astro_keplerian_from_equinoctial(K26AstroKeplerian *out,
                                          const K26AstroEquinoctial *eq)
{
    if (!out || !eq) return;
    double a    = eq->a;
    double e    = sqrt(eq->h * eq->h + eq->k * eq->k);
    double tan_half_i_2 = eq->p * eq->p + eq->q * eq->q;
    double i;
    if (tan_half_i_2 > 0.0) {
        i = 2.0 * atan(sqrt(tan_half_i_2));
    } else {
        i = 0.0;
    }
    double raan = (tan_half_i_2 > 0.0) ? atan2(eq->p, eq->q) : 0.0;
    double argp_plus_raan = (e > 0.0) ? atan2(eq->h, eq->k) : 0.0;
    double argp = argp_plus_raan - raan;
    double M = eq->lambda - argp_plus_raan;

    out->a    = a;
    out->e    = e;
    out->i    = i;
    out->raan = wrap_2pi_(raan);
    out->argp = wrap_2pi_(argp);
    out->M    = wrap_2pi_(M);
    out->t0   = eq->t0;
    out->mu   = eq->mu;
}

/* ---- Anomaly conversions -------------------------------------- */

double k26astro_anomaly_mean_to_true(double M, double e)
{
    if (e < 1.0) {
        /* Elliptic — Newton on E - e sin E = M. */
        double E = (e < 0.8) ? M : K26A_PI;
        for (int iter = 0; iter < 64; iter++) {
            double f  = E - e * sin(E) - M;
            double fp = 1.0 - e * cos(E);
            double dE = f / fp;
            E -= dE;
            if (fabs(dE) < 1.0e-14) break;
        }
        double sE = sin(E), cE = cos(E);
        double nu = atan2(sqrt(1.0 - e * e) * sE, cE - e);
        return wrap_2pi_(nu);
    } else if (e > 1.0) {
        /* Hyperbolic — M = e sinh H - H. */
        double H = M;
        for (int iter = 0; iter < 64; iter++) {
            double f  = e * sinh(H) - H - M;
            double fp = e * cosh(H) - 1.0;
            double dH = f / fp;
            H -= dH;
            if (fabs(dH) < 1.0e-14) break;
        }
        double sH = sinh(H), cH = cosh(H);
        double nu = atan2(sqrt(e * e - 1.0) * sH, e - cH);
        return wrap_2pi_(nu);
    } else {
        /* Parabolic — Barker's equation. */
        double D = cbrt(M + sqrt(M * M + 1.0))
                 - cbrt(-M + sqrt(M * M + 1.0));
        return wrap_2pi_(2.0 * atan(D));
    }
}

double k26astro_anomaly_true_to_mean(double nu, double e)
{
    if (e < 1.0) {
        double sn = sin(nu), cn = cos(nu);
        double E = atan2(sqrt(1.0 - e * e) * sn, e + cn);
        double M = E - e * sin(E);
        return wrap_2pi_(M);
    } else if (e > 1.0) {
        double sn = sin(nu), cn = cos(nu);
        double H = atan2(sqrt(e * e - 1.0) * sn, e + cn);
        H = log(tan(0.5 * H + 0.25 * K26A_PI)) * 2.0;  /* atanh form for stability */
        double M = e * sinh(H) - H;
        return M;
    } else {
        return tan(0.5 * nu) + tan(0.5 * nu) * tan(0.5 * nu) * tan(0.5 * nu) / 3.0;
    }
}

double k26astro_period(const K26AstroKeplerian *k)
{
    if (!k || k->e >= 1.0) return INFINITY;
    return K26A_TWO_PI * sqrt(k->a * k->a * k->a / k->mu);
}

/* Kepler propagation (universal-variable, Stumpff C/S) migrated to
 * libk26astro_conics on 2026-05-22 — see kepler.c + stumpff.c there. */
