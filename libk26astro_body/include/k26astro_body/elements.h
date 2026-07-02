/* k26astro_body/elements.h — orbital elements.
 *
 * Two representations:
 *
 *   Keplerian        (a, e, i, raan, argp, M)
 *     The textbook set. Used by catalogues (Minor Planet Center
 *     publishes orbits in this form), human-readable orbit
 *     specifications, and the on-rails patched-conic mode.
 *     Singularities at e = 0 (raan/argp undefined) and i = 0
 *     (raan undefined).
 *
 *   Equinoctial      (a, h, k, p, q, lambda)
 *     Non-singular form derived by Broucke & Cefola (1972). Used
 *     by perturbation theory and any context where the Keplerian
 *     singularities would cause numerical trouble (low-eccentricity
 *     orbits, low-inclination orbits, transfers through e ≈ 0).
 *
 * Conversion in both directions is provided; the round-trip is
 * exact in IEEE-754 doubles outside the singular cases.
 *
 * The universal-variable Kepler propagator that used to live here
 * migrated to libk26astro_conics/include/k26astro_conics/kepler.h
 * on 2026-05-22; include that header for `k26astro_kepler_propagate`. */
#ifndef K26ASTRO_BODY_ELEMENTS_H
#define K26ASTRO_BODY_ELEMENTS_H

#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Keplerian -------------------------------------------------- */

typedef struct {
    double a;            /* semi-major axis, m (negative for hyperbolic) */
    double e;            /* eccentricity */
    double i;            /* inclination, rad */
    double raan;         /* right ascension of ascending node, rad */
    double argp;         /* argument of periapsis, rad */
    double M;            /* mean anomaly at epoch, rad */
    K26AstroEpoch t0;    /* epoch the elements are valid at */
    double mu;           /* GM of the central body, m^3/s^2 */
} K26AstroKeplerian;

/* ---- Equinoctial ------------------------------------------------ */

typedef struct {
    double a;
    double h;            /* e sin(argp + raan) */
    double k;            /* e cos(argp + raan) */
    double p;            /* tan(i/2) sin(raan) */
    double q;            /* tan(i/2) cos(raan) */
    double lambda;       /* mean longitude = M + argp + raan, rad */
    K26AstroEpoch t0;
    double mu;
} K26AstroEquinoctial;

/* ---- Cartesian state vector ----------------------------------- */
typedef struct {
    K26AstroPos pos;
    K26V3       vel;
    K26AstroEpoch t0;
    double mu;
} K26AstroStateVector;

/* ---- Conversions ---------------------------------------------- */

/* State → Keplerian. Returns 0 on success, non-zero if the
 * relative position to the central body is degenerate (r ≤ 0) or
 * the orbit is parabolic (e ≈ 1) and the caller hasn't supplied
 * elements via the equinoctial form. */
int k26astro_elements_from_state(K26AstroKeplerian *out,
                                  const K26AstroStateVector *s,
                                  const K26AstroPos *central_pos);

/* Keplerian → state. The returned position is in the same inertial
 * frame as `central_pos`; the caller adds `central_pos` to the
 * returned offset to get an absolute K26AstroPos. */
int k26astro_state_from_elements(K26AstroStateVector *out,
                                  const K26AstroKeplerian *k,
                                  const K26AstroPos *central_pos);

/* Keplerian ↔ equinoctial. Round-trip is exact away from the
 * Keplerian singularities. */
void k26astro_equinoctial_from_keplerian(K26AstroEquinoctial *out,
                                          const K26AstroKeplerian *k);
void k26astro_keplerian_from_equinoctial(K26AstroKeplerian *out,
                                          const K26AstroEquinoctial *eq);

/* ---- Mean ↔ true anomaly ------------------------------------- */
/* Standard Newton iteration for E from M = E - e sin E (elliptic);
 * Halley iteration for hyperbolic. Returns the true anomaly ν. */
double k26astro_anomaly_mean_to_true(double M, double e);
double k26astro_anomaly_true_to_mean(double nu, double e);

/* Orbital period for elliptic orbits. Returns infinity for e ≥ 1. */
double k26astro_period(const K26AstroKeplerian *k);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_BODY_ELEMENTS_H */
