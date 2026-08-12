/* chaff.h — chaff cloud RCS statistics.
 *
 * A chaff cloud is a collection of N dipole strips (each cut to
 * λ/2 at the victim radar's frequency). Each strip has a
 * resonant-dipole RCS that, when summed incoherently over random
 * orientations, gives a cloud RCS that follows χ²(2 dof)
 * statistics around its mean (Knott 2004 §6.3 for the dipole
 * model; the χ² distribution arises from the I^2 + Q^2 sum over
 * uncorrelated returns).
 *
 * The lib exposes the mean cloud RCS + the CDF / quantile of the
 * χ²(2) distribution; consumers writing detection gates against
 * a chaff-protected target sample from this distribution rather
 * than the mean alone.
 *
 * Determinism: the χ² CDF / quantile are deterministic. Per-
 * sample stochastic draws use the world's seeded RNG threaded
 * through the K26CRng pointer.
 *
 * Knott et al. 2004 + Skolnik 2008 anchors.
 */
#ifndef K26ASTRO_SOFTKILL_CHAFF_H
#define K26ASTRO_SOFTKILL_CHAFF_H

#include "k26compute.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mean RCS (m^2) of a cloud of n_strips dipoles, each strip
 * with resonant RCS sigma_dipole_m2. Returns 0 on NULL / invalid. */
double k26astro_chaff_mean_rcs(int n_strips, double sigma_dipole_m2);

/* CDF of the χ²(2) distribution evaluated at x.
 *   F(x) = 1 - exp(-x/2)   for x >= 0
 * The cloud RCS divided by half its mean follows this CDF. */
double k26astro_chaff_chi2_2_cdf(double x);

/* Quantile of χ²(2) at probability p in [0, 1):
 *   x = -2 ln(1 - p) */
double k26astro_chaff_chi2_2_quantile(double p);

/* Draw a single cloud RCS sample. Returns (mean_rcs / 2) * χ²(2)
 * where the χ²(2) draw is 2 × Exp(1) (i.e. -2 ln(U), U ~ U(0,1)).
 * Passing NULL rng returns the mean RCS (deterministic). */
double k26astro_chaff_sample_rcs(int n_strips,
                                  double sigma_dipole_m2,
                                  K26CRng *rng);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_SOFTKILL_CHAFF_H */
