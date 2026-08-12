/* chaff.c — dipole-strip chaff cloud RCS statistics.
 *
 * Mean cloud RCS: σ_cloud = N · σ_dipole (sum of incoherent
 * dipole contributions). Individual measurements fluctuate around
 * the mean with χ²(2 dof) statistics (Skolnik 2008 §15.3 noting
 * the Rayleigh-amplitude / exponential-power equivalence).
 *
 * χ²(2 dof) is equivalent to 2 × Exp(1): if U ~ Uniform(0, 1),
 * then -2 ln(U) ~ χ²(2). */

#include "k26astro_softkill/chaff.h"

#include <math.h>

double k26astro_chaff_mean_rcs(int n_strips, double sigma_dipole_m2)
{
    if (n_strips <= 0 || sigma_dipole_m2 <= 0.0) return 0.0;
    return (double)n_strips * sigma_dipole_m2;
}

double k26astro_chaff_chi2_2_cdf(double x)
{
    if (x <= 0.0) return 0.0;
    return 1.0 - exp(-0.5 * x);
}

double k26astro_chaff_chi2_2_quantile(double p)
{
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return INFINITY;
    return -2.0 * log(1.0 - p);
}

double k26astro_chaff_sample_rcs(int n_strips,
                                  double sigma_dipole_m2,
                                  K26CRng *rng)
{
    const double mean = k26astro_chaff_mean_rcs(n_strips, sigma_dipole_m2);
    if (!rng) return mean;
    const double u = k26c_rng_uniform(rng);
    /* Avoid log(0); clamp the uniform draw to a tiny positive
     * value so the χ² sample is finite. */
    const double u_safe = (u <= 0.0) ? 1.0e-300 : u;
    const double chi2 = -2.0 * log(u_safe);
    /* Sample = (mean / 2) · χ²(2). Mean of (mean/2)·Exp(2)
     * equals mean. */
    return 0.5 * mean * chi2;
}
