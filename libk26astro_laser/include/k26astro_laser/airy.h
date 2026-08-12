/* airy.h — diffraction-limited beam geometry.
 *
 * Airy first-null half-angle: θ ≈ 1.22 λ / D. Inflated by beam
 * quality M^2 and (Gaussian-RMS) pointing jitter σ_θ:
 *
 *     θ_total^2 = (M^2 · 1.22 λ / D)^2 + σ_θ^2
 *
 * Spot radius at range R: r = R · θ_total.
 *
 * Encircled-energy fraction inside a disc of radius r at range R:
 *
 *     F(x) = 1 - J_0(x)^2 - J_1(x)^2
 *
 * with x = (π D r) / (λ R) and J_n the Bessel functions of the
 * first kind. F → 0.8378 at the first dark ring (x = 3.8317), the
 * canonical Airy-disc value (Born & Wolf 2002, §8.5.2).
 *
 * The lib provides J_0 and J_1 via the standard library (`j0`,
 * `j1`); for x outside the convergence regime of the asymptotic
 * expansion, the relative error stays below 1e-12 per glibc
 * documentation.
 */
#ifndef K26ASTRO_LASER_AIRY_H
#define K26ASTRO_LASER_AIRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Diffraction-limited + jitter half-angle (rad). */
double k26astro_airy_half_angle(double wavelength_m,
                                double primary_diam_m,
                                double m_squared,
                                double pointing_jitter_rad);

/* Spot diameter (m) at range R. Returns 2 R θ. */
double k26astro_airy_spot_diameter(double range_m,
                                   double wavelength_m,
                                   double primary_diam_m,
                                   double m_squared,
                                   double pointing_jitter_rad);

/* Encircled-energy fraction inside radius r at range R for a
 * diffraction-limited beam with primary D and wavelength λ.
 * Range R > 0, r > 0. */
double k26astro_airy_encircled_fraction(double radius_at_range_m,
                                        double range_m,
                                        double primary_diam_m,
                                        double wavelength_m);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_LASER_AIRY_H */
