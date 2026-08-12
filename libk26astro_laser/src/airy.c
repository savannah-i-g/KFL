/* airy.c — Airy diffraction-limited spot geometry. */

#include "k26astro_laser/airy.h"
#include "k26astro_laser/laser_consts.h"

#include "k26astro_core/consts.h"

#include <math.h>

double k26astro_airy_half_angle(double wavelength_m,
                                double primary_diam_m,
                                double m_squared,
                                double pointing_jitter_rad)
{
    if (wavelength_m <= 0.0 || primary_diam_m <= 0.0) return 0.0;
    if (m_squared < 1.0) m_squared = 1.0;
    if (pointing_jitter_rad < 0.0) pointing_jitter_rad = 0.0;

    const double theta_diff =
        K26ASTRO_LASER_AIRY_FACTOR * m_squared *
        wavelength_m / primary_diam_m;
    /* RSS jitter into the half-angle. */
    return sqrt(theta_diff * theta_diff +
                pointing_jitter_rad * pointing_jitter_rad);
}

double k26astro_airy_spot_diameter(double range_m,
                                   double wavelength_m,
                                   double primary_diam_m,
                                   double m_squared,
                                   double pointing_jitter_rad)
{
    if (range_m <= 0.0) return 0.0;
    const double theta = k26astro_airy_half_angle(
        wavelength_m, primary_diam_m, m_squared, pointing_jitter_rad);
    return 2.0 * range_m * theta;
}

/* Encircled-energy fraction F(x) = 1 - J_0(x)^2 - J_1(x)^2 with
 * x = π D r / (λ R). At x = 3.8317 (first dark ring),
 * F ≈ 0.8378 — the canonical Airy disc value. */
double k26astro_airy_encircled_fraction(double radius_at_range_m,
                                        double range_m,
                                        double primary_diam_m,
                                        double wavelength_m)
{
    if (radius_at_range_m <= 0.0 || range_m <= 0.0) return 0.0;
    if (primary_diam_m <= 0.0 || wavelength_m <= 0.0) return 0.0;

    const double x = K26A_PI * primary_diam_m * radius_at_range_m /
                     (wavelength_m * range_m);
    const double j0v = j0(x);
    const double j1v = j1(x);
    double F = 1.0 - j0v * j0v - j1v * j1v;
    if (F < 0.0) F = 0.0;
    if (F > 1.0) F = 1.0;
    return F;
}
