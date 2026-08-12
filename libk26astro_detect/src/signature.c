/* signature.c — IR and RCS signature math.
 *
 * Implements:
 *   - Bolometric Stefan-Boltzmann emission
 *   - In-band Planck integration via Simpson's rule
 *   - Lambertian projected intensity
 *   - Faceted geometric-optics RCS
 *
 * All functions are pure (no global state). The lib's compilation
 * flags (-ffp-contract=off, -fexcess-precision=standard) bound the
 * floating-point reproducibility across reference CPUs. */

#include "k26astro_detect/signature.h"
#include "k26astro_detect/detect_consts.h"

#include "k26astro_core/consts.h"
#include "k26m3d.h"

#include <math.h>
#include <stddef.h>

/* ---- Bolometric ----------------------------------------------- */

double k26astro_signature_ir_bolometric_power(double area_m2,
                                              double epsilon,
                                              double T_K)
{
    if (area_m2 <= 0.0 || epsilon <= 0.0 || T_K <= 0.0) return 0.0;
    if (epsilon > 1.0) epsilon = 1.0;
    return epsilon * K26ASTRO_DETECT_SIGMA_SB * area_m2 *
           T_K * T_K * T_K * T_K;
}

/* ---- Planck-integrated in-band power -------------------------- */

/* Planck spectral radiance per wavelength (W · m^-2 · sr^-1 · m^-1):
 *     B_λ(T) = (2 h c^2 / λ^5) / (exp(h c / (λ k T)) - 1) */
static double planck_radiance_(double lambda_m, double T_K)
{
    if (lambda_m <= 0.0 || T_K <= 0.0) return 0.0;
    const double h = K26A_H_PLANCK;
    const double c = K26A_C;
    const double k = K26A_K_BOLTZMANN;
    const double lambda5 = lambda_m * lambda_m * lambda_m *
                           lambda_m * lambda_m;
    const double numer = 2.0 * h * c * c / lambda5;
    const double arg   = h * c / (lambda_m * k * T_K);
    /* For very large arg the denominator overflows; the integrand
     * is negligibly small there. Guard against overflow with an
     * exponential cap; physically realistic IR-band integrands at
     * ~300 K with λ in the μm range have arg ~ 50, well within
     * double-precision exp range. */
    if (arg > 700.0) return 0.0;
    return numer / (exp(arg) - 1.0);
}

double k26astro_signature_ir_planck_inband(double area_m2,
                                            double epsilon,
                                            double T_K,
                                            double lambda_lo_m,
                                            double lambda_hi_m,
                                            int n_quad)
{
    if (area_m2 <= 0.0 || epsilon <= 0.0 || T_K <= 0.0) return 0.0;
    if (lambda_hi_m <= lambda_lo_m) return 0.0;
    if (n_quad < 2) n_quad = 2;
    /* Force even n_quad for Simpson's rule. */
    if (n_quad % 2) n_quad++;
    if (epsilon > 1.0) epsilon = 1.0;

    /* Composite Simpson's rule on B_λ over [lambda_lo, lambda_hi]
     * with log-λ change of variable: ∫ B_λ dλ = ∫ B_λ · λ d(log λ).
     * Log-spacing concentrates samples near the Wien peak (which
     * is at λ_peak ≈ 2.898e-3 / T_K) and lets a wide pass band
     * span many decades without undersampling the peak region. */
    const double log_lo = log(lambda_lo_m);
    const double log_hi = log(lambda_hi_m);
    const double h = (log_hi - log_lo) / (double)n_quad;

    double f_lo = planck_radiance_(lambda_lo_m, T_K) * lambda_lo_m;
    double f_hi = planck_radiance_(lambda_hi_m, T_K) * lambda_hi_m;
    double sum = f_lo + f_hi;
    for (int i = 1; i < n_quad; i++) {
        const double lambda = exp(log_lo + h * (double)i);
        const double b = planck_radiance_(lambda, T_K) * lambda;
        sum += ((i & 1) ? 4.0 : 2.0) * b;
    }
    const double radiance = (h / 3.0) * sum;
    /* π comes from Lambertian integration; ε is the greybody
     * emissivity; area is the projected geometric area. */
    return epsilon * area_m2 * K26A_PI * radiance;
}

/* ---- Lambertian projected intensity --------------------------- */

double k26astro_signature_ir_lambertian(double area_m2,
                                        double epsilon,
                                        double T_K,
                                        double cos_view_angle)
{
    if (cos_view_angle <= 0.0) return 0.0;
    if (cos_view_angle > 1.0) cos_view_angle = 1.0;
    const double P = k26astro_signature_ir_bolometric_power(
        area_m2, epsilon, T_K);
    /* I = (P / π) cos θ — Lambertian phase function over the
     * hemisphere. */
    return (P / K26A_PI) * cos_view_angle;
}

/* ---- Faceted RCS ---------------------------------------------- */

double k26astro_signature_rcs_facet_geometry(int n_facets,
                                              const K26V3 *normals_body,
                                              const double *areas_m2,
                                              K26V3 look_unit_body)
{
    if (n_facets <= 0 || !normals_body || !areas_m2) return 0.0;

    /* The geometric sum Σ A_i^2 cos^4 θ_i for facets with
     * cos θ_i = dot(normal_i, -look) > 0. Sign convention:
     * look_unit_body points from observer toward target, so the
     * outward-facing normals that scatter toward the observer
     * satisfy dot(n, -look) > 0. */
    const K26V3 minus_look = k26m3d_v3_neg(look_unit_body);
    double accum = 0.0;
    for (int i = 0; i < n_facets; i++) {
        const double ndot = k26m3d_v3_dot(normals_body[i], minus_look);
        if (ndot <= 0.0) continue;
        const double A = areas_m2[i];
        const double cos4 = ndot * ndot * ndot * ndot;
        accum += A * A * cos4;
    }
    return accum;
}

double k26astro_signature_rcs_monostatic(int n_facets,
                                          const K26V3 *normals_body,
                                          const double *areas_m2,
                                          K26V3 look_unit_body,
                                          double wavelength_m)
{
    if (wavelength_m <= 0.0) return 0.0;
    const double geom = k26astro_signature_rcs_facet_geometry(
        n_facets, normals_body, areas_m2, look_unit_body);
    /* Knott 2004 eq. 5.45 prefactor: 4π / λ^2. */
    return (4.0 * K26A_PI / (wavelength_m * wavelength_m)) * geom;
}
