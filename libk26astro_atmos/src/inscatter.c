/* inscatter.c — single-scatter Rayleigh + Henyey-Greenstein Mie.
 *
 * Computes the in-scattered radiance along a view ray, integrating
 * the per-altitude scattering coefficient against the phase
 * functions. Single-scatter only; LUT-based multi-scattering /
 * shadow integration (Hillaire 2020) is not currently provided.
 *
 * Phase functions:
 *   Rayleigh:  P_R(cos θ) = (3/4) · (1 + cos²θ)
 *   Mie (HG):  P_M(cos θ) = (1 - g²) / (4π (1 + g² - 2g·cos θ)^1.5)
 *
 * Extinction along ray uses the same exp(-h/H) density profile as
 * the refraction code. Sample 32 steps along the view ray clipped
 * to the atmospheric top (h = atmos_top_m).
 *
 * Output: K26V3 triple in arbitrary linear units representing
 * R/G/B-weighted radiance. The Rayleigh coefficient β_R is treated
 * as a single value at 550 nm green; full spectral decomposition
 * (Rayleigh ∝ 1/λ⁴) is not currently provided. The R/G/B output
 * here assumes a fixed blue-skewed weighting: R = 0.4 β,
 * G = 1.0 β, B = 2.5 β (rough approximation of the wavelength
 * dependence). */
#include "k26astro_atmos/atmos.h"
#include "k26astro_atmos/atmos_consts.h"

#include <math.h>

double k26astro_atmos_density_ratio_(const K26AstroAtmos *a, double h_m);

static double v3_norm_(K26V3 v)
{
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static double v3_dot_(K26V3 a, K26V3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static double phase_rayleigh_(double cos_theta)
{
    return 0.75 * (1.0 + cos_theta * cos_theta);
}

static double phase_mie_hg_(double cos_theta, double g)
{
    double denom = 1.0 + g*g - 2.0*g*cos_theta;
    if (!(denom > 0.0)) return 0.0;
    double inv = (1.0 - g*g) / (4.0 * M_PI * pow(denom, 1.5));
    return inv;
}

K26V3 k26astro_atmos_inscatter(const K26AstroAtmos *a,
                                K26V3 view_origin_m,
                                K26V3 view_dir,
                                K26V3 sun_dir)
{
    K26V3 zero = { 0.0, 0.0, 0.0 };
    if (!a) return zero;

    K26AstroAtmosParams p = k26astro_atmos_params(a);
    if (!(p.scale_height_m > 0.0)) return zero;

    /* Normalise direction inputs. */
    double vn = v3_norm_(view_dir);
    double sn = v3_norm_(sun_dir);
    if (!(vn > 0.0) || !(sn > 0.0)) return zero;
    K26V3 vd = { view_dir.x / vn, view_dir.y / vn, view_dir.z / vn };
    K26V3 sd = { sun_dir.x  / sn, sun_dir.y  / sn, sun_dir.z  / sn };

    double cos_theta = v3_dot_(vd, sd);
    double P_R = phase_rayleigh_(cos_theta);
    double P_M = phase_mie_hg_(cos_theta, p.mie_asymmetry_g);

    /* Integrate along view ray from origin to atmosphere exit.
     * The view_origin_m is in planet-centric coordinates; the ray's
     * altitude profile is parameterised by t (metres along ray).
     * h(t) = |view_origin + t · vd| - planet_radius. We don't have
     * an explicit planet radius here, so use the convention that
     * view_origin's magnitude is the planet centre distance and
     * "altitude" h = |position| - reference_R. For Earth this is
     * ~6371 km; for the v0.2 surface we approximate by extracting
     * the radial component along the up direction (zenith ≈
     * view_origin direction). */
    double r0 = v3_norm_(view_origin_m);
    if (!(r0 > 0.0)) return zero;
    K26V3 up = { view_origin_m.x / r0,
                  view_origin_m.y / r0,
                  view_origin_m.z / r0 };

    /* Effective surface radius: r0 minus observer altitude. We
     * don't have that; assume observer is at sea level so
     * R_surface = r0. Caller may set view_origin at an altitude
     * above ground, in which case the integration starts above
     * the surface. */
    double R_surf = r0;
    /* Find ray exit altitude. Solve |view_origin + t·vd|² =
     * (R_surf + atmos_top)² for the positive root. */
    double R_top  = R_surf + p.atmos_top_m;
    double b_dot  = v3_dot_(view_origin_m, vd);
    double disc   = b_dot * b_dot - (r0*r0 - R_top*R_top);
    double t_max;
    if (disc < 0.0) {
        /* Ray doesn't enter atmosphere top — shouldn't happen if
         * observer is inside; cap at scale_height for sanity. */
        t_max = p.scale_height_m;
    } else {
        double sqrt_disc = sqrt(disc);
        double t_plus  = -b_dot + sqrt_disc;
        t_max = (t_plus > 0.0) ? t_plus : p.scale_height_m;
    }
    if (t_max > 1.0e7) t_max = 1.0e7;  /* hard cap 10000 km */

    /* 32-step trapezoidal integration. */
    const int N = 32;
    double dt = t_max / (double)N;
    double accum_R = 0.0;
    double accum_M = 0.0;

    /* (void) up suppresses unused warning */
    (void)up;

    for (int i = 0; i <= N; i++) {
        double t = i * dt;
        K26V3 sample = {
            view_origin_m.x + t * vd.x,
            view_origin_m.y + t * vd.y,
            view_origin_m.z + t * vd.z
        };
        double r_sample = v3_norm_(sample);
        double h_sample = r_sample - R_surf;
        if (h_sample < 0.0) h_sample = 0.0;
        double rho_ratio = k26astro_atmos_density_ratio_(a, h_sample);
        double weight = (i == 0 || i == N) ? 0.5 : 1.0;
        accum_R += rho_ratio * weight;
        accum_M += rho_ratio * weight;
    }
    accum_R *= dt * p.rayleigh_coeff;
    accum_M *= dt * p.mie_coeff;

    double scattered_R = P_R * accum_R;
    double scattered_M = P_M * accum_M;

    /* Wavelength-dependent weighting (placeholder Rayleigh
     * scaling): R, G, B scales of (0.4, 1.0, 2.5) crudely
     * approximate the 1/λ⁴ behaviour. Mie is wavelength-flat
     * to first order. */
    K26V3 out = {
        0.4 * scattered_R + scattered_M,
        1.0 * scattered_R + scattered_M,
        2.5 * scattered_R + scattered_M
    };
    return out;
}
