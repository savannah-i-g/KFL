/* atmos.c — K26AstroAtmos lifecycle + density profile.
 *
 * Holds the K26AstroAtmosParams supplied by the caller (or by the
 * Earth-standard preset) and exposes the barometric density
 * function used by both refraction and inscatter. */
#include "k26astro_atmos/atmos.h"
#include "k26astro_atmos/atmos_consts.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct K26AstroAtmos {
    K26AstroAtmosParams p;
    /* Derived: rho_0 = P_0 / (R_specific · T_0). For Earth standard
     * we use the ICAO 1993 value ρ_0 = 1.225 kg/m³ directly. For
     * non-Earth atmospheres the caller should set ground_pressure_pa
     * and we approximate ρ_0 from pressure scaling. */
    double rho_0_kg_m3;
};

K26AstroAtmos *k26astro_atmos_init(const K26AstroAtmosParams *p)
{
    if (!p) return NULL;
    if (!(p->scale_height_m > 0.0))     return NULL;
    if (!(p->rayleigh_coeff >= 0.0))    return NULL;
    if (!(p->mie_coeff >= 0.0))         return NULL;
    if (!(p->ground_pressure_pa >= 0.0)) return NULL;
    if (!(p->atmos_top_m > 0.0))        return NULL;
    if (p->mie_asymmetry_g <= -1.0 || p->mie_asymmetry_g >= 1.0)
        return NULL;

    K26AstroAtmos *a = (K26AstroAtmos *)calloc(1, sizeof(*a));
    if (!a) return NULL;
    a->p = *p;
    /* ρ_0 from ICAO standard if Earth-like ground pressure, else
     * approximate as P_0/(R_air · T_0) with R_air = 287.05 J/(kg·K)
     * and T_0 = 288.15 K. */
    if (p->ground_pressure_pa > 0.0) {
        a->rho_0_kg_m3 = p->ground_pressure_pa / (287.05 * 288.15);
    } else {
        a->rho_0_kg_m3 = 1.225;
    }
    return a;
}

void k26astro_atmos_destroy(K26AstroAtmos *a)
{
    if (!a) return;
    free(a);
}

K26AstroAtmos *k26astro_atmos_earth_standard(void)
{
    K26AstroAtmosParams p = {
        .scale_height_m              = K26ASTRO_ATMOS_EARTH_H_M,
        .rayleigh_coeff              = K26ASTRO_ATMOS_EARTH_BETA_R,
        .mie_coeff                   = K26ASTRO_ATMOS_EARTH_BETA_M,
        .ground_pressure_pa          = K26ASTRO_ATMOS_P0_PA,
        .mie_asymmetry_g             = K26ASTRO_ATMOS_EARTH_HG_G,
        .n0                          = K26ASTRO_ATMOS_EARTH_N0,
        .atmos_top_m                 = K26ASTRO_ATMOS_EARTH_TOP_M,
        .t_0_k                       = K26ASTRO_ATMOS_EARTH_T0_K,
        .gamma_ratio_specific_heats  = K26ASTRO_ATMOS_EARTH_GAMMA,
        .r_specific_j_per_kg_k       = K26ASTRO_ATMOS_EARTH_R_SPECIFIC
    };
    return k26astro_atmos_init(&p);
}

double k26astro_atmos_density_at(const K26AstroAtmos *a, double h_m)
{
    if (!a) return 0.0;
    /* Below sea level — treat as sea-level density (caller may be
     * inside a planet's surface model; don't synthesise growing
     * density for negative altitude). */
    if (h_m < 0.0) return a->rho_0_kg_m3;
    if (h_m > a->p.atmos_top_m) return 0.0;
    return a->rho_0_kg_m3 * exp(-h_m / a->p.scale_height_m);
}

K26AstroAtmosParams k26astro_atmos_params(const K26AstroAtmos *a)
{
    if (!a) {
        K26AstroAtmosParams zero;
        memset(&zero, 0, sizeof(zero));
        return zero;
    }
    return a->p;
}

/* Internal: refractive index excess (n - 1) at altitude h. Scales
 * with density: (n-1)(h) = n0 · ρ(h)/ρ_0 = n0 · exp(-h/H). */
double k26astro_atmos_n_excess_at_(const K26AstroAtmos *a, double h_m);
double k26astro_atmos_n_excess_at_(const K26AstroAtmos *a, double h_m)
{
    if (!a) return 0.0;
    if (h_m < 0.0)             return a->p.n0;
    if (h_m > a->p.atmos_top_m) return 0.0;
    return a->p.n0 * exp(-h_m / a->p.scale_height_m);
}

double k26astro_atmos_density_ratio_(const K26AstroAtmos *a, double h_m);
double k26astro_atmos_density_ratio_(const K26AstroAtmos *a, double h_m)
{
    if (!a) return 0.0;
    if (h_m < 0.0)             return 1.0;
    if (h_m > a->p.atmos_top_m) return 0.0;
    return exp(-h_m / a->p.scale_height_m);
}
