/* whipple.c — Christiansen 1993 Eq. 11 Whipple ballistic limit
 * + Hayashida-Robinson §2.5 Modified Cour-Palais monolithic
 * penetration. */

#include "k26astro_impactor/whipple.h"
#include "k26astro_impactor/impactor_consts.h"

#include <math.h>

double k26astro_whipple_critical_diameter(double wall_thickness_m,
                                           double bumper_density_kg_per_m3,
                                           double projectile_density_kg_per_m3,
                                           double impact_velocity_m_per_s,
                                           double cos_impact_angle,
                                           double bumper_spacing_m,
                                           double wall_yield_stress_ksi)
{
    if (wall_thickness_m <= 0.0) return 0.0;
    if (bumper_density_kg_per_m3 <= 0.0) return 0.0;
    if (projectile_density_kg_per_m3 <= 0.0) return 0.0;
    if (impact_velocity_m_per_s <= 0.0) return 0.0;
    if (cos_impact_angle <= 0.0) return 0.0;
    if (cos_impact_angle > 1.0) cos_impact_angle = 1.0;
    if (bumper_spacing_m <= 0.0) return 0.0;
    if (wall_yield_stress_ksi <= 0.0) return 0.0;

    /* Christiansen 1993 IJIE 14:145, Eq. 11 (V_n ≥ 7 km/sec):
     *   d_c[m] = K_SI · t_w^(2/3) · ρ_b^(-1/9) · ρ_p^(-1/3) ·
     *            (V · cos θ)^(-2/3) · S^(1/3) · (σ_w / 70)^(1/3)
     */
    const double V_n = impact_velocity_m_per_s * cos_impact_angle;
    return K26ASTRO_IMPACTOR_K_CHRISTIANSEN_SI
         * pow(wall_thickness_m,         2.0 / 3.0)
         * pow(bumper_density_kg_per_m3,    -1.0 / 9.0)
         * pow(projectile_density_kg_per_m3, -1.0 / 3.0)
         * pow(V_n,                         -2.0 / 3.0)
         * pow(bumper_spacing_m,             1.0 / 3.0)
         * pow(wall_yield_stress_ksi / 70.0, 1.0 / 3.0);
}

int k26astro_whipple_penetrates(double wall_thickness_m,
                                 double bumper_density_kg_per_m3,
                                 double projectile_diameter_m,
                                 double projectile_density_kg_per_m3,
                                 double impact_velocity_m_per_s,
                                 double cos_impact_angle,
                                 double bumper_spacing_m,
                                 double wall_yield_stress_ksi)
{
    if (projectile_diameter_m <= 0.0) return 0;
    const double d_c = k26astro_whipple_critical_diameter(
        wall_thickness_m, bumper_density_kg_per_m3,
        projectile_density_kg_per_m3, impact_velocity_m_per_s,
        cos_impact_angle, bumper_spacing_m, wall_yield_stress_ksi);
    if (d_c <= 0.0) return 0;
    return projectile_diameter_m > d_c ? 1 : 0;
}

double k26astro_monolithic_penetration_depth(double projectile_diameter_m,
                                              double target_brinell_hardness,
                                              double projectile_density_kg_per_m3,
                                              double target_density_kg_per_m3,
                                              double impact_velocity_m_per_s,
                                              double cos_impact_angle,
                                              double target_speed_of_sound_m_per_s)
{
    if (projectile_diameter_m <= 0.0) return 0.0;
    if (target_brinell_hardness <= 0.0) return 0.0;
    if (projectile_density_kg_per_m3 <= 0.0) return 0.0;
    if (target_density_kg_per_m3 <= 0.0) return 0.0;
    if (impact_velocity_m_per_s <= 0.0) return 0.0;
    if (cos_impact_angle <= 0.0) return 0.0;
    if (cos_impact_angle > 1.0) cos_impact_angle = 1.0;
    if (target_speed_of_sound_m_per_s <= 0.0) return 0.0;

    /* Hayashida & Robinson 1991 §2.5 Modified Cour-Palais
     * (NASA TM-103565, page 4):
     *   p[cm] = 5.24 · d[cm]^(19/18) · BH^(-0.25) · (ρ_p/ρ_t)^0.5 ·
     *           (V_n / C)^(2/3)
     * Units: d, p in cm; V_n, C in km/s; BH dimensionless. */
    const double V_n_km_s = impact_velocity_m_per_s * cos_impact_angle * 1.0e-3;
    const double C_km_s   = target_speed_of_sound_m_per_s * 1.0e-3;
    const double d_cm     = projectile_diameter_m * 1.0e2;
    const double rho_ratio = projectile_density_kg_per_m3 /
                             target_density_kg_per_m3;
    return 5.24
         * pow(d_cm, 19.0 / 18.0)
         * pow(target_brinell_hardness, -0.25)
         * pow(rho_ratio, 0.5)
         * pow(V_n_km_s / C_km_s, 2.0 / 3.0)
         * 1.0e-2;  /* cm → m */
}
