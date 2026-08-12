/* whipple.h — Whipple-shield + monolithic-plate impact scaling.
 *
 * Christiansen 1993 (IJIE 14:145) Eq. 11 ballistic-limit equation
 * for Whipple shields, NNO (no-perforation) regime V_n ≥ 7 km/s.
 * Native form is cgs/ksi; the SI-equivalent equation absorbs the
 * unit conversion into K_CHRISTIANSEN_SI ≈ 8441:
 *
 *     d_c[m] = K_SI · t_w^(2/3) · ρ_b^(-1/9) · ρ_p^(-1/3) ·
 *              (V · cos θ)^(-2/3) · S^(1/3) · (σ_w / 70)^(1/3)
 *
 * where:
 *   d_c   critical projectile diameter (just-perforates), m
 *   t_w   rear-wall thickness, m
 *   ρ_b   bumper density, kg/m³
 *   ρ_p   projectile density, kg/m³
 *   V     impact velocity, m/s
 *   θ     impact angle from normal (cos θ ∈ (0, 1])
 *   S     bumper-to-wall spacing, m
 *   σ_w   rear-wall yield stress, ksi (Christiansen pivots at
 *         70 ksi for Al 6061-T6; the (σ/70)^(1/3) factor stays
 *         dimensionless)
 *
 * Returning d_c lets the caller compare against the actual
 * projectile diameter: if d_proj > d_c, the projectile
 * perforates the rear wall.
 *
 * Monolithic-plate scaling for unshielded targets is exposed via
 * the separate `_monolithic_penetration_depth` API in this header
 * (Hayashida-Robinson 1991 §2.5 Modified Cour-Palais).
 *
 * References:
 *   Christiansen, E.L. 1993. Design and performance equations for
 *     advanced meteoroid and debris shields. Int. J. Impact
 *     Engng. 14(1-4):145-156. Eq. 11.
 *   Hayashida, K.B., Robinson, J.H. 1991. Single wall penetration
 *     equations. NASA TM-103565, §2.5.
 */
#ifndef K26ASTRO_IMPACTOR_WHIPPLE_H
#define K26ASTRO_IMPACTOR_WHIPPLE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Critical-just-perforates projectile diameter (m) for a Whipple
 * shield, Christiansen 1993 Eq. 11. Caller compares against the
 * actual projectile diameter: d_proj > d_c means perforation. */
double k26astro_whipple_critical_diameter(double wall_thickness_m,
                                           double bumper_density_kg_per_m3,
                                           double projectile_density_kg_per_m3,
                                           double impact_velocity_m_per_s,
                                           double cos_impact_angle,
                                           double bumper_spacing_m,
                                           double wall_yield_stress_ksi);

/* Boolean: does the projectile perforate the Whipple shield?
 * Returns 1 if yes (d_proj > d_c), 0 otherwise. */
int k26astro_whipple_penetrates(double wall_thickness_m,
                                 double bumper_density_kg_per_m3,
                                 double projectile_diameter_m,
                                 double projectile_density_kg_per_m3,
                                 double impact_velocity_m_per_s,
                                 double cos_impact_angle,
                                 double bumper_spacing_m,
                                 double wall_yield_stress_ksi);

/* Monolithic-plate penetration depth (m). Modified Cour-Palais
 * form, Hayashida & Robinson 1991 §2.5. */
double k26astro_monolithic_penetration_depth(double projectile_diameter_m,
                                              double target_brinell_hardness,
                                              double projectile_density_kg_per_m3,
                                              double target_density_kg_per_m3,
                                              double impact_velocity_m_per_s,
                                              double cos_impact_angle,
                                              double target_speed_of_sound_m_per_s);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_IMPACTOR_WHIPPLE_H */
