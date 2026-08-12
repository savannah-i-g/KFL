/* k26astro_impactor/impactor.h — kinetic impactor.
 *
 * libk26astro_impactor composes against libk26astro_vehicle's
 * PAYLOAD slot via the shared K26AstroPayload base. A single
 * impactor carries the projectile parameters (mass, density,
 * diameter) plus the release-pattern selection (single body vs
 * multi-unit swarm). The mission driver dispatches the impactor
 * by creating world bodies from the lib's dispatch helpers; world
 * propagation + intercept geometry then runs through
 * libk26astro_grav.
 *
 * Stated limitations:
 *   - The ballistic-limit equation is the Christiansen high-
 *     velocity regime (V >= 7 km/s); subsonic / hypersonic-
 *     transition regimes are not modelled.
 *   - Modified Cour-Palais monolithic penetration assumes a
 *     homogeneous target; layered / composite targets need
 *     per-layer calls.
 *   - Plume / wake dispersion of the projectile during cruise
 *     is not modelled.
 *
 * References:
 *   Christiansen 1993     IJIE 14:145-156, Eq. 11 ballistic limit
 *   Christiansen 2009     NASA-SP-2009-573, MMOD protection
 *   Hayashida & Robinson 1991   NASA TM-103565, single-wall penetration
 */
#ifndef K26ASTRO_IMPACTOR_IMPACTOR_H
#define K26ASTRO_IMPACTOR_IMPACTOR_H

#include "k26astro_impactor/impactor_consts.h"
#include "k26astro_impactor/whipple.h"
#include "k26astro_impactor/swarm.h"
#include "k26astro_impactor/terminal_phase.h"

#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    K26ASTRO_IMPACTOR_PATTERN_SINGLE    = 1,
    K26ASTRO_IMPACTOR_PATTERN_SWARM = 2
} K26AstroImpactorPattern;

typedef struct K26AstroImpactor K26AstroImpactor;

/* Construct an impactor. */
K26AstroImpactor *k26astro_impactor_new(
    K26AstroImpactorPattern pattern,
    double             projectile_mass_kg,
    double             projectile_density_kg_per_m3,
    double             projectile_diameter_m,
    int                swarm_count,       /* ignored for SINGLE */
    double             swarm_half_angle_rad);

void k26astro_impactor_destroy(K26AstroImpactor *k);
int  k26astro_impactor_bind(K26AstroImpactor *k, struct K26AstroVehicle *v);
struct K26AstroPayload *k26astro_impactor_as_payload(K26AstroImpactor *k);

K26AstroImpactorPattern k26astro_impactor_pattern(const K26AstroImpactor *k);
double k26astro_impactor_projectile_mass(const K26AstroImpactor *k);
double k26astro_impactor_projectile_density(const K26AstroImpactor *k);
double k26astro_impactor_projectile_diameter(const K26AstroImpactor *k);
int    k26astro_impactor_swarm_count(const K26AstroImpactor *k);
double k26astro_impactor_swarm_half_angle(const K26AstroImpactor *k);

/* ---- Impact event -------------------------------------------- */

typedef struct {
    /* Inputs */
    double impact_velocity_m_per_s;
    double cos_impact_angle;
    /* Whipple analysis */
    double critical_diameter_m;
    int    penetrates;
    /* Monolithic-plate fallback (relevant when Whipple parameters
     * are zero / target is monolithic). Uses Hayashida-Robinson
     * §2.5 Modified Cour-Palais form. */
    double monolithic_penetration_m;
} K26AstroImpactEvent;

/* Analyse one impact.
 *
 * Whipple analysis (Christiansen 1993 Eq. 11) runs when
 * target_wall_thickness_m > 0; otherwise skipped. Monolithic
 * analysis (Hayashida-Robinson 1991 §2.5) runs when
 * target_brinell_hardness > 0 and target_density_kg_per_m3 > 0;
 * otherwise skipped. The two analyses are independent and both
 * fields are populated when both sets of parameters are valid. */
K26AstroImpactEvent k26astro_impactor_analyse_impact(
    const K26AstroImpactor *projectile,
    double             impact_velocity_m_per_s,
    double             cos_impact_angle,
    /* Whipple shield parameters: */
    double             target_wall_thickness_m,
    double             target_bumper_density_kg_per_m3,
    double             target_bumper_spacing_m,
    double             target_wall_yield_stress_ksi,
    /* Monolithic-plate parameters: */
    double             target_brinell_hardness,
    double             target_density_kg_per_m3,
    double             target_speed_of_sound_m_per_s);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_IMPACTOR_IMPACTOR_H */
