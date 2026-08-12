/* energy_coupling.h — impact → delivered-energy coupling.
 *
 * The existing impact analysis at k26astro_impactor_analyse_impact
 * returns penetration metrics (critical Whipple diameter, monolithic
 * penetration depth) but does not characterise how much of the
 * projectile's kinetic energy is actually delivered to the target's
 * interior structures. This header adds that final step.
 *
 * The delivered-energy fraction depends on whether the projectile
 * is stopped by the outer bumper (low fraction — most energy is
 * absorbed by outer-skin spallation), breaches the bumper but
 * doesn't reach the inner wall (medium fraction — energy partially
 * deflected and spread by the bumper), or fully penetrates
 * (high fraction — bulk of energy couples to interior structures).
 *
 * For a monolithic target the penetration depth vs. wall thickness
 * ratio drives the coupling: full penetration delivers ~90%,
 * partial penetration scales quadratically with depth ratio.
 *
 * Rationale: a driver that treats delivered energy as raw ½mv²
 * regardless of target structure collapses shielded and unshielded
 * targets into the same outcome and yields discrete per-class
 * damage steps. With this API the delivered energy scales with
 * the penetration state of the target. */
#ifndef K26ASTRO_IMPACTOR_ENERGY_COUPLING_H
#define K26ASTRO_IMPACTOR_ENERGY_COUPLING_H

#include "k26astro_impactor/impactor.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Target structure specification for energy-coupling computation.
 *
 * A target can have Whipple-shielded geometry (outer bumper +
 * stand-off + inner wall), monolithic-plate geometry (single
 * thick wall), or both (the analysis runs both paths and uses the
 * worst-case for the projectile).
 *
 * Coupling efficiencies are the fraction of impact KE delivered
 * to the target interior in each regime; defaults are applied
 * when the caller passes zero. */
typedef struct {
    /* Whipple bumper layer. Zero outer_thickness disables the
     * Whipple analysis branch. */
    double outer_thickness_m;
    double outer_yield_stress_ksi;
    double bumper_density_kg_per_m3;
    double bumper_gap_m;

    /* Inner structural-wall layer for the Whipple
     * stand-off geometry. Optional — zero treats the Whipple
     * geometry as bumper-only. */
    double inner_thickness_m;
    double inner_yield_stress_ksi;

    /* Monolithic-plate fallback parameters. Zero density / hardness
     * disables the monolithic analysis branch (use Whipple only). */
    double monolithic_thickness_m;
    double monolithic_brinell_hardness;
    double monolithic_density_kg_per_m3;
    double monolithic_speed_of_sound_m_per_s;

    /* Coupling efficiencies. If a value is <= 0 the default below
     * is applied:
     *   surface_spallation_eff   default 0.10
     *   bumper_partial_eff       default 0.50
     *   full_penetration_eff     default 0.90 */
    double surface_spallation_eff;
    double bumper_partial_eff;
    double full_penetration_eff;
} K26AstroTargetStructureSpec;

/* Compute the energy delivered to the target's interior in joules.
 *
 * Decision logic:
 *
 *   if Whipple geometry (outer_thickness > 0 && bumper_gap > 0):
 *     - impact analysis already gave us impact.penetrates
 *       (which means the projectile diameter exceeds the Whipple
 *       critical diameter d_c, so the bumper is breached).
 *     - if !penetrates:           E = impact_KE × surface_spallation_eff
 *     - elif inner_thickness > 0: E = impact_KE × bumper_partial_eff
 *       (bumper breached but inner wall stops the spray cone — a
 *       full inner-wall penetration model is a future extension)
 *     - else:                     E = impact_KE × full_penetration_eff
 *       (bumper-only target, projectile reaches interior after
 *       bumper breach)
 *
 *   elif Monolithic geometry (monolithic_thickness > 0):
 *     - if impact.monolithic_penetration_m >= monolithic_thickness:
 *           E = impact_KE × full_penetration_eff
 *     - else:
 *           E = impact_KE × (depth_ratio² × full_penetration_eff)
 *       (quadratic falloff for partial penetration — captures the
 *       crater-volume scaling on hypervelocity impact)
 *
 *   else (no target geometry described):
 *     - E = impact_KE × full_penetration_eff (worst-case full
 *       coupling; consumers should provide target geometry).
 *
 * The returned energy is bounded by 0 <= E <= impact_KE. NULL
 * inputs return 0. */
double k26astro_impactor_energy_delivered_j(
    const K26AstroImpactor *projectile,
    const K26AstroImpactEvent *impact,
    const K26AstroTargetStructureSpec *target);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_IMPACTOR_ENERGY_COUPLING_H */
