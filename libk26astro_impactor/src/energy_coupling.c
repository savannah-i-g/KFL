/* energy_coupling.c — impactor impact → energy delivered to target. */
#include "k26astro_impactor/energy_coupling.h"

#include <stddef.h>

static double resolve_eff_(double user_value, double fallback)
{
    return (user_value > 0.0) ? user_value : fallback;
}

double k26astro_impactor_energy_delivered_j(
    const K26AstroImpactor *projectile,
    const K26AstroImpactEvent *impact,
    const K26AstroTargetStructureSpec *target)
{
    if (!projectile || !impact || !target) return 0.0;

    const double m = k26astro_impactor_projectile_mass(projectile);
    const double v = impact->impact_velocity_m_per_s;
    if (m <= 0.0 || v <= 0.0) return 0.0;
    const double impact_KE = 0.5 * m * v * v;

    /* Resolve coupling efficiencies (use defaults when caller
     * passed <= 0). */
    const double eff_surface =
        resolve_eff_(target->surface_spallation_eff, 0.10);
    const double eff_bumper_partial =
        resolve_eff_(target->bumper_partial_eff, 0.50);
    const double eff_full =
        resolve_eff_(target->full_penetration_eff, 0.90);

    /* Whipple branch. */
    const int has_whipple = (target->outer_thickness_m > 0.0)
                         && (target->bumper_gap_m > 0.0);
    if (has_whipple) {
        if (!impact->penetrates) {
            /* Bumper stopped the projectile (or its critical
             * diameter exceeded the projectile's). Surface
             * spallation dissipates most of the impact KE in the
             * outer skin; only a small fraction couples inward. */
            return impact_KE * eff_surface;
        }
        if (target->inner_thickness_m > 0.0) {
            /* Bumper breached, inner wall present. The Whipple
             * stand-off + inner wall together absorb the spray
             * cone but allow medium energy coupling to the
             * interior. */
            return impact_KE * eff_bumper_partial;
        }
        /* Bumper-only target — projectile reaches interior. */
        return impact_KE * eff_full;
    }

    /* Monolithic branch. */
    if (target->monolithic_thickness_m > 0.0) {
        const double t_wall = target->monolithic_thickness_m;
        const double p_depth = impact->monolithic_penetration_m;
        if (p_depth >= t_wall) {
            /* Full penetration. */
            return impact_KE * eff_full;
        }
        /* Partial penetration: quadratic falloff with depth ratio
         * (captures the crater-volume scaling on hypervelocity
         * impact). */
        const double r = p_depth / t_wall;
        return impact_KE * (r * r) * eff_full;
    }

    /* No target geometry described — return worst-case energy
     * coupling. Consumers should provide a structure spec; this
     * fallback exists so the call doesn't return 0 silently for
     * a caller who forgot to populate target. */
    return impact_KE * eff_full;
}
