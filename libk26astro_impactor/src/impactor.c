/* impactor.c — kinetic-impactor opaque + impact analysis dispatcher. */

#include "k26astro_impactor/impactor.h"

#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#include <stdlib.h>
#include <string.h>

struct K26AstroImpactor {
    K26AstroPayload       base;
    K26AstroImpactorPattern   pattern;
    double               projectile_mass_kg;
    double               projectile_density_kg_per_m3;
    double               projectile_diameter_m;
    int                  swarm_count;
    double               swarm_half_angle_rad;
};

K26AstroImpactor *k26astro_impactor_new(
    K26AstroImpactorPattern pattern,
    double             projectile_mass_kg,
    double             projectile_density_kg_per_m3,
    double             projectile_diameter_m,
    int                swarm_count,
    double             swarm_half_angle_rad)
{
    if (pattern != K26ASTRO_IMPACTOR_PATTERN_SINGLE &&
        pattern != K26ASTRO_IMPACTOR_PATTERN_SWARM) return NULL;
    if (projectile_mass_kg <= 0.0) return NULL;
    if (projectile_density_kg_per_m3 <= 0.0) return NULL;
    if (projectile_diameter_m <= 0.0) return NULL;
    if (pattern == K26ASTRO_IMPACTOR_PATTERN_SWARM) {
        if (swarm_count <= 0) return NULL;
        if (swarm_half_angle_rad <= 0.0) return NULL;
    }

    K26AstroImpactor *k = (K26AstroImpactor *)calloc(1, sizeof(*k));
    if (!k) return NULL;
    k->base.kind             = K26ASTRO_DEFENSE_KIND_IMPACTOR;
    k->base.owner            = NULL;
    k->base.owner_slot       = -1;
    k->base.generation       = 1;
    k->base.on_owner_destroy = NULL;
    k->pattern               = pattern;
    k->projectile_mass_kg    = projectile_mass_kg;
    k->projectile_density_kg_per_m3 = projectile_density_kg_per_m3;
    k->projectile_diameter_m = projectile_diameter_m;
    k->swarm_count       = (pattern == K26ASTRO_IMPACTOR_PATTERN_SWARM)
                                ? swarm_count : 1;
    k->swarm_half_angle_rad = (pattern == K26ASTRO_IMPACTOR_PATTERN_SWARM)
                                ? swarm_half_angle_rad : 0.0;
    return k;
}

void k26astro_impactor_destroy(K26AstroImpactor *k)
{
    if (!k) return;
    k26astro_payload_unlink_(&k->base);
    free(k);
}

int k26astro_impactor_bind(K26AstroImpactor *k, struct K26AstroVehicle *v)
{
    if (!k) return -1;
    return k26astro_payload_bind(&k->base, v);
}

struct K26AstroPayload *k26astro_impactor_as_payload(K26AstroImpactor *k)
{
    return k ? &k->base : NULL;
}

K26AstroImpactorPattern k26astro_impactor_pattern(const K26AstroImpactor *k)
{
    return k ? k->pattern : K26ASTRO_IMPACTOR_PATTERN_SINGLE;
}

double k26astro_impactor_projectile_mass(const K26AstroImpactor *k)
{
    return k ? k->projectile_mass_kg : 0.0;
}

double k26astro_impactor_projectile_density(const K26AstroImpactor *k)
{
    return k ? k->projectile_density_kg_per_m3 : 0.0;
}

double k26astro_impactor_projectile_diameter(const K26AstroImpactor *k)
{
    return k ? k->projectile_diameter_m : 0.0;
}

int k26astro_impactor_swarm_count(const K26AstroImpactor *k)
{
    return k ? k->swarm_count : 0;
}

double k26astro_impactor_swarm_half_angle(const K26AstroImpactor *k)
{
    return k ? k->swarm_half_angle_rad : 0.0;
}

K26AstroImpactEvent k26astro_impactor_analyse_impact(
    const K26AstroImpactor *projectile,
    double             impact_velocity_m_per_s,
    double             cos_impact_angle,
    double             target_wall_thickness_m,
    double             target_bumper_density_kg_per_m3,
    double             target_bumper_spacing_m,
    double             target_wall_yield_stress_ksi,
    double             target_brinell_hardness,
    double             target_density_kg_per_m3,
    double             target_speed_of_sound_m_per_s)
{
    K26AstroImpactEvent ev;
    memset(&ev, 0, sizeof(ev));
    if (!projectile) return ev;

    ev.impact_velocity_m_per_s = impact_velocity_m_per_s;
    ev.cos_impact_angle        = cos_impact_angle;

    if (target_wall_thickness_m > 0.0) {
        ev.critical_diameter_m = k26astro_whipple_critical_diameter(
            target_wall_thickness_m,
            target_bumper_density_kg_per_m3,
            projectile->projectile_density_kg_per_m3,
            impact_velocity_m_per_s,
            cos_impact_angle,
            target_bumper_spacing_m,
            target_wall_yield_stress_ksi);
        ev.penetrates = (ev.critical_diameter_m > 0.0 &&
                         projectile->projectile_diameter_m >
                         ev.critical_diameter_m) ? 1 : 0;
    }

    if (target_brinell_hardness > 0.0 && target_density_kg_per_m3 > 0.0 &&
        target_speed_of_sound_m_per_s > 0.0) {
        ev.monolithic_penetration_m = k26astro_monolithic_penetration_depth(
            projectile->projectile_diameter_m,
            target_brinell_hardness,
            projectile->projectile_density_kg_per_m3,
            target_density_kg_per_m3,
            impact_velocity_m_per_s,
            cos_impact_angle,
            target_speed_of_sound_m_per_s);
    }

    return ev;
}
