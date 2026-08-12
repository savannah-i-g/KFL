/* terminal_phase.c — closing-velocity + impact-geometry math. */

#include "k26astro_impactor/terminal_phase.h"

#include <math.h>

double k26astro_impactor_closing_speed(K26V3 v_projectile, K26V3 v_target)
{
    K26V3 rel = k26m3d_v3_sub(v_projectile, v_target);
    return k26m3d_v3_len(rel);
}

double k26astro_impactor_time_to_closest_approach(
    K26V3 r_projectile, K26V3 v_projectile,
    K26V3 r_target,     K26V3 v_target)
{
    /* Relative state. */
    K26V3 dr = k26m3d_v3_sub(r_projectile, r_target);
    K26V3 dv = k26m3d_v3_sub(v_projectile, v_target);
    const double v2 = k26m3d_v3_dot(dv, dv);
    if (v2 <= 0.0) return 0.0;
    return -k26m3d_v3_dot(dr, dv) / v2;
}

double k26astro_impactor_closest_approach_distance(
    K26V3 r_projectile, K26V3 v_projectile,
    K26V3 r_target,     K26V3 v_target)
{
    K26V3 dr = k26m3d_v3_sub(r_projectile, r_target);
    K26V3 dv = k26m3d_v3_sub(v_projectile, v_target);
    const double v2 = k26m3d_v3_dot(dv, dv);
    if (v2 <= 0.0) return k26m3d_v3_len(dr);
    const double t = -k26m3d_v3_dot(dr, dv) / v2;
    K26V3 closest = k26m3d_v3_add(dr, k26m3d_v3_scale(dv, t));
    return k26m3d_v3_len(closest);
}

double k26astro_impactor_impact_cos_angle(K26V3 closing_velocity_unit,
                                      K26V3 surface_normal_unit)
{
    /* Angle between -closing velocity (the incoming projectile's
     * direction from target's frame) and the outward surface
     * normal. */
    K26V3 neg = k26m3d_v3_neg(closing_velocity_unit);
    double c = k26m3d_v3_dot(neg, surface_normal_unit);
    if (c < 0.0) c = 0.0;
    if (c > 1.0) c = 1.0;
    return c;
}
