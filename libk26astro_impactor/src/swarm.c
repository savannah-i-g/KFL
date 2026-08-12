/* swarm.c — swarm dispatch geometry.
 *
 * Samples uniformly on a spherical cap's solid angle:
 *   φ ~ U(0, 2π), cos θ ~ U(cos α, 1) where α is the cap half-angle.
 * This is the standard cap-sampling result; uniform-cos gives
 * a uniform distribution on the cap's solid-angle measure. */

#include "k26astro_impactor/swarm.h"

#include "k26astro_core/consts.h"

#include <math.h>

/* Build a unit vector that is θ degrees off `axis` with the
 * azimuth pinned at φ. The axis is rotated to a canonical (0,0,1)
 * frame, the cone direction is built in that frame, then rotated
 * back. */
static K26V3 rotate_cone_axis_to_z_back_(K26V3 axis_unit,
                                          K26V3 cone_dir_in_z)
{
    /* If axis is already very close to +z, no rotation needed. */
    const double tol = 1.0e-12;
    if (fabs(axis_unit.x) < tol && fabs(axis_unit.y) < tol) {
        if (axis_unit.z >= 0.0) {
            return cone_dir_in_z;
        }
        /* axis is -z: flip cone_dir's z component. */
        K26V3 r = cone_dir_in_z;
        r.z = -r.z;
        return r;
    }
    /* Build the rotation that maps +z to axis_unit using the
     * Rodrigues formula with axis = z × axis_unit / |.|, angle =
     * acos(z · axis_unit). */
    const double cos_a = axis_unit.z;
    /* sin_a >= 0 since we'll use unit cross-product magnitude. */
    const double sin_a = sqrt(1.0 - cos_a * cos_a);
    /* k = (z × axis) / sin_a. (-axis.y, axis.x, 0) / sin_a. */
    K26V3 k = { -axis_unit.y / sin_a, axis_unit.x / sin_a, 0.0 };
    /* Rodrigues: v_rot = v cos_a + (k × v) sin_a + k (k·v) (1-cos_a). */
    K26V3 v = cone_dir_in_z;
    K26V3 kxv = k26m3d_v3_cross(k, v);
    double kdv = k26m3d_v3_dot(k, v);
    K26V3 term1 = k26m3d_v3_scale(v, cos_a);
    K26V3 term2 = k26m3d_v3_scale(kxv, sin_a);
    K26V3 term3 = k26m3d_v3_scale(k, kdv * (1.0 - cos_a));
    return k26m3d_v3_add(term1, k26m3d_v3_add(term2, term3));
}

K26V3 k26astro_swarm_sample_direction(K26V3 cone_axis_unit,
                                           double half_angle_rad,
                                           K26CRng *rng)
{
    if (!rng || half_angle_rad <= 0.0) return cone_axis_unit;

    /* Sample uniformly on the cap. */
    const double u_cos = k26c_rng_uniform(rng);
    const double u_phi = k26c_rng_uniform(rng);
    const double cos_alpha = cos(half_angle_rad);
    /* cos θ uniform in [cos α, 1]. */
    const double cos_theta = cos_alpha + (1.0 - cos_alpha) * u_cos;
    const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    const double phi = K26A_TWO_PI * u_phi;

    /* Direction in the canonical (axis = +z) frame. */
    K26V3 d_z = {
        sin_theta * cos(phi),
        sin_theta * sin(phi),
        cos_theta
    };
    return rotate_cone_axis_to_z_back_(cone_axis_unit, d_z);
}

double k26astro_swarm_per_unit_mass(int count, double total_mass_kg)
{
    if (count <= 0 || total_mass_kg <= 0.0) return 0.0;
    return total_mass_kg / (double)count;
}

double k26astro_swarm_footprint_m2(double range_m,
                                              double half_angle_rad)
{
    if (range_m <= 0.0 || half_angle_rad <= 0.0) return 0.0;
    const double r = range_m * sin(half_angle_rad);
    return K26A_PI * r * r;
}
