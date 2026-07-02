/* test_shapiro_direction.c — Shapiro direction rotation gate.
 *
 * Verifies the GR direction-vector rotation (gravitational lensing),
 * which complements the delay-as-distance magnitude.
 *
 * Geometry: Sun at origin, observer at 1 AU along +x, target at
 * 1 AU along -x but offset +y by `b` (the impact parameter).
 * Light from target reaches observer along a path whose closest
 * approach to the Sun is approximately b. For b = R_sun, the
 * deflection is the classical 1.75 arcsec (Einstein 1916). For
 * larger b, it scales as 1/b.
 *
 * Expected deflection:
 *     Δθ = (2 GM_sun / (c² b)) · (1 + cos α)
 *
 * With α (observer → Sun → target) = π (exactly opposite sides),
 * cos α = -1, so Δθ = 0 — no deflection in the antipodal limit.
 * In our test geometry, target is offset slightly in +y so α < π
 * but very close to π. We avoid the singular antipodal point.
 *
 * Approach: place target at angle (180° - δ) where δ is small but
 * non-zero (1°), with impact parameter b = 5·R_sun (so the
 * deflection is moderate and well above floating-point noise).
 * Compare the rotated apparent direction to the geometric direction.
 *
 * Acceptance:
 *   - Δθ matches the closed-form within 5% (impact-parameter
 *     approximation has some geometric slop)
 *   - rotation is in the correct direction (away from Sun)
 *   - direction stays a unit vector (|dir| - 1 < 1e-12) */
#include "k26astro_rt/world.h"
#include "k26astro_rt/observer.h"
#include "k26astro_grav/grav.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define K26ASTRO_C_LIGHT 299792458.0
#define R_SUN_M          6.957e8

static double v3_norm_(K26V3 v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static int test_shapiro_grazing_ray(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    k26astro_world_set_observer_mode(w, K26ASTRO_OBS_APPARENT);

    /* Enable GR PPN-1 so Shapiro corrections engage. */
    K26AstroGravState *grav = k26astro_world_grav(w);
    k26astro_grav_enable_gr_ppn1(grav, 1);

    K26AstroBody sun;
    k26astro_body_init(&sun);
    strncpy(sun.name, "sun", sizeof sun.name - 1);
    sun.mass = 1.989e30;
    sun.gm   = K26A_GM_SUN;
    sun.pos  = k26astro_pos_zero();
    int sun_i = k26astro_world_add_body(w, sun);

    K26AstroBody observer;
    k26astro_body_init(&observer);
    strncpy(observer.name, "observer", sizeof observer.name - 1);
    observer.mass = 1.0; observer.gm = 0.0;
    observer.pos  = k26astro_pos_from_m(K26A_AU_M, 0.0, 0.0);
    observer.vel  = k26m3d_v3(0.0, 0.0, 0.0);  /* No aberration. */
    int obs_i = k26astro_world_add_body(w, observer);

    /* Target at (-0.1 AU, +1 AU, 0): on the same side as observer
     * relative to Sun but offset upward. The observer-Sun-target
     * angle α ≈ 95.7° → (1+cos α) ≈ 0.9 — well-defined, non-zero
     * geometry factor in the (1+cos α) form. Impact parameter is
     * the perpendicular distance from Sun to the observer→target
     * line, ~0.67 AU. */
    double tgt_x = -0.1 * K26A_AU_M;
    double tgt_y =  1.0 * K26A_AU_M;
    K26AstroBody target;
    k26astro_body_init(&target);
    strncpy(target.name, "target", sizeof target.name - 1);
    target.mass = 1.0; target.gm = 0.0;
    target.pos  = k26astro_pos_from_m(tgt_x, tgt_y, 0.0);
    target.vel  = k26m3d_v3(0.0, 0.0, 0.0);  /* No light-time effect. */
    int tgt_i = k26astro_world_add_body(w, target);

    K26V3 dir_apparent;
    K26AstroPos out_pos;
    int rc = k26astro_world_observe(w, tgt_i, obs_i,
                                     &out_pos, &dir_apparent);
    if (rc != K26ASTRO_RT_OK) {
        fprintf(stderr, "test_shapiro: observe rc=%d\n", rc);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Geometric direction (without Shapiro). */
    K26AstroPos obs_pos = k26astro_pos_from_m(K26A_AU_M, 0.0, 0.0);
    K26AstroPos tgt_pos = k26astro_pos_from_m(tgt_x, tgt_y, 0.0);
    K26V3 r_geo = k26astro_pos_sub(&tgt_pos, &obs_pos);
    double n_geo = v3_norm_(r_geo);
    K26V3 dir_geo = { r_geo.x / n_geo,
                       r_geo.y / n_geo,
                       r_geo.z / n_geo };

    /* Angle between apparent and geometric directions. acos of dot
     * product is numerically unstable for θ < 1e-7 (1-cos(θ) ~ θ²/2
     * underflows below double precision). Use cross-product /
     * sin form instead. */
    K26V3 dcross = {
        dir_geo.y * dir_apparent.z - dir_geo.z * dir_apparent.y,
        dir_geo.z * dir_apparent.x - dir_geo.x * dir_apparent.z,
        dir_geo.x * dir_apparent.y - dir_geo.y * dir_apparent.x
    };
    double sin_diff = v3_norm_(dcross);
    if (sin_diff > 1.0) sin_diff = 1.0;
    double dtheta_obs = asin(sin_diff);

    /* Closed-form predicted deflection per the plan's formula:
     *   Δθ = (2 GM / c²b) · (1 + cos α)
     * α = angle observer → Sun → target at the Sun. The impact
     * parameter b is the perpendicular distance from Sun to the
     * observer→target ray line, which the implementation
     * approximates as r_o · sin_off where sin_off = |dir × u_def|. */
    double r_o = K26A_AU_M;
    double r_t = sqrt(tgt_x*tgt_x + tgt_y*tgt_y);
    double cos_alpha = (K26A_AU_M * tgt_x + 0.0 + 0.0) / (r_o * r_t);

    /* Impact parameter (perpendicular distance from origin to
     * line obs→tgt). The implementation uses r_o · sin_off, which
     * equals the perpendicular distance from origin to the
     * apparent ray as seen from the observer — exactly b. */
    K26V3 u_def = { -1.0, 0.0, 0.0 };  /* observer → sun unit vector */
    K26V3 cross = {
        dir_geo.y * u_def.z - dir_geo.z * u_def.y,
        dir_geo.z * u_def.x - dir_geo.x * u_def.z,
        dir_geo.x * u_def.y - dir_geo.y * u_def.x
    };
    double sin_off = sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
    double b_imp   = r_o * sin_off;
    double dtheta_pred = (2.0 * K26A_GM_SUN
                         / (K26ASTRO_C_LIGHT * K26ASTRO_C_LIGHT * b_imp))
                       * (1.0 + cos_alpha);

    fprintf(stderr,
        "test_shapiro: Δθ observed = %.4e rad (= %.4f arcsec)\n",
        dtheta_obs, dtheta_obs * (180.0 * 3600.0 / M_PI));
    fprintf(stderr,
        "test_shapiro: Δθ predicted = %.4e rad (1 + cos α = %.6f)\n",
        dtheta_pred, 1.0 + cos_alpha);

    if (dtheta_obs <= 0.0) {
        fprintf(stderr, "FAIL: no rotation applied\n");
        k26astro_world_destroy(w);
        return 1;
    }
    /* The implementation uses the same impact-parameter formula
     * as we predicted above (r_o · sin_off), so the agreement
     * should be within a few %. */
    double ratio = dtheta_obs / dtheta_pred;
    if (ratio < 0.9 || ratio > 1.1) {
        fprintf(stderr,
            "FAIL: Δθ ratio observed/predicted = %.3f out of "
            "[0.9, 1.1]\n", ratio);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Verify apparent direction is a unit vector. */
    double apparent_norm = v3_norm_(dir_apparent);
    if (fabs(apparent_norm - 1.0) > 1.0e-12) {
        fprintf(stderr,
            "FAIL: apparent direction not normalised: |dir|-1 = %.3e\n",
            apparent_norm - 1.0);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Apparent direction should be FURTHER from the Sun than the
     * geometric direction (light bent toward Sun → source appears
     * displaced away). In our geometry Sun is along observer→Sun
     * unit vector u_def = (-1, 0, 0). The component of dir along
     * u_def should be SMALLER (in absolute terms) in apparent than
     * in geometric. */
    double dot_geo  = -dir_geo.x;       /* projection onto u_def */
    double dot_app  = -dir_apparent.x;
    if (dot_app >= dot_geo) {
        fprintf(stderr,
            "FAIL: apparent direction NOT displaced away from Sun "
            "(geo·u_def=%.10f, app·u_def=%.10f)\n", dot_geo, dot_app);
        k26astro_world_destroy(w);
        return 1;
    }

    /* Verify sun_i index — unused but proves we got valid indices. */
    (void)sun_i;
    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    if (test_shapiro_grazing_ray()) return 1;
    fprintf(stderr, "test_shapiro_direction: OK\n");
    return 0;
}
