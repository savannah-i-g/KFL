/* perturb_srp.c — solar radiation pressure with conical shadow.
 *
 * Solar photons impart momentum on any body intercepting them. The
 * acceleration on a body with cross-section A and mass m at
 * heliocentric distance r is:
 *
 *   a_SRP = (L_sun / (4π c r²)) · (A/m) · C_r · r̂
 *
 * where C_r is the radiation-pressure coefficient (~1.0 for a perfect
 * absorber, ~2.0 for a perfect specular reflector, ~1.3 for typical
 * spacecraft surfaces).
 *
 * v0.1 simplifications:
 *  - Central body (body[0]) must be a STAR. Its luminosity is read
 *    from K26A_L_SUN (the Sun's luminosity). For non-solar stars we'd
 *    need a per-star luminosity field; future work.
 *  - Per-body SRP parameters (A/m and C_r) are stored in a parallel
 *    table passed via the perturbation context. The v0.1 default —
 *    when ctx is NULL — uses (A/m = 0.01 m²/kg, C_r = 1.3) for any
 *    body with gm < 1e12 m³/s² (i.e., spacecraft/asteroid scale,
 *    not planets/moons).
 *  - Conical shadow: a body is in shadow if any other body of kind
 *    PLANET/MOON sits between it and the Sun within the geometric
 *    umbra cone. Linear test: skip if the body's heliocentric position
 *    projects through any planet's center within radius R.
 *
 * Reference: Vallado §8.6.4; McMahon & Scheeres (2010) for shadow
 * geometry. */
#include "k26astro_grav/perturb.h"
#include "k26astro_grav/grav.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <math.h>

/* Default SRP parameters when no ctx is supplied. */
#define K26_SRP_DEFAULT_AM   0.01    /* m²/kg, typical spacecraft */
#define K26_SRP_DEFAULT_CR   1.3
#define K26_SRP_BODY_GM_MAX  1.0e12  /* bodies above this don't get SRP */

static int in_shadow_(const K26AstroGravView *view, int sun_idx, int body_idx)
{
    /* Pop quiz: is body_idx in the umbra of any other body, relative
     * to the Sun? For each candidate occluder, compute its angular
     * radius from the body's vantage point and check if the Sun's
     * direction lies within it. */
    const K26AstroBody *sun  = &view->bodies[sun_idx];
    const K26AstroBody *body = &view->bodies[body_idx];

    K26V3 to_sun = k26astro_pos_sub(&sun->pos, &body->pos);
    double r_sun = sqrt(to_sun.x*to_sun.x + to_sun.y*to_sun.y + to_sun.z*to_sun.z);
    if (r_sun <= 0.0) return 0;
    /* Unit vector to Sun. */
    K26V3 dir_sun = { to_sun.x/r_sun, to_sun.y/r_sun, to_sun.z/r_sun };

    for (int j = 0; j < view->n; j++) {
        if (j == sun_idx || j == body_idx) continue;
        const K26AstroBody *occ = &view->bodies[j];
        if (occ->kind != K26ASTRO_BODY_PLANET
         && occ->kind != K26ASTRO_BODY_MOON) continue;
        if (occ->radius <= 0.0) continue;

        K26V3 to_occ = k26astro_pos_sub(&occ->pos, &body->pos);
        double r_occ = sqrt(to_occ.x*to_occ.x + to_occ.y*to_occ.y + to_occ.z*to_occ.z);
        if (r_occ <= 0.0) continue;

        /* Occluder is behind the body? */
        double dot = (to_occ.x*dir_sun.x + to_occ.y*dir_sun.y + to_occ.z*dir_sun.z);
        if (dot < 0.0) continue;
        if (dot > r_sun) continue;  /* occluder past the Sun */

        /* Perpendicular distance from body→Sun line to occluder. */
        K26V3 closest = { dir_sun.x * dot, dir_sun.y * dot, dir_sun.z * dot };
        double dx = to_occ.x - closest.x;
        double dy = to_occ.y - closest.y;
        double dz = to_occ.z - closest.z;
        double perp = sqrt(dx*dx + dy*dy + dz*dz);

        if (perp < occ->radius) {
            return 1;   /* in umbra */
        }
    }
    return 0;
}

int k26astro_srp_shadow_test(const K26AstroGravView *view,
                              int sun_idx, int body_idx)
{
    return in_shadow_(view, sun_idx, body_idx);
}

void k26astro_perturb_srp(const K26AstroGravState *state,
                          const K26AstroGravView  *view,
                          K26V3 *accel_out, void *ctx)
{
    (void)ctx;
    if (!state || !view || !accel_out) return;
    if (view->n < 2) return;

    const K26AstroBody *sun = &view->bodies[0];
    if (sun->kind != K26ASTRO_BODY_STAR) return;

    /* Pre-compute the SRP scalar = L_sun · C_r / (4π c). For non-solar
     * stars this would scale with the star's L, but v0.1 hardwires
     * the Sun. */
    const double L_sun = K26A_L_SUN;
    const double c_light = K26A_C;
    const double four_pi_c = 4.0 * K26A_PI * c_light;

    for (int i = 1; i < view->n; i++) {
        const K26AstroBody *bi = &view->bodies[i];
        if (bi->gm > K26_SRP_BODY_GM_MAX) continue;
        if (bi->mass <= 0.0) continue;

        if (in_shadow_(view, 0, i)) continue;

        K26V3 r = k26astro_pos_sub(&bi->pos, &sun->pos);
        double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
        if (r2 <= 0.0) continue;
        double rmag = sqrt(r2);

        double A_over_m = K26_SRP_DEFAULT_AM;
        double Cr = K26_SRP_DEFAULT_CR;
        double a_mag = (L_sun * Cr) / (four_pi_c * r2) * A_over_m;

        /* Direction: away from Sun (radiation pushes outward). */
        accel_out[i].x += a_mag * r.x / rmag;
        accel_out[i].y += a_mag * r.y / rmag;
        accel_out[i].z += a_mag * r.z / rmag;
    }
}
