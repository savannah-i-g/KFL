/* test_perturb_j2.c — J2 oblateness perturbation.
 *
 * Acceptance:
 *   - J2 disabled: no perturbation acceleration added.
 *   - J2 enabled on a body at z=0 (equator): radial perturbation
 *     only, no z-component.
 *   - J2 enabled on a polar body (z != 0): non-trivial z-component
 *     in the perturbation.
 *   - Newton's 3rd law: central body picks up an equal-and-opposite
 *     reaction (scaled by mass ratio).
 */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_grav/perturb.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroBody bodies[2];
    memset(bodies, 0, sizeof(bodies));

    /* Earth + LEO satellite-like body. */
    bodies[0].kind = K26ASTRO_BODY_PLANET;
    bodies[0].gm   = K26A_GM_EARTH;
    bodies[0].radius = K26A_R_EARTH_EQU;
    bodies[0].j2   = 1.0826e-3;      /* WGS84 */
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = k26m3d_v3(0.0, 0.0, 0.0);
    bodies[0].parent_body_idx = -1;

    /* Satellite at 800 km altitude, on equator. */
    bodies[1].kind = K26ASTRO_BODY_SPACECRAFT;
    bodies[1].gm   = 0.0;
    bodies[1].mass = 1000.0;
    bodies[1].pos  = k26astro_pos_from_m(K26A_R_EARTH_EQU + 800e3, 0.0, 0.0);
    bodies[1].vel  = k26m3d_v3(0.0, 7400.0, 0.0);   /* ~LEO speed */
    bodies[1].parent_body_idx = 0;

    K26AstroGravView view = k26astro_grav_view(bodies, 2);
    K26V3 accel[2];

    /* J2 on equator (z=0): expect no z-component. */
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_j2(NULL, &view, accel, NULL);
    /* state arg is unused in current impl; pass NULL. */
    /* Wait — function reads state for null guard but doesn't use it
     * if not null. Let me supply a minimal state. */
    K26AstroGravState state;
    state.bodies = bodies;
    state.n_bodies = 2;
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_j2(&state, &view, accel, NULL);

    /* On equator: |z| acceleration should be tiny. x/y radial. */
    assert(fabs(accel[1].z) < 1e-12);
    /* Radial perturbation magnitude scales as GM·J2·R²/r⁵. */
    double r = K26A_R_EARTH_EQU + 800e3;
    double expected_mag = 1.5 * K26A_GM_EARTH * bodies[0].j2
                        * K26A_R_EARTH_EQU * K26A_R_EARTH_EQU
                        / pow(r, 4);
    /* At equator with z=0, formula: ax = -common·r·1, |a| = common·r. */
    /* Actually re-derive: common = (3/2)·GM·J2·R²/r⁵.
     * a_x = common * x * (5*z²/r² - 1) = common * r * (-1) = -common*r at z=0. */
    double a_mag = sqrt(accel[1].x*accel[1].x + accel[1].y*accel[1].y);
    fprintf(stderr, "test_perturb_j2: |a_J2| = %.3e (expected ~%.3e)\n",
            a_mag, expected_mag);
    assert(fabs(a_mag - expected_mag) / expected_mag < 0.01);

    /* Polar body (z = r): expect z-acceleration. */
    bodies[1].pos = k26astro_pos_from_m(0.0, 0.0, r);
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_j2(&state, &view, accel, NULL);
    /* At pure pole: x = y = 0, z = r.
     * a_z = common * r * (5 - 3) = 2 * common * r. */
    assert(fabs(accel[1].z) > 1e-8);  /* non-trivial */
    assert(fabs(accel[1].x) < 1e-12); /* no horizontal at pole */

    /* Reaction force on central body should be opposite. */
    /* Mass ratio: bi->gm / GM_central. Satellite has gm=0 so ratio=0
     * → no reaction. Let's give it a mass and check. */
    bodies[1].gm = 1e6;   /* artificial — make it noticeable */
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_j2(&state, &view, accel, NULL);
    /* Reaction on central in z: opposite sign, scaled by mass ratio. */
    double ratio = bodies[1].gm / bodies[0].gm;
    double expected_reaction = -accel[1].z * ratio;
    assert(fabs(accel[0].z - expected_reaction) / fabs(expected_reaction) < 1e-9);

    /* J2 = 0 → no perturbation. */
    bodies[0].j2 = 0.0;
    bodies[1].gm = 0.0;
    memset(accel, 0, sizeof(accel));
    k26astro_perturb_j2(&state, &view, accel, NULL);
    assert(accel[1].x == 0.0 && accel[1].y == 0.0 && accel[1].z == 0.0);

    printf("test_perturb_j2: OK\n");
    return 0;
}
