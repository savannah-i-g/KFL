/* test_vehicle_smoke.c — basic construction + scalar setter /
 * query round-trip gate for libk26astro_vehicle.
 *
 * Validates:
 *   (1) k26astro_vehicle_new returns non-NULL and starts with mass
 *       and MGA at zero, identity inertia, zero COM.
 *   (2) Mass setters round-trip and clamp negatives to zero.
 *   (3) MGA setter feeds into predicted_mass and mev_mass.
 *   (4) Diagonal inertia setter round-trips and zeroes off-diagonals.
 *   (5) Full inertia setter preserves off-diagonals.
 *   (6) Inertia inverse is recomputed by the underlying Ext path.
 *   (7) COM setter round-trips.
 *   (8) bind_body / body accessor.
 *   (9) Destroy is clean on a populated vehicle.
 *  (10) Every public setter is NULL-safe. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_core/epoch.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* (1) Construction defaults. */
    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);
    assert(k26astro_vehicle_mass_now(v) == 0.0);
    assert(k26astro_vehicle_predicted_mass(v) == 0.0);
    assert(k26astro_vehicle_mev_mass(v) == 0.0);

    K26AstroEpoch t0 = k26astro_epoch_j2000_tt();
    K26M3 I0 = k26astro_vehicle_inertia_at(v, t0);
    assert(I0.m[0][0] == 1.0 && I0.m[1][1] == 1.0 && I0.m[2][2] == 1.0);
    assert(I0.m[0][1] == 0.0 && I0.m[1][0] == 0.0);

    K26V3 com0 = k26astro_vehicle_com_at(v, t0);
    assert(com0.x == 0.0 && com0.y == 0.0 && com0.z == 0.0);

    /* (2) + (3) Mass + MGA setters. */
    k26astro_vehicle_set_dry_mass(v, 549054.0);
    k26astro_vehicle_set_mga    (v,  10000.0);
    assert(k26astro_vehicle_mass_now(v)       == 549054.0);
    assert(k26astro_vehicle_predicted_mass(v) == 559054.0);
    assert(k26astro_vehicle_mev_mass(v)       == 559054.0);

    k26astro_vehicle_set_dry_mass(v, -1.0);
    assert(k26astro_vehicle_mass_now(v) == 0.0);
    k26astro_vehicle_set_mga(v, -2.0);
    assert(k26astro_vehicle_predicted_mass(v) == 0.0);
    k26astro_vehicle_set_dry_mass(v, 1000.0);
    k26astro_vehicle_set_mga    (v,   50.0);

    /* (4) Diagonal inertia round-trip. */
    k26astro_vehicle_set_inertia_diag(v, 1000.0, 2000.0, 3000.0);
    K26M3 I1 = k26astro_vehicle_inertia_at(v, t0);
    assert(I1.m[0][0] == 1000.0);
    assert(I1.m[1][1] == 2000.0);
    assert(I1.m[2][2] == 3000.0);
    assert(I1.m[0][1] == 0.0);
    assert(I1.m[1][0] == 0.0);
    assert(I1.m[0][2] == 0.0);
    assert(I1.m[2][0] == 0.0);
    assert(I1.m[1][2] == 0.0);
    assert(I1.m[2][1] == 0.0);

    /* (6) The Ext path's inverse must be recomputed. */
    K26AstroAttitudeStateExt *a = k26astro_vehicle_attitude_ext(v);
    assert(a != NULL);
    assert(fabs(a->inertia_inverse.m[0][0] - 1.0 / 1000.0) < 1e-12);
    assert(fabs(a->inertia_inverse.m[1][1] - 1.0 / 2000.0) < 1e-12);
    assert(fabs(a->inertia_inverse.m[2][2] - 1.0 / 3000.0) < 1e-12);

    /* Negative diagonal rejected (existing values preserved). */
    k26astro_vehicle_set_inertia_diag(v, -1.0, 2000.0, 3000.0);
    K26M3 I1b = k26astro_vehicle_inertia_at(v, t0);
    assert(I1b.m[0][0] == 1000.0);

    /* (5) Full inertia preserves off-diagonals. */
    K26M3 J;
    J.m[0][0] = 1000.0; J.m[0][1] =  100.0; J.m[0][2] =   50.0;
    J.m[1][0] =  100.0; J.m[1][1] = 2000.0; J.m[1][2] =  -25.0;
    J.m[2][0] =   50.0; J.m[2][1] =  -25.0; J.m[2][2] = 3000.0;
    k26astro_vehicle_set_inertia_full(v, J);
    K26M3 I2 = k26astro_vehicle_inertia_at(v, t0);
    assert(I2.m[0][1] == 100.0);
    assert(I2.m[1][0] == 100.0);
    assert(I2.m[1][2] == -25.0);
    assert(I2.m[0][0] == 1000.0);
    assert(I2.m[1][1] == 2000.0);
    assert(I2.m[2][2] == 3000.0);

    /* (7) COM round-trip. */
    k26astro_vehicle_set_com_offset(v, 0.5, -0.25, -1.5);
    K26V3 com = k26astro_vehicle_com_at(v, t0);
    assert(com.x == 0.5);
    assert(com.y == -0.25);
    assert(com.z == -1.5);

    /* (8) bind_body. */
    K26AstroBody b;
    k26astro_body_init(&b);
    k26astro_vehicle_bind_body(v, &b);
    assert(k26astro_vehicle_body(v) == &b);
    /* Re-binding allowed. */
    k26astro_vehicle_bind_body(v, NULL);
    assert(k26astro_vehicle_body(v) == NULL);
    k26astro_vehicle_bind_body(v, &b);

    /* Internal accessors. */
    assert(k26astro_vehicle_mass_accum_get(v) == 0.0);
    k26astro_vehicle_mass_accum_add(v, -2.5);
    assert(k26astro_vehicle_mass_accum_get(v) == -2.5);
    k26astro_vehicle_mass_accum_clear(v);
    assert(k26astro_vehicle_mass_accum_get(v) == 0.0);

    /* (9) Destroy a populated vehicle — clean. */
    k26astro_vehicle_destroy(v);

    /* (10) NULL-safe setters + queries. */
    k26astro_vehicle_destroy(NULL);
    k26astro_vehicle_bind_body(NULL, &b);
    k26astro_vehicle_set_dry_mass(NULL, 100.0);
    k26astro_vehicle_set_mga(NULL, 10.0);
    k26astro_vehicle_set_inertia_diag(NULL, 1, 2, 3);
    K26M3 K = { .m = {{1,0,0},{0,1,0},{0,0,1}} };
    k26astro_vehicle_set_inertia_full(NULL, K);
    k26astro_vehicle_set_com_offset(NULL, 0, 0, 0);
    assert(k26astro_vehicle_mass_now(NULL) == 0.0);
    assert(k26astro_vehicle_predicted_mass(NULL) == 0.0);
    assert(k26astro_vehicle_mev_mass(NULL) == 0.0);
    assert(k26astro_vehicle_mass_at(NULL, t0) == 0.0);
    K26V3 com_null = k26astro_vehicle_com_at(NULL, t0);
    assert(com_null.x == 0.0 && com_null.y == 0.0 && com_null.z == 0.0);
    K26M3 I_null = k26astro_vehicle_inertia_at(NULL, t0);
    assert(I_null.m[0][0] == 1.0 && I_null.m[1][1] == 1.0 && I_null.m[2][2] == 1.0);
    assert(k26astro_vehicle_attitude_ext(NULL) == NULL);
    assert(k26astro_vehicle_body(NULL) == NULL);
    assert(k26astro_vehicle_generation(NULL) == 0);
    k26astro_vehicle_commit_mass_step(NULL, 1.0);

    printf("test_vehicle_smoke: OK\n");
    return 0;
}
