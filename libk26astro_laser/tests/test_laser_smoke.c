/* test_laser_smoke.c — opaque + payload-slot binding lifecycle. */

#include "k26astro_laser/laser.h"

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"
#include "k26astro_body/body.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroBody body;
    memset(&body, 0, sizeof(body));
    K26AstroVehicle *v = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(v, &body);

    /* Construct a laser. */
    K26AstroLaser *l = k26astro_laser_new(
        /* primary_diam_m */ 1.0,
        /* wavelength_nm */ 1064.0,
        /* p_output_W */ 1.0e6,
        /* m_squared */ 1.2,
        /* pointing_jitter_rad */ 1.0e-6,
        /* rms_wavefront_m */ 1064.0e-9 / 20.0,
        /* plasma_attn_k */ 0.5);
    assert(l != NULL);
    assert(k26astro_payload_kind(k26astro_laser_as_payload(l)) ==
           K26ASTRO_DEFENSE_KIND_LASER);

    /* Bind to vehicle. */
    assert(k26astro_laser_bind(l, v) == 0);

    /* Invalid construction returns NULL. */
    assert(k26astro_laser_new(-1.0, 1064.0, 1e6, 1.2, 1e-6, 5e-8, 0.5) == NULL);
    assert(k26astro_laser_new(1.0, -1.0, 1e6, 1.2, 1e-6, 5e-8, 0.5) == NULL);
    assert(k26astro_laser_new(1.0, 1064.0, -1e6, 1.2, 1e-6, 5e-8, 0.5) == NULL);

    /* Bind a second laser to the same vehicle (multi-payload mount). */
    K26AstroLaser *l2 = k26astro_laser_new(
        0.5, 532.0, 5.0e5, 1.3, 2.0e-6, 532.0e-9 / 15.0, 0.5);
    assert(l2 != NULL);
    assert(k26astro_laser_bind(l2, v) == 0);

    /* Sensor-first destruction. */
    k26astro_laser_destroy(l);
    k26astro_laser_destroy(l2);

    /* Vehicle teardown should walk dead slots without issue. */
    k26astro_vehicle_destroy(v);

    /* Reverse order. */
    K26AstroBody body2;
    memset(&body2, 0, sizeof(body2));
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(v2, &body2);

    K26AstroLaser *l3 = k26astro_laser_new(
        1.0, 1064.0, 1e6, 1.2, 1e-6, 5e-8, 0.5);
    assert(k26astro_laser_bind(l3, v2) == 0);
    /* Vehicle dies first — the on_owner_destroy callback in the
     * payload-base hook should null l3's owner back-reference. */
    k26astro_vehicle_destroy(v2);
    k26astro_laser_destroy(l3);

    /* NULL safety. */
    k26astro_laser_destroy(NULL);
    assert(k26astro_laser_as_payload(NULL) == NULL);

    printf("test_laser_smoke: OK\n");
    return 0;
}
