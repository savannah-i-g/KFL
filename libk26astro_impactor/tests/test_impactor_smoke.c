/* test_impactor_smoke.c — impactor opaque + slot binding. */

#include "k26astro_impactor/impactor.h"
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

    /* Single-shot impactor. */
    K26AstroImpactor *k_single = k26astro_impactor_new(
        K26ASTRO_IMPACTOR_PATTERN_SINGLE,
        /* mass */ 10.0,
        /* density */ K26ASTRO_IMPACTOR_RHO_TUNGSTEN_KG_PER_M3,
        /* diameter */ 0.05,
        /* count (ignored) */ 0,
        /* half_angle (ignored) */ 0.0);
    assert(k_single != NULL);
    assert(k26astro_impactor_pattern(k_single) == K26ASTRO_IMPACTOR_PATTERN_SINGLE);
    assert(k26astro_impactor_swarm_count(k_single) == 1);
    assert(k26astro_payload_kind(k26astro_impactor_as_payload(k_single)) ==
           K26ASTRO_DEFENSE_KIND_IMPACTOR);

    /* Swarm pattern. */
    K26AstroImpactor *k_swarm = k26astro_impactor_new(
        K26ASTRO_IMPACTOR_PATTERN_SWARM,
        /* total mass */ 50.0,
        /* density */ K26ASTRO_IMPACTOR_RHO_STEEL_KG_PER_M3,
        /* diameter (each) */ 0.01,
        /* count */ 100,
        /* half_angle */ 0.01);
    assert(k_swarm != NULL);
    assert(k26astro_impactor_pattern(k_swarm) == K26ASTRO_IMPACTOR_PATTERN_SWARM);
    assert(k26astro_impactor_swarm_count(k_swarm) == 100);

    /* Bind both to vehicle. */
    assert(k26astro_impactor_bind(k_single, v) == 0);
    assert(k26astro_impactor_bind(k_swarm, v) == 0);

    /* Invalid: swarm pattern needs positive count. */
    assert(k26astro_impactor_new(K26ASTRO_IMPACTOR_PATTERN_SWARM,
                             50.0, 7850.0, 0.01, 0, 0.01) == NULL);
    /* Invalid: negative mass. */
    assert(k26astro_impactor_new(K26ASTRO_IMPACTOR_PATTERN_SINGLE,
                             -1.0, 7850.0, 0.05, 0, 0.0) == NULL);
    /* Invalid pattern enum. */
    assert(k26astro_impactor_new(99, 10.0, 7850.0, 0.05, 0, 0.0) == NULL);

    /* Sensor-first destruction. */
    k26astro_impactor_destroy(k_single);
    k26astro_impactor_destroy(k_swarm);
    k26astro_vehicle_destroy(v);

    /* Reverse order. */
    K26AstroBody body2;
    memset(&body2, 0, sizeof(body2));
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(v2, &body2);
    K26AstroImpactor *k3 = k26astro_impactor_new(
        K26ASTRO_IMPACTOR_PATTERN_SINGLE, 5.0, 7850.0, 0.03, 0, 0.0);
    assert(k26astro_impactor_bind(k3, v2) == 0);
    k26astro_vehicle_destroy(v2);
    k26astro_impactor_destroy(k3);

    /* NULL safety. */
    k26astro_impactor_destroy(NULL);
    assert(k26astro_impactor_as_payload(NULL) == NULL);

    printf("test_impactor_smoke: OK\n");
    return 0;
}
