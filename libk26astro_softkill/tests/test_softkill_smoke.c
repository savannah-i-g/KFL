/* test_softkill_smoke.c — decoy + jammer opaque + slot binding. */

#include "k26astro_softkill/softkill.h"
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

    /* Decoys. */
    K26AstroDecoy *d_pas = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE,
        /* dry mass */ 10.0,
        /* deploy Δv */ 2.0,
        /* ir match */ 0.3,
        /* rcs match */ 0.2,
        /* accel match */ 0.1);
    assert(d_pas != NULL);
    assert(k26astro_decoy_mode(d_pas) == K26ASTRO_DECOY_MODE_PASSIVE);
    assert(k26astro_payload_kind(k26astro_decoy_as_payload(d_pas)) ==
           K26ASTRO_DEFENSE_KIND_DECOY);

    K26AstroDecoy *d_act = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_ACTIVE,
        /* dry mass */ 50.0, /* Δv */ 2.0,
        /* ir */ 0.9, /* rcs */ 0.7, /* accel */ 0.5);
    assert(d_act != NULL);

    /* Jammer. */
    K26AstroJammer *j = k26astro_jammer_new(
        K26ASTRO_JAMMER_MODE_NOISE,
        /* P_j */ 1.0e3,
        /* G_j */ 25.0,
        /* freq */ 10.0e9,
        /* B */ 1.0e6,
        /* snr_th */ 1.0);
    assert(j != NULL);
    assert(k26astro_jammer_mode(j) == K26ASTRO_JAMMER_MODE_NOISE);
    assert(k26astro_payload_kind(k26astro_jammer_as_payload(j)) ==
           K26ASTRO_DEFENSE_KIND_JAMMER);

    /* Bind to vehicle (mixed kinds in the same PAYLOAD slot). */
    assert(k26astro_decoy_bind(d_pas, v) == 0);
    assert(k26astro_decoy_bind(d_act, v) == 0);
    assert(k26astro_jammer_bind(j, v) == 0);

    /* Invalid constructors return NULL. */
    assert(k26astro_decoy_new(99, 10.0, 2.0, 0.3, 0.2, 0.1) == NULL);
    assert(k26astro_decoy_new(K26ASTRO_DECOY_MODE_PASSIVE, -1.0,
                              2.0, 0.3, 0.2, 0.1) == NULL);
    assert(k26astro_jammer_new(99, 1e3, 25.0, 1e9, 1e6, 1.0) == NULL);
    assert(k26astro_jammer_new(K26ASTRO_JAMMER_MODE_NOISE,
                                -1.0, 25.0, 1e9, 1e6, 1.0) == NULL);

    /* Payload-first destruction. */
    k26astro_decoy_destroy(d_pas);
    k26astro_decoy_destroy(d_act);
    k26astro_jammer_destroy(j);
    k26astro_vehicle_destroy(v);

    /* Reverse order — vehicle dies first. */
    K26AstroBody body2;
    memset(&body2, 0, sizeof(body2));
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(v2, &body2);
    K26AstroDecoy *d2 = k26astro_decoy_new(
        K26ASTRO_DECOY_MODE_PASSIVE, 5.0, 1.0, 0.2, 0.1, 0.1);
    assert(k26astro_decoy_bind(d2, v2) == 0);
    K26AstroJammer *j2 = k26astro_jammer_new(
        K26ASTRO_JAMMER_MODE_DECEPTION, 500.0, 20.0, 5e9, 5e5, 1.0);
    assert(k26astro_jammer_bind(j2, v2) == 0);
    k26astro_vehicle_destroy(v2);
    k26astro_decoy_destroy(d2);
    k26astro_jammer_destroy(j2);

    /* NULL safety. */
    k26astro_decoy_destroy(NULL);
    k26astro_jammer_destroy(NULL);
    assert(k26astro_decoy_as_payload(NULL) == NULL);
    assert(k26astro_jammer_as_payload(NULL) == NULL);

    printf("test_softkill_smoke: OK\n");
    return 0;
}
