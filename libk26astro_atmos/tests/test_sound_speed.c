/* test_sound_speed.c — speed-of-sound accessor gate.
 *
 * Validates:
 *   (1) Earth ISA sea-level sound speed matches the textbook
 *       340.29 m/s within 0.1 m/s.
 *   (2) Isothermal model returns the same value at 1 km, 10 km,
 *       50 km (within the barometric region).
 *   (3) Above the atmos_top_m clip, returns 0.
 *   (4) Zero-initialised AtmosParams (no T₀/γ/R supplied) returns 0
 *       without crashing — designated-initializer callers that
 *       skip the new fields get a defined zero result.
 *   (5) NULL handle returns 0. */
#include "k26astro_atmos/atmos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* ---- (1) ISA sea level: 340.29 m/s ± 0.1 ------------------ */
    K26AstroAtmos *a = k26astro_atmos_earth_standard();
    assert(a != NULL);
    double c0 = k26astro_atmos_sound_speed_at(a, 0.0);
    fprintf(stderr, "test_sound_speed: c(0) = %.4f m/s\n", c0);
    assert(fabs(c0 - 340.29) < 0.1);

    /* ---- (2) Isothermal: same value at 1, 10, 50 km ----------- */
    double c1   = k26astro_atmos_sound_speed_at(a, 1.0e3);
    double c10  = k26astro_atmos_sound_speed_at(a, 1.0e4);
    double c50  = k26astro_atmos_sound_speed_at(a, 5.0e4);
    assert(fabs(c1  - c0) < 1.0e-12);
    assert(fabs(c10 - c0) < 1.0e-12);
    assert(fabs(c50 - c0) < 1.0e-12);

    /* ---- (3) Above atmos top: 0 ------------------------------- */
    double c200 = k26astro_atmos_sound_speed_at(a, 2.0e5);
    assert(c200 == 0.0);

    k26astro_atmos_destroy(a);

    /* ---- (4) Zero AtmosParams: 0 without crash ---------------- */
    K26AstroAtmosParams zero = {
        .scale_height_m     = 8400.0,
        .rayleigh_coeff     = 1.0e-5,
        .mie_coeff          = 1.0e-6,
        .ground_pressure_pa = 101325.0,
        .mie_asymmetry_g    = 0.5,
        .n0                 = 2.93e-4,
        .atmos_top_m        = 100000.0
        /* t_0_k, gamma_ratio_specific_heats, r_specific_j_per_kg_k
         * deliberately omitted — designated-initializer default to 0. */
    };
    K26AstroAtmos *bare = k26astro_atmos_init(&zero);
    assert(bare != NULL);
    double c_bare = k26astro_atmos_sound_speed_at(bare, 0.0);
    assert(c_bare == 0.0);
    k26astro_atmos_destroy(bare);

    /* ---- (5) NULL handle ------------------------------------- */
    assert(k26astro_atmos_sound_speed_at(NULL, 0.0) == 0.0);

    printf("test_sound_speed: OK\n");
    return 0;
}
