/* test_body_set_mass.c — gate for k26astro_body_set_mass.
 *
 * Validates the mass-mutator path used by spacecraft propellant
 * depletion (and by libk26astro_vehicle's stage-event handler):
 *
 *   (1) Round-trip: mass field stores the input.
 *   (2) GM derives as K26A_G * mass.
 *   (3) Negative inputs clamp to zero.
 *   (4) NULL body pointer is a no-op (no crash). */
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    K26AstroBody b;
    k26astro_body_init(&b);

    /* (1) + (2) Round-trip + derived GM. */
    k26astro_body_set_mass(&b, 1000.0);
    assert(b.mass == 1000.0);
    assert(b.gm   == K26A_G * 1000.0);

    /* Subsequent update overwrites both fields. */
    k26astro_body_set_mass(&b, 650.0);
    assert(b.mass == 650.0);
    assert(fabs(b.gm - K26A_G * 650.0) < 1e-25);

    /* (3) Negative clamps to zero — both fields. */
    k26astro_body_set_mass(&b, -1.0);
    assert(b.mass == 0.0);
    assert(b.gm   == 0.0);

    /* (4) NULL safe. */
    k26astro_body_set_mass(NULL, 1.0);

    printf("test_body_set_mass: OK\n");
    return 0;
}
