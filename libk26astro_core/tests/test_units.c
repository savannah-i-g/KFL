/* test_units.c — unit-conversion round-trips. */
#include "k26astro_core/units.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int near_(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void)
{
    /* Length conversions, round-trip. */
    assert(near_(k26a_to_au(k26a_au(1.5)), 1.5, 1e-14));
    assert(near_(k26a_to_pc(k26a_pc(3.26)), 3.26, 1e-14));
    assert(near_(k26a_to_ly(k26a_ly(4.2)), 4.2, 1e-14));
    assert(near_(k26a_to_light_seconds(k26a_light_seconds(8.0)), 8.0, 1e-14));

    /* AU value matches the IAU 2012 definition. */
    assert(k26a_au(1.0) == K26A_AU_M);
    assert(K26A_AU_M == 149597870700.0);

    /* Light-second: c metres. */
    assert(k26a_light_seconds(1.0) == K26A_C);

    /* Parsec via the AU + arcsec definition (within float precision). */
    double pc_expected = K26A_AU_M / K26A_RAD_PER_ARCSEC;
    assert(near_(k26a_pc(1.0), pc_expected, 1e-13));

    /* Time conversions. */
    assert(k26a_day(1.0) == 86400.0);
    assert(k26a_hr(1.0)  == 3600.0);
    assert(near_(k26a_to_year(k26a_year(2.5)), 2.5, 1e-14));

    /* Angle conversions. */
    assert(near_(k26a_to_deg(k26a_deg(45.0)),       45.0,  1e-14));
    assert(near_(k26a_to_arcsec(k26a_arcsec(3600.0)), 3600.0, 1e-12));
    assert(near_(k26a_to_mas(k26a_mas(100.0)),       100.0, 1e-12));

    /* π radians is 180 degrees. */
    assert(near_(k26a_to_deg(K26A_PI), 180.0, 1e-13));

    /* Mass via GM ratio. */
    /* GM_sun in solar units = K26A_GM_SUN; in kg via /G. */
    double m_sun_kg = k26a_msun(1.0);
    assert(m_sun_kg > 1.9e30 && m_sun_kg < 2.0e30);
    assert(near_(k26a_to_msun(m_sun_kg), 1.0, 1e-14));

    /* Sanity check on speed of light. */
    assert(K26A_C == 299792458.0);

    printf("test_units: OK\n");
    return 0;
}
