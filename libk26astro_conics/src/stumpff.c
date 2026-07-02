/* stumpff.c — universal-variable Stumpff functions.
 *
 * C(z) =  (1 - cos√z) / z          for z > 0  (elliptic)
 *      =  (cosh√(-z) - 1) / (-z)    for z < 0  (hyperbolic)
 *      ≈  1/2 - z/24 + z²/720       for |z| ≪ 1
 *
 * S(z) =  (√z - sin√z) / z^(3/2)
 *      =  (sinh√(-z) - √(-z)) / (-z)^(3/2)
 *      ≈  1/6 - z/120 + z²/5040
 *
 * The Maclaurin band (|z| < 1e-6) avoids cancellation in the
 * (1 - cos√z) / z evaluation as z → 0; both summands of the numerator
 * approach 1 there. */
#include "stumpff_internal.h"

#include <math.h>

double k26astro_conics_stumpff_C(double z)
{
    if (z > 1.0e-6) {
        double sqrt_z = sqrt(z);
        return (1.0 - cos(sqrt_z)) / z;
    }
    if (z < -1.0e-6) {
        double sqrt_mz = sqrt(-z);
        return (cosh(sqrt_mz) - 1.0) / (-z);
    }
    return 0.5 - z / 24.0 + z * z / 720.0;
}

double k26astro_conics_stumpff_S(double z)
{
    if (z > 1.0e-6) {
        double sqrt_z = sqrt(z);
        return (sqrt_z - sin(sqrt_z)) / (sqrt_z * sqrt_z * sqrt_z);
    }
    if (z < -1.0e-6) {
        double sqrt_mz = sqrt(-z);
        return (sinh(sqrt_mz) - sqrt_mz) / (sqrt_mz * sqrt_mz * sqrt_mz);
    }
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
}
