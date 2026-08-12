/* mirror_strehl.c — Maréchal Strehl ratio for wavefront error. */

#include "k26astro_laser/mirror_strehl.h"

#include "k26astro_core/consts.h"

#include <math.h>

double k26astro_mirror_strehl(double rms_wavefront_m,
                              double wavelength_m)
{
    if (wavelength_m <= 0.0) return 1.0;
    if (rms_wavefront_m <= 0.0) return 1.0;
    const double sigma_waves = rms_wavefront_m / wavelength_m;
    const double arg = 2.0 * K26A_PI * sigma_waves;
    return exp(-arg * arg);
}
