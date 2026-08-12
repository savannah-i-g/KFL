/* plasma_plug.c — ablation-plume self-attenuation. */

#include "k26astro_laser/plasma_plug.h"

#include <math.h>

double k26astro_plasma_plug_transmissivity(double fluence_J_per_m2,
                                            double threshold_J_per_m2,
                                            double k_attenuation)
{
    if (threshold_J_per_m2 <= 0.0) return 1.0;
    if (fluence_J_per_m2 < threshold_J_per_m2) return 1.0;
    if (k_attenuation <= 0.0) return 1.0;

    const double x = fluence_J_per_m2 / threshold_J_per_m2 - 1.0;
    const double T = exp(-k_attenuation * x);
    if (T < 0.0) return 0.0;
    if (T > 1.0) return 1.0;
    return T;
}
