/* soi_radius.c — Hill and Laplace SOI radius calculators. */
#include "k26astro_conics/soi.h"

#include <math.h>

double k26astro_soi_radius_hill(double parent_gm, double child_gm,
                                double semi_major_axis)
{
    if (parent_gm <= 0.0 || child_gm <= 0.0 || semi_major_axis <= 0.0)
        return 0.0;
    return semi_major_axis * cbrt(child_gm / (3.0 * parent_gm));
}

double k26astro_soi_radius_laplace(double parent_gm, double child_gm,
                                   double semi_major_axis)
{
    if (parent_gm <= 0.0 || child_gm <= 0.0 || semi_major_axis <= 0.0)
        return 0.0;
    return semi_major_axis * pow(child_gm / parent_gm, 0.4);
}
