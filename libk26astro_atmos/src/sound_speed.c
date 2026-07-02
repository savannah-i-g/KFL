/* sound_speed.c — isothermal sound-speed accessor.
 *
 * The barometric density profile is isothermal: T(h) = T_0 across
 * the barometric region. The ideal-gas adiabatic sound speed
 *
 *   c = sqrt(γ · R_specific · T)
 *
 * is therefore constant in h within that region. The signature
 * carries `h_m` so a non-isothermal model (e.g., NRLMSISE-00's
 * thermospheric T(altitude)) can be substituted later without
 * changing call sites.
 *
 * Reference: ICAO 1993 Standard Atmosphere (Document 7488/3) §2.5;
 * Liepmann & Roshko, Elements of Gasdynamics (1957), §1.4. */
#include "k26astro_atmos/atmos.h"

#include <math.h>

double k26astro_atmos_sound_speed_at(const K26AstroAtmos *a, double h_m)
{
    if (!a) return 0.0;
    K26AstroAtmosParams p = k26astro_atmos_params(a);
    if (h_m > p.atmos_top_m) return 0.0;
    double T     = p.t_0_k;
    double gamma = p.gamma_ratio_specific_heats;
    double R     = p.r_specific_j_per_kg_k;
    if (!(T > 0.0 && gamma > 0.0 && R > 0.0)) return 0.0;
    return sqrt(gamma * R * T);
}
