/* refraction.c — atmospheric refraction (Bennett 1982 closed-form).
 *
 * Reference: G. G. Bennett, "The calculation of astronomical
 * refraction in marine navigation", Journal of Navigation 35:255
 * (1982). The closed-form for refraction angle R (arcminutes) given
 * true elevation h (degrees) is
 *
 *     R(h) = 1 / tan(h + 7.31 / (h + 4.4))    arcminutes
 *
 * Valid for h in [0°, 90°]. R → 34.5' at h=0 (horizon), R → 0 at
 * h=90° (zenith). The formula is asymptotic at the horizon, so we
 * clamp to the Bennett value at h=0; below-horizon targets use the
 * same refraction until the apparent direction lifts above the
 * geometric horizon.
 *
 * Non-Earth atmospheres: scale the Bennett refraction by the
 * ratio of n₀_atmos / n₀_earth. For Venus (n₀ ≈ 6e-3) this gives
 * ~20× the Earth refraction, qualitatively correct for the dense
 * Venusian atmosphere. Ray-traced Snell integration over an
 * NRLMSISE-class atmosphere is not currently provided. */
#include "k26astro_atmos/atmos.h"
#include "k26astro_atmos/atmos_consts.h"

#include <math.h>

static double bennett_arcmin_(double elev_deg)
{
    if (elev_deg < 0.0)   elev_deg = 0.0;
    if (elev_deg >= 90.0) return 0.0;
    double arg_deg = elev_deg + 7.31 / (elev_deg + 4.4);
    double arg_rad = arg_deg * (M_PI / 180.0);
    double t = tan(arg_rad);
    if (!(t > 0.0)) return 34.5;
    return 1.0 / t;
}

double k26astro_atmos_refraction_rad(const K26AstroAtmos *a,
                                      double true_elevation_rad)
{
    if (!a) return 0.0;
    double elev_deg = true_elevation_rad * (180.0 / M_PI);
    double R_arcmin = bennett_arcmin_(elev_deg);
    K26AstroAtmosParams p = k26astro_atmos_params(a);
    double scale = (p.n0 > 0.0)
                 ? p.n0 / K26ASTRO_ATMOS_EARTH_N0
                 : 1.0;
    return R_arcmin * (M_PI / (180.0 * 60.0)) * scale;
}

K26V3 k26astro_atmos_apparent(const K26AstroAtmos *a,
                               K26V3 geometric_dir,
                               K26V3 zenith_dir)
{
    /* Normalise inputs. */
    double gn = sqrt(geometric_dir.x*geometric_dir.x
                   + geometric_dir.y*geometric_dir.y
                   + geometric_dir.z*geometric_dir.z);
    double zn = sqrt(zenith_dir.x*zenith_dir.x
                   + zenith_dir.y*zenith_dir.y
                   + zenith_dir.z*zenith_dir.z);
    if (!(gn > 0.0) || !(zn > 0.0) || !a) return geometric_dir;

    K26V3 g = { geometric_dir.x / gn,
                geometric_dir.y / gn,
                geometric_dir.z / gn };
    K26V3 z = { zenith_dir.x / zn,
                zenith_dir.y / zn,
                zenith_dir.z / zn };

    /* True elevation: arcsin of dot product. */
    double sin_elev = g.x*z.x + g.y*z.y + g.z*z.z;
    if (sin_elev >  1.0) sin_elev =  1.0;
    if (sin_elev < -1.0) sin_elev = -1.0;
    double elev = asin(sin_elev);

    /* At/above zenith: no refraction. */
    if (elev >= 0.5 * M_PI - 1.0e-12) return g;

    /* Bennett refraction in radians. */
    double R = k26astro_atmos_refraction_rad(a, elev);
    if (R <= 0.0) return g;

    /* Apparent elevation = true + R, clamped to apparent horizon. */
    double apparent = elev + R;
    if (apparent >= 0.5 * M_PI) apparent = 0.5 * M_PI;

    /* Rotate g toward z by (apparent - elev) radians in the plane
     * spanned by g and z. Build perpendicular component p (in plane,
     * perpendicular to g) pointing in the direction of z's
     * projection onto the perpendicular subspace:
     *   p_raw = z - (z·g) g
     *   p = p_raw / |p_raw|
     * Then apparent_dir = g cos(R) + p sin(R). */
    double zg = z.x*g.x + z.y*g.y + z.z*g.z;
    K26V3 p_raw = { z.x - zg * g.x,
                    z.y - zg * g.y,
                    z.z - zg * g.z };
    double pn = sqrt(p_raw.x*p_raw.x + p_raw.y*p_raw.y + p_raw.z*p_raw.z);
    if (!(pn > 0.0)) return g;  /* g parallel to z — at zenith */
    K26V3 p = { p_raw.x / pn, p_raw.y / pn, p_raw.z / pn };

    double cR = cos(R);
    double sR = sin(R);
    K26V3 out = {
        g.x * cR + p.x * sR,
        g.y * cR + p.y * sR,
        g.z * cR + p.z * sR
    };
    /* Renormalise (against numerical drift). */
    double on = sqrt(out.x*out.x + out.y*out.y + out.z*out.z);
    if (!(on > 0.0)) return g;
    out.x /= on; out.y /= on; out.z /= on;
    return out;
}
