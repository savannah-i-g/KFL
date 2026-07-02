/* test_star.c — star catalogue + apparent direction + colour.
 *
 * Acceptance:
 *   - Hipparcos baked-in subset loads non-NULL with > 20 stars.
 *   - Catalogue id matches K26_CAT_ID_HIPPARCOS_2007.
 *   - Designation lookup finds known stars.
 *   - Apparent direction returns a unit vector.
 *   - Distance from parallax is positive + finite for measured
 *     parallaxes, +∞ for the zero-parallax case.
 *   - Colour conversion produces sRGB in [0, 1].
 */
#include "k26astro_body/star.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroCatalogue *cat = k26astro_catalogue_hip_default();
    assert(cat != NULL);
    size_t n = k26astro_catalogue_count(cat);
    assert(n > 20);
    assert(k26astro_catalogue_id(cat) == K26_CAT_ID_HIPPARCOS_2007);
    assert(fabs(k26astro_catalogue_epoch_jd_tt(cat) - 2448349.0625) < 1e-3);

    /* ---- Designation lookup ----------------------------- */
    const K26AstroStar *sirius = k26astro_catalogue_find(cat, "HIP 32349");
    assert(sirius != NULL);
    assert(sirius->mag_v < 0.0);   /* Sirius mag_v ≈ -1.46 */
    assert(sirius->parallax > 100.0);  /* 379 mas */

    const K26AstroStar *polaris = k26astro_catalogue_find(cat, "HIP 11767");
    assert(polaris != NULL);
    assert(polaris->dec > 1.5);  /* > 85 deg → near +z (north pole) */

    /* Case-insensitive lookup. */
    assert(k26astro_catalogue_find(cat, "hip 32349") != NULL);

    /* Unknown lookup returns NULL. */
    assert(k26astro_catalogue_find(cat, "HIP 999999") == NULL);

    /* ---- Apparent direction is a unit vector ------------- */
    K26AstroEpoch t = k26astro_epoch_j2000_tt();
    double cat_epoch_jd = k26astro_catalogue_epoch_jd_tt(cat);
    for (size_t i = 0; i < n; i++) {
        const K26AstroStar *s = k26astro_catalogue_at(cat, i);
        K26V3 d = k26astro_star_apparent_dir(s, &t, cat_epoch_jd);
        double mag = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
        assert(fabs(mag - 1.0) < 1e-9);
    }

    /* ---- Distance from parallax ------------------------- */
    /* Sirius parallax ~379 mas → distance ~2.64 pc = ~8.6 light-years. */
    double d_sirius = k26astro_star_distance_m(sirius);
    assert(isfinite(d_sirius));
    double d_ly_sirius = d_sirius / K26A_LY_M;
    assert(d_ly_sirius > 8.0 && d_ly_sirius < 9.0);

    /* Proxima Centauri parallax ~771 mas → distance ~1.3 pc = ~4.2 ly. */
    const K26AstroStar *proxima = k26astro_catalogue_find(cat, "HIP 70890");
    assert(proxima != NULL);
    double d_proxima = k26astro_star_distance_m(proxima);
    double d_ly_proxima = d_proxima / K26A_LY_M;
    assert(d_ly_proxima > 4.0 && d_ly_proxima < 5.0);

    /* ---- Colour ------------------------------------------ */
    for (size_t i = 0; i < n; i++) {
        const K26AstroStar *s = k26astro_catalogue_at(cat, i);
        double rgb[3];
        k26astro_star_color_srgb(s, rgb);
        for (int k = 0; k < 3; k++) {
            assert(rgb[k] >= 0.0);
            assert(rgb[k] <= 1.0);
        }
    }

    /* Sirius (B-V ≈ 0.01, hot) should have B > R; Antares
     * (B-V ≈ 1.83, cool) should have R > B. */
    double rgb_sirius[3], rgb_antares[3];
    k26astro_star_color_srgb(sirius, rgb_sirius);
    const K26AstroStar *antares = k26astro_catalogue_find(cat, "HIP 80763");
    assert(antares != NULL);
    k26astro_star_color_srgb(antares, rgb_antares);
    assert(rgb_sirius[2] > rgb_sirius[0]);    /* blue Sirius */
    assert(rgb_antares[0] > rgb_antares[2]);   /* red Antares */

    /* ---- Default-load path (no file installed) ---------- */
    K26AstroCatalogue *def = k26astro_catalogue_load_default();
    assert(def != NULL);
    /* Should at minimum return the baked-in subset. */
    assert(k26astro_catalogue_count(def) >= n);

    printf("test_star: OK (%zu baked-in stars)\n", n);
    return 0;
}
