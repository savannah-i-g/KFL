/* k26astro_body/star.h — astrometric star records + catalogue API.
 *
 * Why stars aren't K26AstroBody: a typical star catalogue has 10⁴-
 * 10⁹ entries; integrating each one is pointless (stars don't
 * perturb each other on human timescales). Stars are queried en
 * masse per frame for rendering; the natural representation isn't
 * Cartesian state but astrometric coordinates (RA, Dec, parallax,
 * proper motion) at a catalogue epoch.
 *
 * K26 ships:
 *   - Hipparcos naked-eye subset (~1 MB, baked into the ISO via
 *     star_catalogue_hip.c — provides a default starfield without
 *     external files)
 *   - Hipparcos full (~5 MB) as apk add-on
 *   - Gaia DR3 bright subset (~120 MB) as apk add-on
 *   - Gaia DR4 (post-2026-12-02 release) as future apk add-on
 *
 * Catalogue file format /usr/share/k26astro/stars/[name].k26stars
 * (mmap'd, little-endian):
 *
 *   header:
 *     bytes 0-7  : magic "K26STARS"
 *     bytes 8-11 : version (uint32)
 *     bytes 12-15: endian marker (uint32 = 0x01020304)
 *     bytes 16-23: epoch_jd (TT scale) (double)
 *     bytes 24-27: n_stars (uint32)
 *     bytes 28-31: record_size in bytes (uint32)
 *     bytes 32-39: catalogue_id (uint64) — K26_CAT_ID_*
 *     bytes 40-63: reserved
 *   records:
 *     n_stars × record_size bytes; layout matches K26AstroStarRecord
 *
 * Catalogue ID values let readers identify which catalogue produced
 * the file without ambiguity (Hipparcos vs van Leeuwen 2007 vs
 * Gaia DR3 vs Gaia DR4); each new release gets a new ID. */
#ifndef K26ASTRO_BODY_STAR_H
#define K26ASTRO_BODY_STAR_H

#include <stddef.h>
#include <stdint.h>

#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Catalogue file format constants ----------------------- */
#define K26_STARS_MAGIC          "K26STARS"
#define K26_STARS_MAGIC_LEN      8
#define K26_STARS_FORMAT_VERSION 1u
#define K26_STARS_ENDIAN_MARKER  0x01020304u

/* Catalogue identifiers — new release → new ID; never reuse. */
#define K26_CAT_ID_HIPPARCOS_1997   1u   /* original Hipparcos 1997 */
#define K26_CAT_ID_HIPPARCOS_2007   2u   /* van Leeuwen 2007 New Reduction */
#define K26_CAT_ID_TYCHO_2          3u
#define K26_CAT_ID_GAIA_DR3         4u
#define K26_CAT_ID_GAIA_DR4         5u   /* post-2026-12-02 release */

/* ---- Star record ----------------------------------------- */

#define K26ASTRO_STAR_DESIGNATION_MAX 32

typedef struct {
    /* identity */
    uint64_t   id;
    char       designation[K26ASTRO_STAR_DESIGNATION_MAX];  /* HIP/Tycho/Gaia id */

    /* astrometric data, at the catalogue's stored epoch */
    double     ra;                /* radians at catalogue epoch */
    double     dec;               /* radians */
    double     parallax;          /* mas */
    double     pm_ra_cos_dec;     /* mas/yr */
    double     pm_dec;            /* mas/yr */
    double     rv;                /* km/s, radial velocity */

    /* photometry */
    double     mag_v;             /* apparent V-band magnitude */
    double     b_v;               /* B - V colour index */
} K26AstroStar;

/* ---- Star operations ----------------------------------------- */

/* Apparent direction (ICRF unit vector) at simulation epoch t.
 * Propagates RA/Dec from catalogue epoch via proper motion. */
K26V3  k26astro_star_apparent_dir(const K26AstroStar *s,
                                   const K26AstroEpoch *t,
                                   double catalogue_epoch_jd_tt);

/* Distance in metres, computed from parallax. Returns +infinity
 * (HUGE_VAL) when parallax ≤ 0 (catalogue records with unmeasured
 * parallax). */
double k26astro_star_distance_m(const K26AstroStar *s);

/* B-V → CIE chromaticity → display sRGB approximation. Used by the
 * starfield renderer to give cool stars a redder tint, hot stars a
 * bluer one. Output is non-premultiplied sRGB in [0, 1] (caller
 * applies HDR exposure curve). */
void   k26astro_star_color_srgb(const K26AstroStar *s,
                                double out_rgb[3]);

/* ---- Catalogue opaque ---------------------------------------- */
typedef struct K26AstroCatalogue K26AstroCatalogue;

K26AstroCatalogue *k26astro_catalogue_load(const char *path);
K26AstroCatalogue *k26astro_catalogue_load_default(void);
                  /* loads /usr/share/k26astro/stars/hip_naked.k26stars
                   * if installed; falls back to the baked-in subset
                   * via k26astro_catalogue_hip_default() */
void               k26astro_catalogue_close(K26AstroCatalogue *c);

size_t  k26astro_catalogue_count(const K26AstroCatalogue *c);
const K26AstroStar *
        k26astro_catalogue_at(const K26AstroCatalogue *c, size_t idx);

/* Catalogue's epoch — Julian Date (TT). Hipparcos uses 1991.25;
 * van Leeuwen 2007 uses 1991.25; Gaia DR3 uses 2016.0; Gaia DR4
 * uses 2017.5. */
double  k26astro_catalogue_epoch_jd_tt(const K26AstroCatalogue *c);

/* Catalogue id — see K26_CAT_ID_* above. */
uint64_t k26astro_catalogue_id(const K26AstroCatalogue *c);

/* Iteration. */
typedef void (*K26AstroStarFn)(const K26AstroStar *s, void *user);
void k26astro_catalogue_each(const K26AstroCatalogue *c,
                              K26AstroStarFn fn, void *user);

/* Lookup by designation (HIP/Tycho/Gaia id). Returns NULL if not
 * found. */
const K26AstroStar *
     k26astro_catalogue_find(const K26AstroCatalogue *c,
                              const char *designation);

/* ---- Hipparcos baked-in subset --------------------------- *
 *
 * Returns a fixed catalogue containing the major naked-eye stars
 * baked into the binary. Used by k26astro_catalogue_load_default
 * when the default file isn't installed. */
K26AstroCatalogue *k26astro_catalogue_hip_default(void);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_BODY_STAR_H */
