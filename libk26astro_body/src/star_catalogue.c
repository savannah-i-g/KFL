/* star_catalogue.c — K26AstroCatalogue (mmap'd binary loader) +
 * proper-motion + colour + distance utilities.
 *
 * File format defined in star.h. The loader mmaps the file
 * read-only; closing the catalogue is a single munmap. The baked-in
 * Hipparcos subset (star_catalogue_hip.c) wraps a const array in
 * the same opaque type so callers don't branch on "is this file-
 * backed or in-memory?".
 */
#include "k26astro_body/star.h"
#include "k26astro_core/consts.h"
#include "internal.h"

#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

/* ---- File header layout ---------------------------------- */

#define K26_STARS_HEADER_BYTES 64
#define K26_STARS_RECORD_BYTES sizeof(K26AstroStar)

/* ---- Loader ---------------------------------------------- */

K26AstroCatalogue *k26astro_catalogue_load(const char *path)
{
    if (!path) return NULL;
    int fd = open(path, O_RDONLY);
    if (fd < 0) return NULL;
    struct stat st;
    if (fstat(fd, &st) != 0) { close(fd); return NULL; }
    if (st.st_size < K26_STARS_HEADER_BYTES) { close(fd); return NULL; }

    void *map = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED) return NULL;

    const uint8_t *base = (const uint8_t *)map;
    if (memcmp(base, K26_STARS_MAGIC, K26_STARS_MAGIC_LEN) != 0) {
        munmap(map, (size_t)st.st_size); return NULL;
    }
    uint32_t version, endian, n_stars, rec_size;
    double   epoch_jd;
    uint64_t cat_id;
    memcpy(&version,  base + 8,  4);
    memcpy(&endian,   base + 12, 4);
    memcpy(&epoch_jd, base + 16, 8);
    memcpy(&n_stars,  base + 24, 4);
    memcpy(&rec_size, base + 28, 4);
    memcpy(&cat_id,   base + 32, 8);

    if (version != K26_STARS_FORMAT_VERSION) {
        munmap(map, (size_t)st.st_size); return NULL;
    }
    if (endian != K26_STARS_ENDIAN_MARKER) {
        munmap(map, (size_t)st.st_size); return NULL;
    }
    if (rec_size != (uint32_t)K26_STARS_RECORD_BYTES) {
        munmap(map, (size_t)st.st_size); return NULL;
    }
    size_t expected = K26_STARS_HEADER_BYTES + (size_t)n_stars * rec_size;
    if ((size_t)st.st_size < expected) {
        munmap(map, (size_t)st.st_size); return NULL;
    }

    K26AstroCatalogue *c = (K26AstroCatalogue *)calloc(1, sizeof(*c));
    if (!c) { munmap(map, (size_t)st.st_size); return NULL; }
    c->map          = map;
    c->map_size     = (size_t)st.st_size;
    c->records      = (const K26AstroStar *)(base + K26_STARS_HEADER_BYTES);
    c->n_records    = n_stars;
    c->epoch_jd_tt  = epoch_jd;
    c->catalogue_id = cat_id;
    c->file_backed  = 1;
    return c;
}

K26AstroCatalogue *k26astro_catalogue_load_default(void)
{
    K26AstroCatalogue *c =
        k26astro_catalogue_load("/usr/share/k26astro/stars/hip_naked.k26stars");
    if (c) return c;
    return k26astro_catalogue_hip_default();
}

void k26astro_catalogue_close(K26AstroCatalogue *c)
{
    if (!c) return;
    if (c->file_backed && c->map) munmap(c->map, c->map_size);
    free(c);
}

size_t k26astro_catalogue_count(const K26AstroCatalogue *c)
{
    return c ? c->n_records : 0;
}

const K26AstroStar *
k26astro_catalogue_at(const K26AstroCatalogue *c, size_t idx)
{
    if (!c || idx >= c->n_records) return NULL;
    return &c->records[idx];
}

double k26astro_catalogue_epoch_jd_tt(const K26AstroCatalogue *c)
{
    return c ? c->epoch_jd_tt : 0.0;
}

uint64_t k26astro_catalogue_id(const K26AstroCatalogue *c)
{
    return c ? c->catalogue_id : 0;
}

void k26astro_catalogue_each(const K26AstroCatalogue *c,
                              K26AstroStarFn fn, void *user)
{
    if (!c || !fn) return;
    for (size_t i = 0; i < c->n_records; i++) {
        fn(&c->records[i], user);
    }
}

const K26AstroStar *
k26astro_catalogue_find(const K26AstroCatalogue *c, const char *designation)
{
    if (!c || !designation) return NULL;
    for (size_t i = 0; i < c->n_records; i++) {
        if (strncasecmp(c->records[i].designation, designation,
                        K26ASTRO_STAR_DESIGNATION_MAX) == 0) {
            return &c->records[i];
        }
    }
    return NULL;
}

/* ---- Apparent direction --------------------------------- */

K26V3 k26astro_star_apparent_dir(const K26AstroStar *s,
                                  const K26AstroEpoch *t,
                                  double catalogue_epoch_jd_tt)
{
    K26V3 zero = { 0.0, 0.0, 0.0 };
    if (!s || !t) return zero;

    /* Years since catalogue epoch (in TT). */
    K26AstroEpoch tt = *t;
    if (tt.scale != K26A_TS_TT) k26astro_epoch_convert(&tt, K26A_TS_TT);
    double t_jd_tt = (double)tt.days_since_J2000 + tt.seconds_of_day / 86400.0
                   + 2451545.0;
    double dt_yr = (t_jd_tt - catalogue_epoch_jd_tt) / 365.25;

    /* Propagate RA, Dec via proper motion. Convert mas/yr → rad/yr. */
    const double mas_to_rad = 1.0e-3 * K26A_RAD_PER_ARCSEC;
    /* pm_ra_cos_dec is the RA proper motion already × cos(Dec) per
     * catalogue convention, so we divide back out to get dRA/dt. */
    double cos_dec = cos(s->dec);
    if (fabs(cos_dec) < 1.0e-12) cos_dec = 1.0e-12;
    double d_ra  = (s->pm_ra_cos_dec / cos_dec) * mas_to_rad * dt_yr;
    double d_dec =  s->pm_dec               * mas_to_rad * dt_yr;

    double ra_app  = s->ra  + d_ra;
    double dec_app = s->dec + d_dec;

    /* Build unit direction vector (ICRF). */
    double cd = cos(dec_app), sd = sin(dec_app);
    double cr = cos(ra_app),  sr = sin(ra_app);
    K26V3 out;
    out.x = cd * cr;
    out.y = cd * sr;
    out.z = sd;
    return out;
}

double k26astro_star_distance_m(const K26AstroStar *s)
{
    if (!s || s->parallax <= 0.0) return HUGE_VAL;
    /* parallax in mas → parsecs: d = 1000 / parallax_mas. */
    double d_pc = 1000.0 / s->parallax;
    return d_pc * K26A_PC_M;
}

void k26astro_star_color_srgb(const K26AstroStar *s, double out_rgb[3])
{
    if (!out_rgb) return;
    out_rgb[0] = out_rgb[1] = out_rgb[2] = 1.0;
    if (!s) return;
    /* Approximate B-V → temperature, then temperature → CIE → sRGB.
     * For starfield render we don't need photometric precision; a
     * piecewise polynomial mapping suffices.
     *
     * Effective temperature from B-V (Ballesteros 2012):
     *   T = 4600 K * (1 / (0.92 BV + 1.7) + 1 / (0.92 BV + 0.62))
     */
    double bv = s->b_v;
    double T = 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62));
    if (T < 1000.0)   T = 1000.0;
    if (T > 40000.0)  T = 40000.0;

    /* Approximate blackbody → sRGB via Krystek 1985 piecewise fits.
     * Output non-premultiplied; caller applies HDR exposure. */
    double r, g, b;
    /* Red */
    if (T <= 6600.0) {
        r = 1.0;
    } else {
        double x = (T - 6000.0) * 0.01;
        r = 329.698727446 * pow(x, -0.1332047592) / 255.0;
        if (r > 1.0) r = 1.0;
        if (r < 0.0) r = 0.0;
    }
    /* Green */
    if (T <= 6600.0) {
        double x = T * 0.01;
        g = (99.4708025861 * log(x) - 161.1195681661) / 255.0;
    } else {
        double x = (T - 6000.0) * 0.01;
        g = 288.1221695283 * pow(x, -0.0755148492) / 255.0;
    }
    if (g > 1.0) g = 1.0;
    if (g < 0.0) g = 0.0;
    /* Blue */
    if (T >= 6600.0) {
        b = 1.0;
    } else if (T <= 1900.0) {
        b = 0.0;
    } else {
        double x = (T - 1000.0) * 0.01;
        b = (138.5177312231 * log(x) - 305.0447927307) / 255.0;
        if (b > 1.0) b = 1.0;
        if (b < 0.0) b = 0.0;
    }
    out_rgb[0] = r;
    out_rgb[1] = g;
    out_rgb[2] = b;
}
