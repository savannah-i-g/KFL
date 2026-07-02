/* OSM tag value parsers + the highway-class width table.
 *
 * Rules from the OSM wiki:
 *   - Key:height                — metres by default; "10 m", "10.5", "0,8", "7'4\"".
 *   - Key:building:levels       — positive integer count of above-ground floors.
 *   - Key:building:min_level    — same shape, defaults to 0.
 *   - Key:width (on highway=*)  — metres.
 *   - Key:lanes                 — count; no canonical metres/lane (we use 3.5).
 *   - Key:highway               — class enum; widths set by the renderer.
 *
 * The class-default width table below is our own. Sources cross-checked:
 *   - OSM Carto default stylesheet (motorway/trunk on the wider side).
 *   - UK DfT manual for urban arterials (~14 m motorway, ~9 m primary).
 *   - Footway/path 1.5 m matches OSM "narrow" convention.
 *
 * These are *visual* widths for a 3-D renderer at city scale; not survey
 * data. If reality calls for 28 m motorways the user tags `width=*` and
 * we honour it.
 */

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

double k26geo_default_building_height_m(void)
{
    return K26GEO_DEFAULT_BUILDING_HEIGHT_M;
}
double k26geo_metres_per_level(void)  { return K26GEO_METRES_PER_LEVEL;  }
double k26geo_metres_per_lane (void)  { return K26GEO_METRES_PER_LANE;   }

/* ---- string helpers --------------------------------------------- */

static const char *skip_ws(const char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

/* Strip leading/trailing whitespace by writing into `buf` (size cap).
 * Returns the new in-place length. */
static size_t strip_into(char *buf, size_t cap, const char *s)
{
    if (!buf || cap == 0) return 0;
    s = skip_ws(s);
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) n--;
    if (n >= cap) n = cap - 1;
    memcpy(buf, s, n);
    buf[n] = '\0';
    return n;
}

/* Replace any commas with periods (continental decimal separator). */
static void normalise_decimal(char *buf)
{
    for (char *p = buf; *p; p++) if (*p == ',') *p = '.';
}

/* ---- height parser ---------------------------------------------- *
 *
 * Decision table:
 *   1. Strip whitespace; lowercase tail.
 *   2. If contains "'", treat as feet-and-inches.
 *   3. Else strtod the numeric prefix, then inspect the tail unit.
 *      - "" / "m"        → metres
 *      - "ft"            → feet  → metres
 *      - others          → reject
 *   4. Reject finite-zero / negative / non-finite.
 */
static K26GeoStatus parse_feet_inches(const char *s, double *out_m)
{
    /* Accept patterns: 7'4", 7' 4", 7' (no inches), 7'4 (no closing quote) */
    char *endp = NULL;
    const double ft = strtod(s, &endp);
    if (endp == s) return K26GEO_ERR_PARSE;
    while (*endp && isspace((unsigned char)*endp)) endp++;
    if (*endp != '\'') return K26GEO_ERR_PARSE;
    endp++;
    while (*endp && isspace((unsigned char)*endp)) endp++;

    double in = 0.0;
    if (*endp && *endp != '"') {
        char *endp2 = NULL;
        in = strtod(endp, &endp2);
        if (endp2 == endp) return K26GEO_ERR_PARSE;
        endp = endp2;
        while (*endp && isspace((unsigned char)*endp)) endp++;
    }
    if (*endp == '"') endp++;
    while (*endp && isspace((unsigned char)*endp)) endp++;
    if (*endp) return K26GEO_ERR_PARSE;

    if (!isfinite(ft) || !isfinite(in)) return K26GEO_ERR_PARSE;
    if (ft < 0.0 || in < 0.0)           return K26GEO_ERR_PARSE;

    const double m = ft * 0.3048 + in * 0.0254;
    if (m <= 0.0)                       return K26GEO_ERR_PARSE;
    *out_m = m;
    return K26GEO_OK;
}

K26GeoStatus k26geo_parse_height_m(const char *s, double *out_m)
{
    if (!s || !out_m) return K26GEO_ERR_INVAL;

    char buf[64];
    size_t n = strip_into(buf, sizeof buf, s);
    if (n == 0) return K26GEO_ERR_PARSE;
    normalise_decimal(buf);

    /* Reject explicit range / inequality forms. */
    for (size_t i = 0; i < n; i++) {
        if (buf[i] == '-' && i > 0 && (isdigit((unsigned char)buf[i - 1]) ||
                                       buf[i - 1] == '.'))
            return K26GEO_ERR_PARSE;
        if (buf[i] == '<' || buf[i] == '>' || buf[i] == '~')
            return K26GEO_ERR_PARSE;
    }

    if (strchr(buf, '\'')) return parse_feet_inches(buf, out_m);

    char *endp = NULL;
    double v = strtod(buf, &endp);
    if (endp == buf || !isfinite(v)) return K26GEO_ERR_PARSE;
    if (v <= 0.0)                     return K26GEO_ERR_PARSE;

    /* Optional unit tail. Tolerate `3m` (no space) and `3 m`. */
    while (*endp && isspace((unsigned char)*endp)) endp++;
    if (*endp == '\0') { *out_m = v; return K26GEO_OK; }

    char tail[8] = {0};
    size_t ti = 0;
    while (*endp && ti < sizeof tail - 1) {
        tail[ti++] = (char)tolower((unsigned char)*endp++);
    }
    if (*endp) return K26GEO_ERR_PARSE;

    if (strcmp(tail, "m") == 0)   { *out_m = v;          return K26GEO_OK; }
    if (strcmp(tail, "ft") == 0)  { *out_m = v * 0.3048; return K26GEO_OK; }
    return K26GEO_ERR_PARSE;
}

K26GeoStatus k26geo_parse_width_m(const char *s, double *out_m)
{
    /* Same shape as height but no foot-inches form expected (rarely seen
     * for road widths). Reuse the parser; if you tag width=7'4" it just
     * works. */
    return k26geo_parse_height_m(s, out_m);
}

K26GeoStatus k26geo_parse_levels_m(const char *s, double *out_m)
{
    if (!s || !out_m) return K26GEO_ERR_INVAL;

    char buf[32];
    size_t n = strip_into(buf, sizeof buf, s);
    if (n == 0) return K26GEO_ERR_PARSE;

    char *endp = NULL;
    long v = strtol(buf, &endp, 10);
    if (endp == buf) return K26GEO_ERR_PARSE;
    while (*endp && isspace((unsigned char)*endp)) endp++;
    if (*endp) return K26GEO_ERR_PARSE;
    if (v <= 0 || v > 200) return K26GEO_ERR_PARSE;   /* 200 floors is the cap */

    *out_m = (double)v * K26GEO_METRES_PER_LEVEL;
    return K26GEO_OK;
}

/* ---- highway class width table --------------------------------- */

struct ClassW { const char *k; double w; };

static const struct ClassW HW[] = {
    { "motorway",          14.0 },
    { "motorway_link",      8.0 },
    { "trunk",             12.0 },
    { "trunk_link",         7.0 },
    { "primary",            9.0 },
    { "primary_link",       6.0 },
    { "secondary",          7.0 },
    { "secondary_link",     5.0 },
    { "tertiary",           6.0 },
    { "tertiary_link",      5.0 },
    { "unclassified",       5.0 },
    { "residential",        5.0 },
    { "living_street",      4.0 },
    { "service",            4.0 },
    { "pedestrian",         5.0 },
    { "track",              3.0 },
    { "bus_guideway",       3.5 },
    { "raceway",            8.0 },
    { "road",               5.0 },
    { "footway",            1.5 },
    { "path",               1.5 },
    { "cycleway",           2.0 },
    { "bridleway",          1.8 },
    { "steps",              1.2 },
    { "corridor",           1.5 },
    { NULL,                 0.0 }
};

double k26geo_highway_class_width_m(const char *value)
{
    if (!value) return 0.0;
    for (size_t i = 0; HW[i].k; i++) {
        if (strcmp(HW[i].k, value) == 0) return HW[i].w;
    }
    return 0.0;
}

double k26geo_resolve_road_width_m(const char *highway_value,
                                   K26GeoTagLookupFn lookup, void *user)
{
    /* 1. Explicit width=*. */
    if (lookup) {
        const char *w_str = lookup("width", user);
        if (w_str) {
            double w = 0.0;
            if (k26geo_parse_width_m(w_str, &w) == K26GEO_OK && w > 0.0) return w;
        }
        /* 2. lanes=* fallback. */
        const char *l_str = lookup("lanes", user);
        if (l_str) {
            char *endp = NULL;
            long lanes = strtol(l_str, &endp, 10);
            if (endp != l_str && lanes > 0 && lanes < 32) {
                return (double)lanes * K26GEO_METRES_PER_LANE;
            }
        }
    }
    /* 3. Class default. */
    const double cls = k26geo_highway_class_width_m(highway_value);
    if (cls > 0.0) return cls;
    return 5.0;   /* "road" — most generous fallback */
}
