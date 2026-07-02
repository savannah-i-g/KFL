/* .k26atlas — bespoke text session-file reader + emitter.
 *
 * Format spec lives in include/k26geo.h above the K26AtlasSession
 * declaration. Hand-rolled lexer; indentation is decorative. Comments
 * begin with `#` and run to end-of-line. Strings use C-style double
 * quotes (no escapes for v1; the only string fields are title/source
 * paths and we forbid double-quotes in them).
 *
 * Numbers: doubles parsed with strtod. Booleans: `on`/`off`/`true`/
 * `false` (case-insensitive). Blocks (`camera`, `layers`, `route`,
 * `annotation`) end with the matching `end` keyword. */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

void k26atlas_session_defaults(K26AtlasSession *s)
{
    if (!s) return;
    memset(s, 0, sizeof *s);
    snprintf(s->title, sizeof s->title, "Untitled");
    s->cam_eye    = (K26V3){ 0.0, 200.0,  400.0 };
    s->cam_target = (K26V3){ 0.0,   0.0,    0.0 };
    s->cam_up     = (K26V3){ 0.0,   1.0,    0.0 };
    s->cam_fov_y_deg     = 55.0;
    s->layer_mask        = K26ATLAS_DEFAULT_LAYER_MASK;
    s->route_origin_lon  = (double)NAN;
    s->route_origin_lat  = (double)NAN;
    s->route_dest_lon    = (double)NAN;
    s->route_dest_lat    = (double)NAN;
}

/* ---- writer ----------------------------------------------------- */

static int write_str_field(FILE *f, const char *indent,
                           const char *key, const char *value)
{
    /* Strings may not contain double quotes in v1. Replace them with
     * a placeholder rather than fail; preserves round-trip on the
     * common case. */
    if (fprintf(f, "%s%s \"", indent, key) < 0) return -1;
    for (const char *p = value ? value : ""; *p; p++) {
        char c = (*p == '"' ? '\'' : *p);
        if (fputc(c, f) == EOF) return -1;
    }
    if (fputs("\"\n", f) == EOF) return -1;
    return 0;
}

K26GeoStatus k26atlas_session_save(const char *path,
                                   const K26AtlasSession *s,
                                   K26AtlasAnnCountFn count_fn,
                                   K26AtlasAnnEmitFn  emit_fn,
                                   void *user)
{
    if (!path || !s) return K26GEO_ERR_INVAL;

    FILE *f = fopen(path, "wb");
    if (!f) return K26GEO_ERR_IO;

    if (fprintf(f, "atlas v%d\n", K26ATLAS_FORMAT_VERSION) < 0)
        goto io_err;

    if (write_str_field(f, "    ", "title",  s->title))            goto io_err;
    if (write_str_field(f, "    ", "source", s->source_path))      goto io_err;
    if (s->scenario_path[0])
        if (write_str_field(f, "    ", "scenario_path",
                            s->scenario_path))                     goto io_err;
    if (fprintf(f, "    origin %.9f %.9f\n",
                s->origin_lon, s->origin_lat) < 0)                 goto io_err;

    /* Camera block. */
    if (fprintf(f, "\n    camera\n") < 0)                          goto io_err;
    if (fprintf(f, "        eye    %.6f %.6f %.6f\n",
                s->cam_eye.x, s->cam_eye.y, s->cam_eye.z) < 0)     goto io_err;
    if (fprintf(f, "        target %.6f %.6f %.6f\n",
                s->cam_target.x, s->cam_target.y, s->cam_target.z) < 0) goto io_err;
    if (fprintf(f, "        up     %.6f %.6f %.6f\n",
                s->cam_up.x, s->cam_up.y, s->cam_up.z) < 0)        goto io_err;
    if (fprintf(f, "        fov_y  %.4f\n", s->cam_fov_y_deg) < 0) goto io_err;
    if (fprintf(f, "    end\n") < 0)                               goto io_err;

    /* Layers block. */
    if (fprintf(f, "\n    layers\n") < 0)                          goto io_err;
    static const struct { const char *name; unsigned bit; } LAYERS[] = {
        { "buildings",  K26ATLAS_LAYER_BUILDINGS  },
        { "roads",      K26ATLAS_LAYER_ROADS      },
        { "water",      K26ATLAS_LAYER_WATER      },
        { "path",       K26ATLAS_LAYER_PATH       },
        { "railways",   K26ATLAS_LAYER_RAILWAYS   },
        { "boundaries", K26ATLAS_LAYER_BOUNDARIES },
        { "pois",       K26ATLAS_LAYER_POIS       },
        { "labels",     K26ATLAS_LAYER_LABELS     },
        { "route",      K26ATLAS_LAYER_ROUTE      },
        { NULL, 0 }
    };
    for (size_t i = 0; LAYERS[i].name; i++) {
        const char *on = (s->layer_mask & LAYERS[i].bit) ? "on" : "off";
        if (fprintf(f, "        %-11s %s\n", LAYERS[i].name, on) < 0) goto io_err;
    }
    if (fprintf(f, "    end\n") < 0)                               goto io_err;

    /* Route block. */
    if (isfinite(s->route_origin_lon) && isfinite(s->route_origin_lat) &&
        isfinite(s->route_dest_lon)   && isfinite(s->route_dest_lat)) {
        if (fprintf(f, "\n    route\n") < 0)                       goto io_err;
        if (fprintf(f, "        origin %.9f %.9f\n",
                    s->route_origin_lon, s->route_origin_lat) < 0) goto io_err;
        if (fprintf(f, "        dest   %.9f %.9f\n",
                    s->route_dest_lon, s->route_dest_lat) < 0)     goto io_err;
        if (fprintf(f, "    end\n") < 0)                           goto io_err;
    }

    /* Annotations. The library wraps each one with `annotation <kind>` /
     * `end`; the caller's emit_fn writes the body. */
    if (count_fn && emit_fn) {
        int n = count_fn(user);
        for (int i = 0; i < n; i++) {
            if (fprintf(f, "\n    annotation\n") < 0)              goto io_err;
            if (emit_fn(f, (size_t)i, user) != 0)                  goto io_err;
            if (fprintf(f, "    end\n") < 0)                       goto io_err;
        }
    }

    if (fprintf(f, "end\n") < 0)                                   goto io_err;
    if (fclose(f) != 0)        return K26GEO_ERR_IO;
    return K26GEO_OK;

io_err:
    fclose(f);
    return K26GEO_ERR_IO;
}

/* ---- reader ----------------------------------------------------- */

/* The reader is a simple line-oriented scanner. We strip `#` comments,
 * trim whitespace, then dispatch by the first token. Quoted strings
 * may span at most one line. */

static char *strip(char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = '\0';
    return s;
}

static void strip_comment(char *s)
{
    int in_q = 0;
    for (char *p = s; *p; p++) {
        if (*p == '"') in_q = !in_q;
        else if (*p == '#' && !in_q) { *p = '\0'; return; }
    }
}

/* Parse a quoted string starting at *src. On success, returns a
 * pointer into `out_buf`, advances *src past the closing quote, and
 * sets *out_len. On failure returns NULL. */
static const char *parse_quoted(const char **src, char *out_buf, size_t cap)
{
    const char *p = *src;
    while (*p && isspace((unsigned char)*p)) p++;
    if (*p != '"') return NULL;
    p++;
    size_t n = 0;
    while (*p && *p != '"' && n + 1 < cap) {
        out_buf[n++] = *p++;
    }
    if (*p != '"') return NULL;
    p++;
    out_buf[n] = '\0';
    *src = p;
    return out_buf;
}

static int parse_bool(const char *s, int *out)
{
    if (!strcasecmp(s, "on") || !strcasecmp(s, "true")  ||
        !strcasecmp(s, "yes")|| !strcmp(s, "1")) { *out = 1; return 1; }
    if (!strcasecmp(s, "off")|| !strcasecmp(s, "false") ||
        !strcasecmp(s, "no") || !strcmp(s, "0")) { *out = 0; return 1; }
    return 0;
}

K26GeoStatus k26atlas_session_load(const char *path,
                                   K26AtlasSession *s,
                                   K26AtlasAnnParseFn parse_fn,
                                   void *user)
{
    if (!path || !s) return K26GEO_ERR_INVAL;
    k26atlas_session_defaults(s);

    FILE *f = fopen(path, "rb");
    if (!f) return K26GEO_ERR_IO;

    char line[2048];
    int  saw_magic   = 0;
    enum { BLK_NONE, BLK_CAMERA, BLK_LAYERS, BLK_ROUTE, BLK_ANN } block = BLK_NONE;
    char ann_kind[64] = "";
    K26GeoStatus rc = K26GEO_OK;

    while (fgets(line, sizeof line, f)) {
        strip_comment(line);
        char *t = strip(line);
        if (*t == '\0') continue;

        /* Tokenise: take first whitespace-separated word as keyword. */
        char *kw = t;
        char *rest = t;
        while (*rest && !isspace((unsigned char)*rest)) rest++;
        if (*rest) { *rest++ = '\0'; while (*rest && isspace((unsigned char)*rest)) rest++; }

        if (!saw_magic) {
            /* Expect "atlas v1" */
            if (strcmp(kw, "atlas") != 0) { rc = K26GEO_ERR_PARSE; break; }
            /* rest is "v1" or "v<N>" */
            if (rest[0] != 'v') { rc = K26GEO_ERR_PARSE; break; }
            long v = strtol(rest + 1, NULL, 10);
            if (v != K26ATLAS_FORMAT_VERSION) { rc = K26GEO_ERR_UNSUPPORTED; break; }
            saw_magic = 1;
            continue;
        }

        if (block == BLK_ANN && parse_fn) {
            if (strcmp(kw, "end") == 0) {
                block = BLK_NONE;
                continue;
            }
            char body[1536];
            snprintf(body, sizeof body, "%s %s", kw, rest);
            if (parse_fn(ann_kind, body, user) != 0) { rc = K26GEO_ERR_PARSE; break; }
            continue;
        }

        if (block == BLK_CAMERA) {
            if (strcmp(kw, "end") == 0) { block = BLK_NONE; continue; }
            if (strcmp(kw, "eye")    == 0) sscanf(rest, "%lf %lf %lf",
                                                  &s->cam_eye.x, &s->cam_eye.y, &s->cam_eye.z);
            else if (strcmp(kw, "target") == 0) sscanf(rest, "%lf %lf %lf",
                                                  &s->cam_target.x, &s->cam_target.y, &s->cam_target.z);
            else if (strcmp(kw, "up") == 0) sscanf(rest, "%lf %lf %lf",
                                                  &s->cam_up.x, &s->cam_up.y, &s->cam_up.z);
            else if (strcmp(kw, "fov_y") == 0) s->cam_fov_y_deg = strtod(rest, NULL);
            continue;
        }

        if (block == BLK_LAYERS) {
            if (strcmp(kw, "end") == 0) { block = BLK_NONE; continue; }
            int on = 1; if (!parse_bool(rest, &on)) on = 1;
            unsigned bit = 0;
            if      (strcmp(kw, "buildings") == 0)  bit = K26ATLAS_LAYER_BUILDINGS;
            else if (strcmp(kw, "roads")     == 0)  bit = K26ATLAS_LAYER_ROADS;
            else if (strcmp(kw, "water")     == 0)  bit = K26ATLAS_LAYER_WATER;
            else if (strcmp(kw, "path")      == 0)  bit = K26ATLAS_LAYER_PATH;
            else if (strcmp(kw, "railways")  == 0)  bit = K26ATLAS_LAYER_RAILWAYS;
            else if (strcmp(kw, "boundaries")== 0)  bit = K26ATLAS_LAYER_BOUNDARIES;
            else if (strcmp(kw, "pois")      == 0)  bit = K26ATLAS_LAYER_POIS;
            else if (strcmp(kw, "labels")    == 0)  bit = K26ATLAS_LAYER_LABELS;
            else if (strcmp(kw, "route")     == 0)  bit = K26ATLAS_LAYER_ROUTE;
            if (bit) {
                if (on) s->layer_mask |=  bit;
                else    s->layer_mask &= ~bit;
            }
            continue;
        }

        if (block == BLK_ROUTE) {
            if (strcmp(kw, "end") == 0) { block = BLK_NONE; continue; }
            if (strcmp(kw, "origin") == 0)
                sscanf(rest, "%lf %lf", &s->route_origin_lon, &s->route_origin_lat);
            else if (strcmp(kw, "dest") == 0)
                sscanf(rest, "%lf %lf", &s->route_dest_lon, &s->route_dest_lat);
            continue;
        }

        /* Top-level statements / block openers. */
        if (strcmp(kw, "title") == 0) {
            char buf[64];
            const char *p = rest;
            if (parse_quoted(&p, buf, sizeof buf))
                snprintf(s->title, sizeof s->title, "%s", buf);
        } else if (strcmp(kw, "source") == 0) {
            char buf[1024];
            const char *p = rest;
            if (parse_quoted(&p, buf, sizeof buf))
                snprintf(s->source_path, sizeof s->source_path, "%s", buf);
        } else if (strcmp(kw, "scenario_path") == 0) {
            char buf[1024];
            const char *p = rest;
            if (parse_quoted(&p, buf, sizeof buf))
                snprintf(s->scenario_path, sizeof s->scenario_path, "%s", buf);
        } else if (strcmp(kw, "origin") == 0) {
            sscanf(rest, "%lf %lf", &s->origin_lon, &s->origin_lat);
        } else if (strcmp(kw, "camera") == 0) {
            block = BLK_CAMERA;
        } else if (strcmp(kw, "layers") == 0) {
            block = BLK_LAYERS;
            s->layer_mask = 0;          /* explicit toggles fill it */
        } else if (strcmp(kw, "route") == 0) {
            block = BLK_ROUTE;
        } else if (strcmp(kw, "annotation") == 0) {
            block = BLK_ANN;
            /* `rest` may carry an inline kind like "annotation line".
             * If empty, the body's first token is the kind. */
            if (*rest)
                strncpy(ann_kind, rest, sizeof ann_kind - 1);
            else
                ann_kind[0] = '\0';
        } else if (strcmp(kw, "end") == 0) {
            /* Closing the whole document. */
            break;
        }
        /* Unknown top-level keywords are silently skipped — forward compat. */
    }

    fclose(f);
    return saw_magic ? rc : K26GEO_ERR_PARSE;
}
