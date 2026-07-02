/* .osm XML loader — single-pass expat SAX parser.
 *
 * The OSM XML format is simple:
 *
 *     <osm version="0.6" generator="...">
 *       <bounds minlat=".." minlon=".." maxlat=".." maxlon=".."/>
 *       <node id="N" lat="L" lon="L" .../>
 *         <tag k="K" v="V"/>
 *       </node>
 *       <way id="W" ...>
 *         <nd ref="N"/>
 *         <tag k="K" v="V"/>
 *       </way>
 *       <relation id="R" ...>
 *         <member type="..." ref="..." role="..."/>
 *         <tag k="K" v="V"/>
 *       </relation>
 *     </osm>
 *
 * Strings are arena-owned (one big bump allocator). Nodes are sorted
 * by id post-parse to support binary-search lookup from ways.
 *
 * Relations are captured during parse with their member refs, roles,
 * and tags; the public surface exposes them via
 * k26geo_osm_for_each_relation. */

#include <ctype.h>
#include <expat.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>     /* transparent gzip-stream + raw-file reader */

#include "k26geo.h"

/* ---- arena ------------------------------------------------------ */

#define ARENA_CHUNK   (256 * 1024)

typedef struct ArenaChunk {
    struct ArenaChunk *next;
    size_t             used, cap;
    char               data[];
} ArenaChunk;

typedef struct {
    ArenaChunk *head;
} Arena;

static void arena_init(Arena *a) { a->head = NULL; }

static void arena_free(Arena *a)
{
    ArenaChunk *c = a->head;
    while (c) { ArenaChunk *n = c->next; free(c); c = n; }
    a->head = NULL;
}

static void *arena_alloc(Arena *a, size_t n)
{
    if (n == 0) n = 1;
    n = (n + 7) & ~(size_t)7;   /* 8-byte align */
    if (!a->head || a->head->used + n > a->head->cap) {
        size_t cap = ARENA_CHUNK;
        if (n > cap) cap = n + 64;
        ArenaChunk *c = (ArenaChunk *)malloc(sizeof *c + cap);
        if (!c) return NULL;
        c->next = a->head;
        c->used = 0;
        c->cap  = cap;
        a->head = c;
    }
    void *p = a->head->data + a->head->used;
    a->head->used += n;
    return p;
}

static char *arena_strdup(Arena *a, const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *p = (char *)arena_alloc(a, n);
    if (!p) return NULL;
    memcpy(p, s, n);
    return p;
}

/* ---- doc structures --------------------------------------------- */

typedef struct {
    const char *k;
    const char *v;
} OsmTag;

typedef struct {
    int64_t     id;
    double      lon, lat;
    OsmTag     *tags;
    size_t      n_tags;
} OsmNode;

typedef struct {
    int64_t     id;
    int64_t    *refs;
    size_t      n_refs;
    OsmTag     *tags;
    size_t      n_tags;
} OsmWay;

typedef struct {
    int64_t       id;
    int64_t      *member_refs;
    const char  **member_roles;
    const char  **member_types;
    size_t        n_members;
    OsmTag       *tags;
    size_t        n_tags;
} OsmRelation;

struct K26GeoOsmDoc {
    Arena    arena;

    OsmNode *nodes;
    size_t   n_nodes, cap_nodes;

    OsmWay  *ways;
    size_t   n_ways, cap_ways;

    OsmRelation *relations;
    size_t       n_relations, cap_relations;

    double   bbox_min_lon, bbox_min_lat;
    double   bbox_max_lon, bbox_max_lat;
    int      bbox_set;

    /* Scratch state during parse. */
    enum { IN_NONE = 0, IN_NODE, IN_WAY, IN_RELATION } cur;
    OsmNode      cur_node;     /* id+lon+lat valid only when cur == IN_NODE */
    OsmWay       cur_way;      /* id valid only when cur == IN_WAY */
    OsmRelation  cur_rel;      /* id valid only when cur == IN_RELATION */
    OsmTag       tmp_tags[256];
    size_t       n_tmp_tags;
    int64_t      tmp_refs[16384];
    size_t       n_tmp_refs;
    /* Per-relation member buffers — separate from way-node-refs so a
     * <relation> can contain its own <nd> (it can't, but defensive). */
    int64_t      tmp_member_refs[16384];
    const char  *tmp_member_roles[16384];
    const char  *tmp_member_types[16384];
    size_t       n_tmp_members;

    int      had_error;
};

/* ---- helpers ---------------------------------------------------- */

static int64_t parse_i64(const char *s)
{
    if (!s) return 0;
    return (int64_t)strtoll(s, NULL, 10);
}

static double parse_d(const char *s, double dflt)
{
    if (!s) return dflt;
    char *e = NULL;
    double v = strtod(s, &e);
    return (e == s) ? dflt : v;
}

static const char *attr(const char **atts, const char *key)
{
    for (size_t i = 0; atts[i]; i += 2) {
        if (strcmp(atts[i], key) == 0) return atts[i + 1];
    }
    return NULL;
}

static void push_node(K26GeoOsmDoc *d, OsmNode *n)
{
    if (d->n_nodes == d->cap_nodes) {
        size_t cap = d->cap_nodes ? d->cap_nodes * 2 : 1024;
        OsmNode *p = (OsmNode *)realloc(d->nodes, cap * sizeof *p);
        if (!p) { d->had_error = 1; return; }
        d->nodes = p;
        d->cap_nodes = cap;
    }
    d->nodes[d->n_nodes++] = *n;
}

static void push_way(K26GeoOsmDoc *d, OsmWay *w)
{
    if (d->n_ways == d->cap_ways) {
        size_t cap = d->cap_ways ? d->cap_ways * 2 : 1024;
        OsmWay *p = (OsmWay *)realloc(d->ways, cap * sizeof *p);
        if (!p) { d->had_error = 1; return; }
        d->ways = p;
        d->cap_ways = cap;
    }
    d->ways[d->n_ways++] = *w;
}

static void push_relation(K26GeoOsmDoc *d, OsmRelation *r)
{
    if (d->n_relations == d->cap_relations) {
        size_t cap = d->cap_relations ? d->cap_relations * 2 : 256;
        OsmRelation *p = (OsmRelation *)realloc(d->relations, cap * sizeof *p);
        if (!p) { d->had_error = 1; return; }
        d->relations = p;
        d->cap_relations = cap;
    }
    d->relations[d->n_relations++] = *r;
}

/* Order-by-id comparators for bsearch + qsort. */
static int cmp_node_id(const void *a, const void *b)
{
    const OsmNode *na = (const OsmNode *)a;
    const OsmNode *nb = (const OsmNode *)b;
    if (na->id < nb->id) return -1;
    if (na->id > nb->id) return  1;
    return 0;
}

static int cmp_way_id(const void *a, const void *b)
{
    const OsmWay *wa = (const OsmWay *)a;
    const OsmWay *wb = (const OsmWay *)b;
    if (wa->id < wb->id) return -1;
    if (wa->id > wb->id) return  1;
    return 0;
}

/* ---- expat callbacks -------------------------------------------- */

static void XMLCALL on_start(void *user, const XML_Char *name, const XML_Char **atts)
{
    K26GeoOsmDoc *d = (K26GeoOsmDoc *)user;

    if (strcmp(name, "node") == 0) {
        d->cur = IN_NODE;
        d->cur_node.id       = parse_i64(attr(atts, "id"));
        d->cur_node.lon      = parse_d(attr(atts, "lon"), 0.0);
        d->cur_node.lat      = parse_d(attr(atts, "lat"), 0.0);
        d->cur_node.tags     = NULL;
        d->cur_node.n_tags   = 0;
        d->n_tmp_tags = 0;
    } else if (strcmp(name, "way") == 0) {
        d->cur = IN_WAY;
        d->cur_way.id       = parse_i64(attr(atts, "id"));
        d->cur_way.refs     = NULL;
        d->cur_way.n_refs   = 0;
        d->cur_way.tags     = NULL;
        d->cur_way.n_tags   = 0;
        d->n_tmp_tags = 0;
        d->n_tmp_refs = 0;
    } else if (strcmp(name, "relation") == 0) {
        d->cur = IN_RELATION;
        d->cur_rel.id          = parse_i64(attr(atts, "id"));
        d->cur_rel.member_refs = NULL;
        d->cur_rel.member_roles= NULL;
        d->cur_rel.member_types= NULL;
        d->cur_rel.n_members   = 0;
        d->cur_rel.tags        = NULL;
        d->cur_rel.n_tags      = 0;
        d->n_tmp_tags    = 0;
        d->n_tmp_members = 0;
    } else if (strcmp(name, "tag") == 0) {
        if (d->cur == IN_NODE || d->cur == IN_WAY || d->cur == IN_RELATION) {
            const char *k = attr(atts, "k");
            const char *v = attr(atts, "v");
            if (k && v && d->n_tmp_tags < sizeof d->tmp_tags / sizeof d->tmp_tags[0]) {
                d->tmp_tags[d->n_tmp_tags].k = arena_strdup(&d->arena, k);
                d->tmp_tags[d->n_tmp_tags].v = arena_strdup(&d->arena, v);
                if (!d->tmp_tags[d->n_tmp_tags].k || !d->tmp_tags[d->n_tmp_tags].v) {
                    d->had_error = 1;
                } else {
                    d->n_tmp_tags++;
                }
            }
        }
    } else if (strcmp(name, "nd") == 0) {
        if (d->cur == IN_WAY) {
            const char *ref = attr(atts, "ref");
            if (ref && d->n_tmp_refs < sizeof d->tmp_refs / sizeof d->tmp_refs[0]) {
                d->tmp_refs[d->n_tmp_refs++] = parse_i64(ref);
            }
        }
    } else if (strcmp(name, "member") == 0) {
        if (d->cur == IN_RELATION) {
            const char *type = attr(atts, "type");
            const char *ref  = attr(atts, "ref");
            const char *role = attr(atts, "role");
            if (type && ref &&
                d->n_tmp_members < sizeof d->tmp_member_refs / sizeof d->tmp_member_refs[0]) {
                d->tmp_member_refs [d->n_tmp_members] = parse_i64(ref);
                d->tmp_member_types[d->n_tmp_members] = arena_strdup(&d->arena, type);
                d->tmp_member_roles[d->n_tmp_members] = role ? arena_strdup(&d->arena, role) : "";
                if (!d->tmp_member_types[d->n_tmp_members]) {
                    d->had_error = 1;
                } else {
                    d->n_tmp_members++;
                }
            }
        }
    } else if (strcmp(name, "bounds") == 0) {
        d->bbox_min_lon = parse_d(attr(atts, "minlon"),  180.0);
        d->bbox_min_lat = parse_d(attr(atts, "minlat"),   90.0);
        d->bbox_max_lon = parse_d(attr(atts, "maxlon"), -180.0);
        d->bbox_max_lat = parse_d(attr(atts, "maxlat"),  -90.0);
        d->bbox_set     = 1;
    }
    /* member, osm, etc — ignored. */
}

static void XMLCALL on_end(void *user, const XML_Char *name)
{
    K26GeoOsmDoc *d = (K26GeoOsmDoc *)user;

    if (strcmp(name, "node") == 0 && d->cur == IN_NODE) {
        /* Copy any tags into arena. */
        if (d->n_tmp_tags > 0) {
            OsmTag *t = (OsmTag *)arena_alloc(&d->arena, d->n_tmp_tags * sizeof *t);
            if (!t) { d->had_error = 1; }
            else {
                memcpy(t, d->tmp_tags, d->n_tmp_tags * sizeof *t);
                d->cur_node.tags   = t;
                d->cur_node.n_tags = d->n_tmp_tags;
            }
        }
        push_node(d, &d->cur_node);
        d->cur = IN_NONE;

        /* Accumulate bbox from nodes if no explicit <bounds>. */
        if (!d->bbox_set) {
            if (d->n_nodes == 1) {
                d->bbox_min_lon = d->bbox_max_lon = d->cur_node.lon;
                d->bbox_min_lat = d->bbox_max_lat = d->cur_node.lat;
            } else {
                if (d->cur_node.lon < d->bbox_min_lon) d->bbox_min_lon = d->cur_node.lon;
                if (d->cur_node.lon > d->bbox_max_lon) d->bbox_max_lon = d->cur_node.lon;
                if (d->cur_node.lat < d->bbox_min_lat) d->bbox_min_lat = d->cur_node.lat;
                if (d->cur_node.lat > d->bbox_max_lat) d->bbox_max_lat = d->cur_node.lat;
            }
        }
    } else if (strcmp(name, "way") == 0 && d->cur == IN_WAY) {
        if (d->n_tmp_tags > 0) {
            OsmTag *t = (OsmTag *)arena_alloc(&d->arena, d->n_tmp_tags * sizeof *t);
            if (!t) { d->had_error = 1; }
            else {
                memcpy(t, d->tmp_tags, d->n_tmp_tags * sizeof *t);
                d->cur_way.tags   = t;
                d->cur_way.n_tags = d->n_tmp_tags;
            }
        }
        if (d->n_tmp_refs > 0) {
            int64_t *r = (int64_t *)arena_alloc(&d->arena, d->n_tmp_refs * sizeof *r);
            if (!r) { d->had_error = 1; }
            else {
                memcpy(r, d->tmp_refs, d->n_tmp_refs * sizeof *r);
                d->cur_way.refs   = r;
                d->cur_way.n_refs = d->n_tmp_refs;
            }
        }
        push_way(d, &d->cur_way);
        d->cur = IN_NONE;
    } else if (strcmp(name, "relation") == 0 && d->cur == IN_RELATION) {
        /* Stash tags + members into arena-owned storage. */
        if (d->n_tmp_tags > 0) {
            OsmTag *t = (OsmTag *)arena_alloc(&d->arena, d->n_tmp_tags * sizeof *t);
            if (!t) { d->had_error = 1; }
            else {
                memcpy(t, d->tmp_tags, d->n_tmp_tags * sizeof *t);
                d->cur_rel.tags   = t;
                d->cur_rel.n_tags = d->n_tmp_tags;
            }
        }
        if (d->n_tmp_members > 0) {
            int64_t      *r = (int64_t *)arena_alloc(&d->arena,
                d->n_tmp_members * sizeof *r);
            const char  **rl = (const char **)arena_alloc(&d->arena,
                d->n_tmp_members * sizeof *rl);
            const char  **tp = (const char **)arena_alloc(&d->arena,
                d->n_tmp_members * sizeof *tp);
            if (!r || !rl || !tp) { d->had_error = 1; }
            else {
                memcpy(r,  d->tmp_member_refs,  d->n_tmp_members * sizeof *r);
                memcpy(rl, d->tmp_member_roles, d->n_tmp_members * sizeof *rl);
                memcpy(tp, d->tmp_member_types, d->n_tmp_members * sizeof *tp);
                d->cur_rel.member_refs  = r;
                d->cur_rel.member_roles = rl;
                d->cur_rel.member_types = tp;
                d->cur_rel.n_members    = d->n_tmp_members;
            }
        }
        push_relation(d, &d->cur_rel);
        d->cur = IN_NONE;
    }
}

/* ---- public API ------------------------------------------------- */

K26GeoStatus k26geo_osm_load_file(const char *path, K26GeoOsmDoc **out)
{
    if (!path || !out) return K26GEO_ERR_INVAL;
    *out = NULL;

    /* zlib's gzopen reads .gz files AND plain files transparently:
     * gzread on a non-gzipped stream returns the raw bytes. So a single
     * code path handles both .osm and .osm.gz. */
    gzFile fp = gzopen(path, "rb");
    if (!fp) return K26GEO_ERR_IO;

    K26GeoOsmDoc *d = (K26GeoOsmDoc *)calloc(1, sizeof *d);
    if (!d) { gzclose(fp); return K26GEO_ERR_OOM; }
    arena_init(&d->arena);

    XML_Parser p = XML_ParserCreate(NULL);
    if (!p) { free(d); gzclose(fp); return K26GEO_ERR_OOM; }
    XML_SetUserData(p, d);
    XML_SetElementHandler(p, on_start, on_end);

    char buf[64 * 1024];
    int done = 0;
    K26GeoStatus rc = K26GEO_OK;
    while (!done && !d->had_error) {
        int r = gzread(fp, buf, (int)sizeof buf);
        if (r < 0) { rc = K26GEO_ERR_IO; break; }
        done = (r < (int)sizeof buf);
        if (XML_Parse(p, buf, r, done) == XML_STATUS_ERROR) {
            rc = K26GEO_ERR_PARSE;
            break;
        }
    }
    if (rc == K26GEO_OK && d->had_error) rc = K26GEO_ERR_OOM;

    XML_ParserFree(p);
    gzclose(fp);

    if (rc != K26GEO_OK) { k26geo_osm_free(d); return rc; }

    /* Sort nodes + ways for binary-search lookup. Relations are walked
     * by iteration only so we leave their order untouched (preserving
     * the file's authoring sequence helps debug bad multipolygon
     * reconstructions when they happen). */
    if (d->n_nodes > 1) qsort(d->nodes, d->n_nodes, sizeof *d->nodes, cmp_node_id);
    if (d->n_ways  > 1) qsort(d->ways,  d->n_ways,  sizeof *d->ways,  cmp_way_id);

    *out = d;
    return K26GEO_OK;
}

void k26geo_osm_free(K26GeoOsmDoc *d)
{
    if (!d) return;
    arena_free(&d->arena);
    free(d->nodes);
    free(d->ways);
    free(d->relations);
    free(d);
}

size_t k26geo_osm_n_nodes    (const K26GeoOsmDoc *d) { return d ? d->n_nodes     : 0; }
size_t k26geo_osm_n_ways     (const K26GeoOsmDoc *d) { return d ? d->n_ways      : 0; }
size_t k26geo_osm_n_relations(const K26GeoOsmDoc *d) { return d ? d->n_relations : 0; }

int k26geo_osm_node_lookup(const K26GeoOsmDoc *d, int64_t id,
                           double *out_lon, double *out_lat)
{
    if (!d) return 0;
    /* Binary search over sorted-by-id nodes. */
    size_t lo = 0, hi = d->n_nodes;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (d->nodes[mid].id < id)      lo = mid + 1;
        else if (d->nodes[mid].id > id) hi = mid;
        else {
            if (out_lon) *out_lon = d->nodes[mid].lon;
            if (out_lat) *out_lat = d->nodes[mid].lat;
            return 1;
        }
    }
    return 0;
}

void k26geo_osm_bbox(const K26GeoOsmDoc *d,
                     double *out_min_lon, double *out_min_lat,
                     double *out_max_lon, double *out_max_lat)
{
    if (!d) return;
    if (out_min_lon) *out_min_lon = d->bbox_min_lon;
    if (out_min_lat) *out_min_lat = d->bbox_min_lat;
    if (out_max_lon) *out_max_lon = d->bbox_max_lon;
    if (out_max_lat) *out_max_lat = d->bbox_max_lat;
}

void k26geo_osm_for_each_node(const K26GeoOsmDoc *d, K26GeoOsmNodeFn fn, void *user)
{
    if (!d || !fn) return;
    /* Build the keys[] / vals[] aliases on the fly into a small stack buf
     * (most nodes carry 0-5 tags). */
    enum { TAG_BUF = 64 };
    const char *keys[TAG_BUF];
    const char *vals[TAG_BUF];
    for (size_t i = 0; i < d->n_nodes; i++) {
        const OsmNode *n = &d->nodes[i];
        size_t nt = n->n_tags;
        if (nt > TAG_BUF) nt = TAG_BUF;
        for (size_t j = 0; j < nt; j++) {
            keys[j] = n->tags[j].k;
            vals[j] = n->tags[j].v;
        }
        if (fn(n->id, n->lon, n->lat, keys, vals, nt, user)) return;
    }
}

void k26geo_osm_for_each_way(const K26GeoOsmDoc *d, K26GeoOsmWayFn fn, void *user)
{
    if (!d || !fn) return;
    enum { TAG_BUF = 128 };
    const char *keys[TAG_BUF];
    const char *vals[TAG_BUF];
    for (size_t i = 0; i < d->n_ways; i++) {
        const OsmWay *w = &d->ways[i];
        size_t nt = w->n_tags;
        if (nt > TAG_BUF) nt = TAG_BUF;
        for (size_t j = 0; j < nt; j++) {
            keys[j] = w->tags[j].k;
            vals[j] = w->tags[j].v;
        }
        if (fn(w->id, w->refs, w->n_refs, keys, vals, nt, user)) return;
    }
}

int k26geo_osm_way_lookup(const K26GeoOsmDoc *d, int64_t id,
                          const int64_t **out_refs, size_t *out_n_refs,
                          const char *const **out_keys,
                          const char *const **out_vals,
                          size_t *out_n_tags)
{
    if (!d) return 0;
    size_t lo = 0, hi = d->n_ways;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if      (d->ways[mid].id < id) lo = mid + 1;
        else if (d->ways[mid].id > id) hi = mid;
        else {
            const OsmWay *w = &d->ways[mid];
            if (out_refs)   *out_refs   = w->refs;
            if (out_n_refs) *out_n_refs = w->n_refs;
            /* Tag accessors fold into a thread-local alias buf at the
             * caller side; for relation-assembly we only need refs. To
             * keep this hot path branch-free we return NULL for the
             * key/val aliases — callers who need tags call the iterator
             * instead. */
            if (out_keys)   *out_keys   = NULL;
            if (out_vals)   *out_vals   = NULL;
            if (out_n_tags) *out_n_tags = w->n_tags;
            return 1;
        }
    }
    return 0;
}

void k26geo_osm_for_each_relation(const K26GeoOsmDoc *d,
                                  K26GeoOsmRelationFn fn,
                                  void *user)
{
    if (!d || !fn) return;
    enum { TAG_BUF = 128 };
    const char *keys[TAG_BUF];
    const char *vals[TAG_BUF];
    for (size_t i = 0; i < d->n_relations; i++) {
        const OsmRelation *r = &d->relations[i];
        size_t nt = r->n_tags;
        if (nt > TAG_BUF) nt = TAG_BUF;
        for (size_t j = 0; j < nt; j++) {
            keys[j] = r->tags[j].k;
            vals[j] = r->tags[j].v;
        }
        if (fn(r->id, r->member_refs, r->member_roles, r->member_types,
               r->n_members, keys, vals, nt, user)) return;
    }
}
