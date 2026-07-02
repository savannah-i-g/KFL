/* libk26geo — geographic primitives for K26 mapping / visualisation apps.
 *
 * This library provides:
 *   - Web Mercator tangent-plane projection (forward + inverse).
 *   - Ramer-Douglas-Peucker polyline simplification.
 *   - Uniform GxG spatial grid index.
 *   - OSM tag parsers: `height`, `building:levels`, `min_height`,
 *     `building:min_level`, `width`, `lanes`, highway-class table.
 *   - `.osm` XML loader (expat SAX, single-pass) with node, way, and
 *     relation iterators plus id-based lookup.
 *   - `.k26atlas` bespoke text session-file reader + emitter.
 *   - `.k26bake` opaque binary mesh cache reader + emitter.
 *
 * Determinism contract:
 *   - All numeric routines are bit-reproducible under `-ffp-contract=off`
 *     `-fexcess-precision=standard`, matching libk26m3d.
 *   - Mercator round-trip is below 1e-6 m at lat |phi| <= 60 deg on
 *     standard IEEE-754.
 *
 * Not provided:
 *   - `.osm.pbf` (Protocol Buffers) loader.
 *   - Quadtree / R-tree spatial indices.
 *   - WGS84 ellipsoid (the implementation uses spherical Mercator with
 *     R=6378137).
 *   - Roof shape geometry (gabled / hipped / etc.; owned by k26atlas).
 *   - `building:part` relation resolution.
 *
 * Threading: every function below is reentrant. Opaque handles are
 * owned by exactly one thread at a time; the library performs no
 * internal locking. */
#ifndef K26GEO_H
#define K26GEO_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

#define K26GEO_LIB_VERSION         "0.1.0"
#define K26ATLAS_FORMAT_VERSION    1
#define K26BAKE_FORMAT_VERSION     1
#define K26GEO_EARTH_RADIUS_M      6378137.0
#define K26GEO_MAX_ABS_LAT_DEG     70.0  /* polar Mercator distortion cap */

/* ---- Status codes ----------------------------------------------- */
typedef enum {
    K26GEO_OK            = 0,
    K26GEO_ERR_INVAL,
    K26GEO_ERR_OOM,
    K26GEO_ERR_IO,
    K26GEO_ERR_PARSE,
    K26GEO_ERR_RANGE,
    K26GEO_ERR_UNSUPPORTED,
    K26GEO_ERR_INTERNAL
} K26GeoStatus;

const char *k26geo_status_str(K26GeoStatus s);

/* ================================================================
 * Mercator projection — tangent-plane linearisation at a fixed origin
 * ================================================================ */

typedef struct {
    double lon_deg;     /* origin longitude (degrees) */
    double lat_deg;     /* origin latitude  (degrees) */
    double merc_y0;     /* cached log(tan(π/4 + lat·π/360)) */
    double cos_lat0;    /* cached cos(lat·π/180); zero/near-zero at poles */
    double scale;       /* R · cos(lat0) — common Mercator scale factor   */
} K26GeoOrigin;

/* Initialise a tangent-plane origin. Returns K26GEO_ERR_RANGE if
 * |lat| > K26GEO_MAX_ABS_LAT_DEG. */
K26GeoStatus k26geo_origin_init(K26GeoOrigin *out, double lon_deg, double lat_deg);

/* Forward projection: lon/lat (degrees) → world-space metres.
 * World axes match libk26m3d's right-handed convention:
 *   world.x =  +east_metres
 *   world.y =   0 (caller sets altitude)
 *   world.z = -north_metres   (north into the screen with the default cam) */
K26V3 k26geo_lonlat_to_world(const K26GeoOrigin *o, double lon_deg, double lat_deg);

/* Inverse projection: world (x, z) → lon/lat (degrees). The y
 * coordinate is ignored. Outputs are clamped to ±180° / ±85°. */
void  k26geo_world_to_lonlat(const K26GeoOrigin *o, double world_x, double world_z,
                             double *out_lon_deg, double *out_lat_deg);

/* Convenience: spherical Mercator y-component for a bare latitude. */
double k26geo_merc_y(double lat_deg);

/* ================================================================
 * Ramer-Douglas-Peucker simplification (in-place)
 * ================================================================ */

/* Operates on packed (x, z) pairs in `coords[2*n]`. Endpoints always
 * preserved. Iterative stack-based; never recurses. Epsilon in the
 * same units as the input (metres if fed from k26geo_lonlat_to_world).
 * Returns the new point count (≤ n). epsilon ≤ 0 is a no-op. */
size_t k26geo_rdp_xz(double *coords, size_t n, double epsilon);

/* Preset epsilons (world metres) matching SentinelMapper LOSSLESS /
 * HIGH / MEDIUM / LOW. Documented constants — change deliberately. */
#define K26GEO_RDP_LOSSLESS_M  0.0
#define K26GEO_RDP_HIGH_M      0.5
#define K26GEO_RDP_MEDIUM_M    2.0
#define K26GEO_RDP_LOW_M       8.0

/* ================================================================
 * Sutherland-Hodgman polygon clipping against an axis-aligned bbox
 * ================================================================ */

/* Clip a closed polygon (packed `n_in` (x, z) pairs in `in_coords`)
 * against the bbox rect. Writes the clipped polygon to `out_coords`
 * (a caller-owned buffer holding up to `cap_out` doubles, i.e.
 * `cap_out / 2` vertices). Returns the new vertex count, or 0 if
 * the polygon is fully outside the bbox or the output buffer is
 * too small.
 *
 * Safe to call with `in_coords == out_coords` — the routine
 * scratches into an internal stack buffer. Worst-case output size is
 * 2 × n_in vertices (each clip edge can double the polygon's
 * vertex count); a `cap_out` of `2 * n_in * 2` doubles is always
 * sufficient. */
size_t k26geo_clip_polygon_xz(const double *in_coords, size_t n_in,
                              double *out_coords, size_t cap_out,
                              double min_x, double min_z,
                              double max_x, double max_z);

/* ================================================================
 * Polyline clipping (Liang-Barsky, segment-by-segment)
 *
 * A polyline may enter and exit the bbox multiple times; the
 * callback fires once per in-bbox sub-polyline. The data passed to
 * the callback is owned by the clipper (transient); copy if you
 * need to keep it. Returns the count of sub-polylines emitted.
 * ================================================================ */

typedef void (*K26GeoPolylineCb)(const double *coords_xz, size_t n_pts, void *user);

size_t k26geo_clip_polyline_xz(const double *in_coords, size_t n_in,
                               double min_x, double min_z,
                               double max_x, double max_z,
                               K26GeoPolylineCb cb, void *user);

/* ================================================================
 * Uniform G×G spatial grid
 * ================================================================ */

typedef struct K26GeoSpGrid K26GeoSpGrid;

/* Allocate a grid covering the inclusive rect (min_x, min_z .. max_x,
 * max_z) with n_cells_x × n_cells_z buckets. n_cells_* ≥ 1. */
K26GeoSpGrid *k26geo_spgrid_new (double min_x, double min_z,
                                 double max_x, double max_z,
                                 int n_cells_x, int n_cells_z);
void          k26geo_spgrid_free(K26GeoSpGrid *g);

/* Add a feature id to every cell its AABB overlaps. The feature is
 * stored by id only; coordinates are caller-owned. */
K26GeoStatus k26geo_spgrid_add(K26GeoSpGrid *g, int feature_id,
                               double minx, double minz,
                               double maxx, double maxz);

/* Per-cell visit callback. `feature_ids` may contain duplicates if a
 * caller added the same id to multiple cells; consumers dedup if needed. */
typedef void (*K26GeoSpGridVisitFn)(int cell_x, int cell_z,
                                    const int *feature_ids, size_t n,
                                    void *user);

void k26geo_spgrid_visit_rect(const K26GeoSpGrid *g,
                              double minx, double minz,
                              double maxx, double maxz,
                              K26GeoSpGridVisitFn fn, void *user);

/* Bounds of a single cell. */
void k26geo_spgrid_cell_bounds(const K26GeoSpGrid *g, int cx, int cz,
                               double *out_min_x, double *out_min_z,
                               double *out_max_x, double *out_max_z);

int k26geo_spgrid_n_x(const K26GeoSpGrid *g);
int k26geo_spgrid_n_z(const K26GeoSpGrid *g);

/* ================================================================
 * OSM tag parsing
 * ================================================================ */

/* `height=*` → metres. Accepts:
 *    "10"          → 10.0
 *    "10 m"        → 10.0
 *    "10.5"        → 10.5
 *    "0,8"         → 0.8     (continental comma decimal)
 *    "3m"          → 3.0     (missing space tolerated)
 *    "7'4\""       → 2.2352  (feet + inches, typewriter quotes)
 *    "11' 4\""     → ditto with intra-quote space
 * Rejects: "0.0", "0", negatives, ranges ("3-5"), open inequalities
 * ("<5"). Returns K26GEO_OK + writes *out_m on success. */
K26GeoStatus k26geo_parse_height_m(const char *s, double *out_m);

/* `building:levels=*` → metres at the documented community convention
 * (`k26geo_metres_per_level`, default 3.0 m). Positive integers only;
 * fractional / zero / negative rejected. */
K26GeoStatus k26geo_parse_levels_m(const char *s, double *out_m);

/* Default building height when both `height` and `building:levels`
 * are absent. Compile-time constant K26GEO_DEFAULT_BUILDING_HEIGHT_M. */
#define K26GEO_DEFAULT_BUILDING_HEIGHT_M  6.0
#define K26GEO_METRES_PER_LEVEL           3.0
#define K26GEO_METRES_PER_LANE            3.5

double k26geo_default_building_height_m(void);
double k26geo_metres_per_level(void);
double k26geo_metres_per_lane(void);

/* `width=*` on a highway → metres. Same parser as height with the
 * same tolerance for `"3 m"` etc. Rejects zero/negative. */
K26GeoStatus k26geo_parse_width_m(const char *s, double *out_m);

/* Highway-class width table. Returns metres for known classes, 0.0
 * for unknown.  Values are project-judgement defaults documented in
 * osmtags.c (motorway 14, primary 9, secondary 7, tertiary 6,
 * residential 5, service 4, footway 1.5, ...). */
double k26geo_highway_class_width_m(const char *highway_value);

/* Resolve a road's render width. Precedence:
 *   1. explicit `width=*` if parseable as positive metres
 *   2. `lanes=*` × K26GEO_METRES_PER_LANE
 *   3. class default from k26geo_highway_class_width_m
 * The caller supplies a tag-lookup function (returns the tag value
 * string or NULL). highway_value is the value of the `highway` tag. */
typedef const char *(*K26GeoTagLookupFn)(const char *key, void *user);

double k26geo_resolve_road_width_m(const char *highway_value,
                                   K26GeoTagLookupFn lookup, void *user);

/* ================================================================
 * OSM XML (.osm) loader
 * ================================================================ */

typedef struct K26GeoOsmDoc K26GeoOsmDoc;

/* Stream-parse an .osm XML file via expat. The returned document owns
 * all string data through an internal arena; one free releases all.
 * Returns OOM/IO/PARSE on failure with *out unchanged. */
K26GeoStatus k26geo_osm_load_file(const char *path, K26GeoOsmDoc **out);
void         k26geo_osm_free     (K26GeoOsmDoc *d);

/* Per-record iterators. The arrays passed to the callback are owned
 * by `d` and remain valid until k26geo_osm_free. Return non-zero from
 * the callback to short-circuit the walk. */
typedef int (*K26GeoOsmNodeFn)(int64_t id, double lon, double lat,
                               const char *const *keys,
                               const char *const *vals, size_t n_tags,
                               void *user);
typedef int (*K26GeoOsmWayFn) (int64_t id,
                               const int64_t *node_refs, size_t n_refs,
                               const char *const *keys,
                               const char *const *vals, size_t n_tags,
                               void *user);

void k26geo_osm_for_each_node(const K26GeoOsmDoc *d, K26GeoOsmNodeFn fn, void *user);
void k26geo_osm_for_each_way (const K26GeoOsmDoc *d, K26GeoOsmWayFn  fn, void *user);

/* Relation iterator. Captures only the members and tags needed for
 * multipolygon building reconstruction: member type, ref, and role.
 * Way-typed members are the useful ones for the k26atlas consumer;
 * node and nested-relation members are passed through but are
 * typically ignored by the geometry path. */
typedef int (*K26GeoOsmRelationFn)(int64_t id,
                                   const int64_t   *member_refs,
                                   const char *const *member_roles,
                                   const char *const *member_types,
                                   size_t n_members,
                                   const char *const *keys,
                                   const char *const *vals,
                                   size_t n_tags,
                                   void *user);
void k26geo_osm_for_each_relation(const K26GeoOsmDoc *d,
                                  K26GeoOsmRelationFn fn,
                                  void *user);

/* Resolve a node id to its (lon, lat). Returns 1 on hit, 0 on miss. */
int  k26geo_osm_node_lookup(const K26GeoOsmDoc *d, int64_t id,
                            double *out_lon, double *out_lat);

/* Resolve a way id to its node-ref array and tag table. Out params
 * may be null. Returns 1 on hit, 0 on miss. The returned pointers are
 * owned by `d` and stay valid until k26geo_osm_free. */
int  k26geo_osm_way_lookup(const K26GeoOsmDoc *d, int64_t id,
                           const int64_t **out_refs, size_t *out_n_refs,
                           const char *const **out_keys,
                           const char *const **out_vals,
                           size_t *out_n_tags);

void k26geo_osm_bbox(const K26GeoOsmDoc *d,
                     double *out_min_lon, double *out_min_lat,
                     double *out_max_lon, double *out_max_lat);

size_t k26geo_osm_n_nodes(const K26GeoOsmDoc *d);
size_t k26geo_osm_n_ways (const K26GeoOsmDoc *d);
size_t k26geo_osm_n_relations(const K26GeoOsmDoc *d);

/* ================================================================
 * .k26atlas — bespoke text session-file
 *
 *   atlas v1
 *       title    "Bristol BS1"
 *       source   "data/bristol-bs1.geojson"
 *       origin   51.4544 -2.5879
 *
 *       camera
 *           eye      -120  90  240
 *           target     0    0    0
 *           up         0    1    0
 *           fov_y    55
 *       end
 *
 *       layers
 *           buildings  on
 *           roads      on
 *           water      on
 *           ...
 *       end
 *
 *       route
 *           origin   51.4544 -2.5879
 *           dest     51.4500 -2.5800
 *       end
 *
 *       annotation line
 *           from -42.1  0.0  31.4
 *           to   -28.6  0.0  44.8
 *           color #ff8030  width 2.5
 *       end
 *   end
 *
 * `#` comments allowed. Indentation is cosmetic; the lexer skips
 * leading whitespace.
 * ================================================================ */

#define K26ATLAS_LAYER_BUILDINGS    (1u << 0)
#define K26ATLAS_LAYER_ROADS        (1u << 1)
#define K26ATLAS_LAYER_WATER        (1u << 2)
#define K26ATLAS_LAYER_PATH         (1u << 3)
#define K26ATLAS_LAYER_RAILWAYS     (1u << 4)
#define K26ATLAS_LAYER_BOUNDARIES   (1u << 5)
#define K26ATLAS_LAYER_POIS         (1u << 6)
#define K26ATLAS_LAYER_LABELS       (1u << 7)
#define K26ATLAS_LAYER_ROUTE        (1u << 8)

#define K26ATLAS_DEFAULT_LAYER_MASK 0x1FFu  /* all 9 layers on */

typedef struct {
    char     title[64];
    char     source_path[1024];

    double   origin_lon, origin_lat;

    K26V3    cam_eye;
    K26V3    cam_target;
    K26V3    cam_up;
    double   cam_fov_y_deg;

    uint32_t layer_mask;

    /* (NaN, NaN) → unset */
    double   route_origin_lon, route_origin_lat;
    double   route_dest_lon,   route_dest_lat;

    /* Optional sibling .k26sc scenario reference. Empty string means
     * no scenario. Readers that do not understand the key skip it. */
    char     scenario_path[1024];
} K26AtlasSession;

void         k26atlas_session_defaults(K26AtlasSession *s);

/* Annotation read/write hooks. The library owns the surrounding
 * `annotation <kind>` blocks; the caller's emit_fn writes the body
 * lines, and parse_fn consumes them. Pass NULL hooks to skip
 * annotations entirely (sessions-without-annotations are legal). */
typedef int  (*K26AtlasAnnEmitFn) (FILE *f, size_t index, void *user);
typedef int  (*K26AtlasAnnCountFn)(void *user);
typedef int  (*K26AtlasAnnParseFn)(const char *kind, const char *body, void *user);

K26GeoStatus k26atlas_session_save(const char *path,
                                   const K26AtlasSession *s,
                                   K26AtlasAnnCountFn count_fn,
                                   K26AtlasAnnEmitFn  emit_fn,
                                   void *user);

K26GeoStatus k26atlas_session_load(const char *path,
                                   K26AtlasSession *s,
                                   K26AtlasAnnParseFn parse_fn,
                                   void *user);

/* ================================================================
 * .k26bake — derived binary mesh cache
 *
 * Magic   "K26BAKE\0"
 * Version u32                              (currently 1)
 * Source key 16 bytes                       (SHA-128 of size+mtime+first-4KiB)
 * Origin  double lon, double lat
 * BBox    4 doubles (min_lon, min_lat, max_lon, max_lat)
 * Grid    int32 grid_x, int32 grid_z
 * Entries u32 n_entries
 *   per entry:
 *     u32 tile_id, u32 kind, u64 blob_len, blob bytes
 * ================================================================ */

typedef struct {
    uint32_t   tile_id;
    uint32_t   kind;       /* opaque to libk26geo; app schemas it */
    const void *blob;
    size_t      blob_n;
} K26BakeTileEntry;

typedef struct {
    uint8_t           src_key[16];
    double            origin_lon, origin_lat;
    double            bbox[4];          /* min_lon, min_lat, max_lon, max_lat */
    int32_t           grid_x, grid_z;
    const K26BakeTileEntry *entries;
    size_t            n_entries;
    void             *_internal;        /* arena owning blobs */
} K26BakeFile;

/* Cache key — 128-bit hash of (file size, mtime, first 4 KiB). */
K26GeoStatus k26bake_source_key(const char *src_path, uint8_t out_key[16]);

/* Quick freshness check: opens header, compares key. Returns 1 if
 * matched, 0 if mismatched / unreadable / missing. Cheap; no full read. */
int          k26bake_is_fresh(const char *bake_path, const uint8_t key[16]);

/* Full read into memory (blobs live in an internal arena freed by
 * k26bake_close). */
K26GeoStatus k26bake_open (const char *bake_path, K26BakeFile *out);
void         k26bake_close(K26BakeFile *f);

/* Streaming write: begin_write reserves header + descriptor table by
 * leaving a hole at the front, each write_entry appends an entry,
 * end_write seeks back and patches in the descriptor table. */
typedef struct K26BakeWriter K26BakeWriter;

K26GeoStatus k26bake_begin_write(const char *path,
                                 const uint8_t src_key[16],
                                 double origin_lon, double origin_lat,
                                 const double bbox[4],
                                 int32_t grid_x, int32_t grid_z,
                                 K26BakeWriter **out_w);
K26GeoStatus k26bake_write_entry(K26BakeWriter *w, const K26BakeTileEntry *entry);
K26GeoStatus k26bake_end_write  (K26BakeWriter *w);

#ifdef __cplusplus
}
#endif

#endif /* K26GEO_H */
