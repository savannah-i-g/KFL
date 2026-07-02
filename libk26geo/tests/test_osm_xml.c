/* Smoke test for the .osm XML loader.
 *
 * Builds a tiny .osm document in /tmp, parses it via libk26geo, and
 * walks the resulting K26GeoOsmDoc to confirm nodes + ways + tags
 * round-trip cleanly. No external fixtures required. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "k26geo.h"

static const char FIXTURE[] =
    "<?xml version='1.0' encoding='UTF-8'?>\n"
    "<osm version='0.6' generator='test'>\n"
    "  <bounds minlat='51.500' minlon='-0.130' maxlat='51.510' maxlon='-0.120'/>\n"
    "  <node id='1' lat='51.5010' lon='-0.1280'/>\n"
    "  <node id='2' lat='51.5012' lon='-0.1278'/>\n"
    "  <node id='3' lat='51.5014' lon='-0.1276'/>\n"
    "  <node id='4' lat='51.5016' lon='-0.1274'/>\n"
    "  <node id='100' lat='51.5020' lon='-0.1270'>\n"
    "    <tag k='amenity' v='cafe'/>\n"
    "    <tag k='name' v='Test Cafe'/>\n"
    "  </node>\n"
    "  <way id='10'>\n"
    "    <nd ref='1'/>\n"
    "    <nd ref='2'/>\n"
    "    <nd ref='3'/>\n"
    "    <tag k='highway' v='residential'/>\n"
    "    <tag k='name' v='Test Road'/>\n"
    "  </way>\n"
    "  <way id='11'>\n"
    "    <nd ref='2'/>\n"
    "    <nd ref='3'/>\n"
    "    <nd ref='4'/>\n"
    "    <nd ref='2'/>\n"
    "    <tag k='building' v='yes'/>\n"
    "  </way>\n"
    "</osm>\n";

static int fail = 0;

static int saw_road    = 0;
static int saw_bldg    = 0;
static int saw_cafe    = 0;
static int n_node_calls = 0;
static int n_way_calls  = 0;

static int on_node(int64_t id, double lon, double lat,
                   const char *const *keys, const char *const *vals,
                   size_t n_tags, void *user)
{
    (void)user;
    n_node_calls++;
    if (id == 100) {
        if (lat < 51.5019 || lat > 51.5021) { fprintf(stderr, "FAIL node 100 lat\n"); fail++; }
        if (lon > -0.1269 || lon < -0.1271) { fprintf(stderr, "FAIL node 100 lon\n"); fail++; }
        for (size_t i = 0; i < n_tags; i++) {
            if (strcmp(keys[i], "amenity") == 0 && strcmp(vals[i], "cafe") == 0) saw_cafe = 1;
        }
    }
    return 0;
}

static int on_way(int64_t id, const int64_t *refs, size_t n_refs,
                  const char *const *keys, const char *const *vals,
                  size_t n_tags, void *user)
{
    (void)user;
    n_way_calls++;
    int has_highway = 0, has_building = 0;
    for (size_t i = 0; i < n_tags; i++) {
        if (strcmp(keys[i], "highway")  == 0) has_highway  = 1;
        if (strcmp(keys[i], "building") == 0) has_building = 1;
    }
    if (id == 10) {
        if (n_refs != 3) { fprintf(stderr, "FAIL way 10 refs %zu\n", n_refs); fail++; }
        if (has_highway) saw_road = 1;
    }
    if (id == 11) {
        if (n_refs != 4) { fprintf(stderr, "FAIL way 11 refs %zu\n", n_refs); fail++; }
        if (has_building) saw_bldg = 1;
    }
    return 0;
}

int main(void)
{
    const char *path = "/tmp/test_osm_xml.osm";
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "fopen %s failed\n", path); return 1; }
    fwrite(FIXTURE, 1, sizeof FIXTURE - 1, f);
    fclose(f);

    K26GeoOsmDoc *d = NULL;
    K26GeoStatus s = k26geo_osm_load_file(path, &d);
    if (s != K26GEO_OK || !d) {
        fprintf(stderr, "FAIL load: %s\n", k26geo_status_str(s));
        return 1;
    }

    if (k26geo_osm_n_nodes(d) != 5) {
        fprintf(stderr, "FAIL node count: got %zu want 5\n", k26geo_osm_n_nodes(d)); fail++;
    }
    if (k26geo_osm_n_ways(d) != 2) {
        fprintf(stderr, "FAIL way count: got %zu want 2\n",  k26geo_osm_n_ways(d));  fail++;
    }

    /* Sort-by-id should let us look up node 2. */
    double lon, lat;
    if (!k26geo_osm_node_lookup(d, 2, &lon, &lat)) {
        fprintf(stderr, "FAIL node 2 lookup\n"); fail++;
    } else if (lat < 51.5011 || lat > 51.5013) {
        fprintf(stderr, "FAIL node 2 lat: %.6f\n", lat); fail++;
    }

    /* bbox from <bounds>. */
    double mnlon, mnlat, mxlon, mxlat;
    k26geo_osm_bbox(d, &mnlon, &mnlat, &mxlon, &mxlat);
    if (mnlat > 51.501 || mxlat < 51.509) {
        fprintf(stderr, "FAIL bbox lat: %.3f..%.3f\n", mnlat, mxlat); fail++;
    }

    k26geo_osm_for_each_node(d, on_node, NULL);
    k26geo_osm_for_each_way (d, on_way,  NULL);

    if (n_node_calls != 5) {
        fprintf(stderr, "FAIL node-iter count %d\n", n_node_calls); fail++;
    }
    if (n_way_calls != 2) {
        fprintf(stderr, "FAIL way-iter count %d\n", n_way_calls); fail++;
    }
    if (!saw_road) { fprintf(stderr, "FAIL did not see highway way\n"); fail++; }
    if (!saw_bldg) { fprintf(stderr, "FAIL did not see building way\n"); fail++; }
    if (!saw_cafe) { fprintf(stderr, "FAIL did not see cafe tag\n"); fail++; }

    k26geo_osm_free(d);
    remove(path);

    if (fail) {
        fprintf(stderr, "%d failure(s)\n", fail);
        return 1;
    }
    puts("test_osm_xml: OK");
    return 0;
}
