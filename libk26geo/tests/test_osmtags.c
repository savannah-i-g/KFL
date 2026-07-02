#include <math.h>
#include <stdio.h>
#include <string.h>

#include "k26geo.h"

static int fail = 0;

static void ok_h(const char *in, double want_m, double tol)
{
    double got = 0.0;
    K26GeoStatus s = k26geo_parse_height_m(in, &got);
    if (s != K26GEO_OK) {
        fprintf(stderr, "FAIL height \"%s\": status %s\n", in, k26geo_status_str(s));
        fail++; return;
    }
    if (fabs(got - want_m) > tol) {
        fprintf(stderr, "FAIL height \"%s\": got %.6f want %.6f\n", in, got, want_m);
        fail++;
    }
}

static void bad_h(const char *in)
{
    double got = 0.0;
    if (k26geo_parse_height_m(in, &got) == K26GEO_OK) {
        fprintf(stderr, "FAIL height \"%s\" should reject (got %.6f)\n", in, got);
        fail++;
    }
}

static void ok_l(const char *in, double want_m)
{
    double got = 0.0;
    K26GeoStatus s = k26geo_parse_levels_m(in, &got);
    if (s != K26GEO_OK || got != want_m) {
        fprintf(stderr, "FAIL levels \"%s\": status=%s got=%g want=%g\n",
                in, k26geo_status_str(s), got, want_m);
        fail++;
    }
}

static void bad_l(const char *in)
{
    double got = 0.0;
    if (k26geo_parse_levels_m(in, &got) == K26GEO_OK) {
        fprintf(stderr, "FAIL levels \"%s\" should reject (got %.6f)\n", in, got);
        fail++;
    }
}

/* Tag lookup helper for resolve_road_width_m */
typedef struct { const char *k, *v; } KV;
static const char *kv_lookup(const char *key, void *user)
{
    KV *arr = (KV *)user;
    for (size_t i = 0; arr[i].k; i++) {
        if (strcmp(arr[i].k, key) == 0) return arr[i].v;
    }
    return NULL;
}

int main(void)
{
    /* --- height --- */
    ok_h("10",            10.0,    1e-9);
    ok_h("10 m",          10.0,    1e-9);
    ok_h("10m",           10.0,    1e-9);
    ok_h("10.5",          10.5,    1e-9);
    ok_h("0,8",            0.8,    1e-9);
    ok_h("3m",             3.0,    1e-9);
    ok_h("  4.2  ",        4.2,    1e-9);
    ok_h("7'4\"",          2.2352, 1e-4);
    ok_h("11' 4\"",        3.4544, 1e-4);
    ok_h("7'",             2.1336, 1e-4);
    ok_h("5 ft",           1.524,  1e-4);

    bad_h("");
    bad_h("0");
    bad_h("0.0");
    bad_h("-5");
    bad_h("3-5");
    bad_h(">20");
    bad_h("<5");
    bad_h("ten");
    bad_h("5 inches");

    /* --- levels --- */
    ok_l("1",  3.0);
    ok_l("2",  6.0);
    ok_l("12", 36.0);
    ok_l(" 4 ", 12.0);

    bad_l("");
    bad_l("0");
    bad_l("-1");
    bad_l("1.5");      /* fractional rejected per wiki */
    bad_l("two");

    /* --- highway class widths --- */
    if (k26geo_highway_class_width_m("motorway") != 14.0) {
        fprintf(stderr, "FAIL motorway width\n"); fail++;
    }
    if (k26geo_highway_class_width_m("footway") != 1.5) {
        fprintf(stderr, "FAIL footway width\n"); fail++;
    }
    if (k26geo_highway_class_width_m("__not_real__") != 0.0) {
        fprintf(stderr, "FAIL unknown class should be 0.0\n"); fail++;
    }

    /* --- resolve_road_width_m --- */
    {
        KV nothing[] = { { NULL, NULL } };
        if (k26geo_resolve_road_width_m("primary", kv_lookup, nothing) != 9.0) {
            fprintf(stderr, "FAIL primary fallback\n"); fail++;
        }
    }
    {
        KV with_width[] = { { "width", "6.5" }, { NULL, NULL } };
        if (k26geo_resolve_road_width_m("primary", kv_lookup, with_width) != 6.5) {
            fprintf(stderr, "FAIL explicit width override\n"); fail++;
        }
    }
    {
        KV with_lanes[] = { { "lanes", "4" }, { NULL, NULL } };
        double w = k26geo_resolve_road_width_m("primary", kv_lookup, with_lanes);
        if (fabs(w - 14.0) > 1e-9) {
            fprintf(stderr, "FAIL lanes*3.5: got %g\n", w); fail++;
        }
    }
    {
        /* width wins over lanes */
        KV both[] = { { "width", "8" }, { "lanes", "10" }, { NULL, NULL } };
        if (k26geo_resolve_road_width_m("primary", kv_lookup, both) != 8.0) {
            fprintf(stderr, "FAIL width takes precedence over lanes\n"); fail++;
        }
    }

    if (fail) {
        fprintf(stderr, "%d failure(s)\n", fail);
        return 1;
    }
    puts("test_osmtags: OK");
    return 0;
}
