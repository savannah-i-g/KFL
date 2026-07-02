/* test_body_binding_emit.c: parse-time body bindings.
 *
 * Compiles kflc/examples/astro_body_parse_time.kfl via `kflc --emit`
 * (captures the emitted C to a tmpfile), greps for the expected
 * patterns:
 *
 *   PRESENT:   `int _kfl_body_earth_idx = -1;`        (prologue)
 *              `_kfl_body_earth_idx = k26astro_world_add_body`
 *              `int _kfl_idx = _kfl_body_earth_idx;`  (propagate)
 *              `int _kfl_t = _kfl_body_earth_idx;`    (observe target)
 *
 *   ABSENT:    `k26astro_world_find_body(world, "earth"`  for known names
 *
 * Unknown names (e.g. via runtime input strings) still resolve via
 * find_body — but the parse-time-known fixture has no such names,
 * so the emitted C should have zero find_body calls in the
 * relevant fn body. */
#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *slurp_(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    if (sz < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char *buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return NULL; }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf); fclose(f); return NULL;
    }
    buf[sz] = '\0';
    fclose(f);
    return buf;
}

int main(void)
{
    /* Run kflc --emit against the fixture and capture the output. */
    const char *cmd =
        "./bin/kflc --emit examples/astro_body_parse_time.kfl "
        "> /tmp/test_body_binding_emit.cc 2>&1";
    int rc = system(cmd);
    assert(rc == 0);

    char *src = slurp_("/tmp/test_body_binding_emit.cc");
    assert(src != NULL);

    /* PRESENT — prologue idx decls + assignments + idx-deref uses. */
    assert(strstr(src, "int _kfl_body_sun_idx = -1;") != NULL);
    assert(strstr(src, "int _kfl_body_earth_idx = -1;") != NULL);
    assert(strstr(src, "int _kfl_body_mars_idx = -1;") != NULL);

    assert(strstr(src, "_kfl_body_sun_idx = k26astro_world_add_body") != NULL);
    assert(strstr(src, "_kfl_body_earth_idx = k26astro_world_add_body") != NULL);
    assert(strstr(src, "_kfl_body_mars_idx = k26astro_world_add_body") != NULL);

    /* propagate earth + propagate mars → both deref via idx. */
    assert(strstr(src, "int _kfl_idx = _kfl_body_earth_idx;") != NULL);
    assert(strstr(src, "int _kfl_idx = _kfl_body_mars_idx;") != NULL);

    /* observe earth from sun + observe mars from earth → both endpoints
     * resolve via idx. */
    assert(strstr(src, "int _kfl_t = _kfl_body_earth_idx;") != NULL);
    assert(strstr(src, "int _kfl_o = _kfl_body_sun_idx;") != NULL);
    assert(strstr(src, "int _kfl_t = _kfl_body_mars_idx;") != NULL);
    assert(strstr(src, "int _kfl_o = _kfl_body_earth_idx;") != NULL);

    /* parent=sun for earth + mars → resolved via idx. */
    assert(strstr(src, "_kfl_b.parent_body_idx = _kfl_body_sun_idx;") != NULL);

    /* ABSENT — no runtime find_body for parse-time-known names. */
    assert(strstr(src, "k26astro_world_find_body(world, \"sun\")") == NULL);
    assert(strstr(src, "k26astro_world_find_body(world, \"earth\")") == NULL);
    assert(strstr(src, "k26astro_world_find_body(world, \"mars\")") == NULL);

    free(src);
    printf("test_body_binding_emit: OK\n");
    return 0;
}
