/* test_w8_arena_emit.c: KFL memory model emit gate.
 *
 * Compiles kflc/examples/v4_arena.kfl via `kflc --emit` (captures the
 * emitted C++ to a tmpfile), greps for the expected patterns:
 *
 *   PRESENT (when the form declares any arena):
 *     `#include "kflc_runtime.h"`             (conditional include)
 *     `static K26KflArena *_kfl_arena_world_arena = NULL;`
 *     `static K26KflArena *_kfl_arena_scratch = NULL;`
 *     `static K26KflArena *_kfl_arena_small_buffer = NULL;`
 *     `static void _kfl_arena_reset_cleanup_(K26KflArena **a)`
 *     `static void _kfl_init_arenas(void) {`  (one-shot factory)
 *     `_kfl_arena_world_arena = k26kfl_arena_create((size_t)67108864L);`
 *     `_kfl_arena_scratch = k26kfl_arena_create((size_t)16384L);`
 *     `_kfl_arena_small_buffer = k26kfl_arena_create((size_t)4096L);`
 *     `_kfl_init_arenas();`                   (form ctor call)
 *
 *     fn-body `allocator = X` (no reset_mode, default fn_exit) emits
 *     with cleanup attribute for auto-reset on every fn return:
 *       `K26KflArena *_kfl_active_arena `
 *       `__attribute__((cleanup(_kfl_arena_reset_cleanup_))) = ...`
 *
 *   ABSENT (forms without arena decls keep the slim include surface):
 *     A non-arena example does NOT pull in kflc_runtime.h.
 *
 * Compares against the 64 MB / 16 KB / 4096 unit-suffix decoding
 * to confirm parse and round-trip preserve byte counts exactly.
 */
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
    /* Run kflc --emit against the arena fixture and capture the output. */
    int rc = system("./bin/kflc --emit examples/v4_arena.kfl "
                    "> /tmp/test_w8_arena_emit.cc 2>&1");
    assert(rc == 0);

    char *v4 = slurp_("/tmp/test_w8_arena_emit.cc");
    assert(v4 != NULL);

    /* PRESENT: runtime header included exactly once. */
    assert(strstr(v4, "#include \"kflc_runtime.h\"") != NULL);

    /* Static handles for every arena declaration. */
    assert(strstr(v4, "static K26KflArena *_kfl_arena_world_arena = NULL;")
           != NULL);
    assert(strstr(v4, "static K26KflArena *_kfl_arena_scratch = NULL;")
           != NULL);
    assert(strstr(v4, "static K26KflArena *_kfl_arena_small_buffer = NULL;")
           != NULL);

    /* Cleanup helper emitted alongside the arena decls. */
    assert(strstr(v4,
        "static void _kfl_arena_reset_cleanup_(K26KflArena **a)") != NULL);

    /* One-shot init factory. */
    assert(strstr(v4, "static void _kfl_init_arenas(void) {") != NULL);

    /* Capacity-decode check: confirms parser handled KB/MB suffixes
     * and raw bytes correctly:
     *   64 MB → 64 * 1024 * 1024 = 67108864
     *   16 KB → 16 * 1024        = 16384
     *   4096  → 4096
     */
    assert(strstr(v4,
        "_kfl_arena_world_arena = k26kfl_arena_create((size_t)67108864L);")
           != NULL);
    assert(strstr(v4,
        "_kfl_arena_scratch = k26kfl_arena_create((size_t)16384L);")
           != NULL);
    assert(strstr(v4,
        "_kfl_arena_small_buffer = k26kfl_arena_create((size_t)4096L);")
           != NULL);

    /* Form ctor calls the init. */
    assert(strstr(v4, "_kfl_init_arenas();") != NULL);

    /* fn-body `allocator = world_arena` installs the active-arena alias.
     * Default reset_mode (fn_exit) emits with cleanup attr. */
    assert(strstr(v4,
        "K26KflArena *_kfl_active_arena "
        "__attribute__((cleanup(_kfl_arena_reset_cleanup_))) "
        "= _kfl_arena_world_arena;") != NULL);

    free(v4);

    /* ABSENT: a non-arena example should not pull in kflc_runtime.h.
     * Run kflc --emit on a fixture without arena decls and confirm
     * the slim include surface. */
    rc = system("./bin/kflc --emit examples/astro_apophis_2029.kfl "
                "> /tmp/test_w8_arena_emit_v3.cc 2>&1");
    assert(rc == 0);

    char *v3 = slurp_("/tmp/test_w8_arena_emit_v3.cc");
    assert(v3 != NULL);
    assert(strstr(v3, "kflc_runtime.h") == NULL);
    /* And no arena artifacts. */
    assert(strstr(v3, "_kfl_arena_") == NULL);
    assert(strstr(v3, "_kfl_init_arenas") == NULL);
    free(v3);

    printf("test_w8_arena_emit: OK\n");
    return 0;
}
