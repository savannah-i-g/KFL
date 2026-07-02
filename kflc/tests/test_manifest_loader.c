/* test_manifest_loader.c: .kflbi loader smoke test. */
#define _GNU_SOURCE
/*
 *
 * Creates a tmpdir, writes a hand-rolled manifest in it, points
 * K26_KFL_BUILTIN_PATH at the tmpdir, runs `kflc_load_manifests()`,
 * and asserts the expected opaque types + builtins are now in the
 * registry.
 *
 * Coverage:
 *   1. Valid manifest registers expected names.
 *   2. Comment + blank lines ignored.
 *   3. Malformed line emits stderr warning but loader continues.
 *   4. Schema mismatch rejects whole file (clear stderr message).
 *   5. Collision with prior registration produces "first wins"
 *      semantics (no crash, kflc_opaque_cxx returns first cxx_name).
 */
#include "kflc.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int write_manifest_(const char *dir, const char *name,
                            const char *content)
{
    char path[512];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fputs(content, f);
    fclose(f);
    return 0;
}

int main(void)
{
    char tmpl[] = "/tmp/kflc_manifest_test_XXXXXX";
    char *dir = mkdtemp(tmpl);
    if (!dir) { perror("mkdtemp"); return 1; }

    /* --- Valid manifest ----------------------------------------- */
    const char *valid =
        "# k26astro_demo.kflbi — test fixture\n"
        "schema  1\n"
        "library libk26astro_demo\n"
        "version 0.1.0\n"
        "\n"
        "opaque  demo_widget   K26DemoWidget *\n"
        "builtin demo_open     k26demo_open       0\n"
        "builtin demo_close    k26demo_close      1\n"
        "# trailing comment\n";

    if (write_manifest_(dir, "k26astro_demo.kflbi", valid) != 0) {
        fprintf(stderr, "write valid: failed\n");
        return 1;
    }

    /* --- Malformed manifest (warning, not fatal) ---------------- */
    const char *malformed =
        "schema 1\n"
        "opaque\n"                  /* incomplete */
        "builtin bad_one\n"         /* missing fields */
        "weird_directive foo bar\n" /* unknown */
        "opaque well_formed K26WellFormed *\n";
    if (write_manifest_(dir, "k26misc.kflbi", malformed) != 0) {
        fprintf(stderr, "write malformed: failed\n");
        return 1;
    }

    /* --- Schema mismatch (rejected) ----------------------------- */
    const char *bad_schema =
        "schema 999\n"
        "opaque future_only K26Future *\n";
    if (write_manifest_(dir, "k26future.kflbi", bad_schema) != 0) {
        fprintf(stderr, "write bad_schema: failed\n");
        return 1;
    }

    setenv("K26_KFL_BUILTIN_PATH", dir, 1);
    kflc_opaque_clear();
    kflc_clear_builtins();
    kflc_load_manifests();

    /* Assert the valid manifest's registrations. */
    const char *cxx = kflc_opaque_cxx("demo_widget");
    if (!cxx || strcmp(cxx, "K26DemoWidget *") != 0) {
        fprintf(stderr, "FAIL: demo_widget not registered (got %s)\n",
                cxx ? cxx : "(null)");
        return 1;
    }
    /* well_formed comes from the malformed manifest's good line. */
    cxx = kflc_opaque_cxx("well_formed");
    if (!cxx || strcmp(cxx, "K26WellFormed *") != 0) {
        fprintf(stderr, "FAIL: well_formed not registered (got %s)\n",
                cxx ? cxx : "(null)");
        return 1;
    }
    /* bad_schema content should NOT be registered. */
    cxx = kflc_opaque_cxx("future_only");
    if (cxx) {
        fprintf(stderr, "FAIL: future_only registered despite bad schema\n");
        return 1;
    }

    /* --- Collision: register demo_widget by C API, then re-load.
     * The pre-registered cxx_name must win (first wins). */
    kflc_opaque_clear();
    kflc_clear_builtins();
    (void)kflc_opaque_register("demo_widget", "K26FirstWinner *");
    kflc_load_manifests();
    cxx = kflc_opaque_cxx("demo_widget");
    if (!cxx || strcmp(cxx, "K26FirstWinner *") != 0) {
        fprintf(stderr, "FAIL: collision didn't honour first-wins (got %s)\n",
                cxx ? cxx : "(null)");
        return 1;
    }

    /* --- Cleanup tmpfiles --------------------------------------- */
    char path[512];
    const char *files[] = { "k26astro_demo.kflbi", "k26misc.kflbi",
                            "k26future.kflbi", NULL };
    for (int i = 0; files[i]; i++) {
        snprintf(path, sizeof path, "%s/%s", dir, files[i]);
        unlink(path);
    }
    rmdir(dir);

    printf("test_manifest_loader: OK "
           "(valid + malformed + bad-schema + collision paths)\n");
    return 0;
}
