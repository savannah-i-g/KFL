#define _GNU_SOURCE  /* mkstemps + fdopen */

/* tests/round_trip.c
 *
 * Round-trip property test for KFL forms. Promised by
 * Docs/KFL/01-overview.md §3.3: `parse / serialize / parse` yields a
 * semantically identical AST.
 *
 * For each input .kfl file:
 *   1. Parse to AST_A (with diag_A).
 *   2. Serialize AST_A to a temp file.
 *   3. Parse the temp file to AST_B (with diag_B).
 *   4. Assert kflc_ast_equal(AST_A, AST_B) and that diagnostic counts
 *      match.
 *
 * Exits 0 if every input round-trips cleanly, 1 otherwise. The first
 * failing example halts the run with a diagnostic identifying which
 * check failed.
 *
 * Wired into `kflc/Makefile` as the `test` target.
 */

#include "kflc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int run_one(const char *path)
{
    /* Pass 1: parse the original. */
    KflcArena *arena_a = kflc_arena_create();
    KflcDiag   diag_a;
    kflc_diag_init(&diag_a, path, stderr);
    KflcNode  *form_a  = kflc_parse_file(path, arena_a, &diag_a);
    if (!form_a) {
        fprintf(stderr, "round_trip: %s: pass-1 parse returned NULL "
                "(%d error(s))\n", path, diag_a.errors);
        kflc_arena_release(arena_a);
        return 1;
    }

    /* Serialize to a temp file. mkstemp + fdopen so the file gets
     * cleaned up via unlink even if a later step fails. */
    char tmp_path[] = "/tmp/kflc_roundtrip_XXXXXX.kfl";
    int  tmp_fd = mkstemps(tmp_path, 4);
    if (tmp_fd < 0) {
        perror("round_trip: mkstemps");
        kflc_arena_release(arena_a);
        return 1;
    }
    FILE *tmp_out = fdopen(tmp_fd, "w");
    if (!tmp_out) {
        perror("round_trip: fdopen");
        close(tmp_fd);
        unlink(tmp_path);
        kflc_arena_release(arena_a);
        return 1;
    }
    if (kflc_serialize(tmp_out, form_a) != 0) {
        fprintf(stderr, "round_trip: %s: serialize failed\n", path);
        fclose(tmp_out);
        unlink(tmp_path);
        kflc_arena_release(arena_a);
        return 1;
    }
    fclose(tmp_out);

    /* Pass 2: parse the serialized form. */
    KflcArena *arena_b = kflc_arena_create();
    KflcDiag   diag_b;
    kflc_diag_init(&diag_b, tmp_path, stderr);
    KflcNode  *form_b  = kflc_parse_file(tmp_path, arena_b, &diag_b);

    int rc = 0;

    if (!form_b) {
        fprintf(stderr, "round_trip: %s: pass-2 parse returned NULL "
                "(serialized form unparseable: %d error(s))\n"
                "  preserved at %s for inspection\n",
                path, diag_b.errors, tmp_path);
        rc = 1;
        goto cleanup_b;  /* don't unlink — leave for debugging */
    }

    if (diag_a.errors != diag_b.errors) {
        fprintf(stderr, "round_trip: %s: error-count drift "
                "(pass-1=%d, pass-2=%d)\n"
                "  preserved at %s for inspection\n",
                path, diag_a.errors, diag_b.errors, tmp_path);
        rc = 1;
        goto cleanup_b;
    }
    if (diag_a.warnings != diag_b.warnings) {
        fprintf(stderr, "round_trip: %s: warning-count drift "
                "(pass-1=%d, pass-2=%d)\n"
                "  preserved at %s for inspection\n",
                path, diag_a.warnings, diag_b.warnings, tmp_path);
        rc = 1;
        goto cleanup_b;
    }
    if (!kflc_ast_equal(form_a, form_b)) {
        fprintf(stderr, "round_trip: %s: AST mismatch after round-trip\n"
                "  serialized form preserved at %s for inspection\n",
                path, tmp_path);
        rc = 1;
        goto cleanup_b;
    }

    /* Clean. */
    unlink(tmp_path);

cleanup_b:
    kflc_arena_release(arena_b);
    kflc_arena_release(arena_a);
    return rc;
}

extern void kflc_register_astro_builtins(void);

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: round_trip <file.kfl> [<file.kfl> ...]\n");
        return 2;
    }

    /* Register astro opaques (body, world, starfield, ...) so
     * fixtures using opaque-typed bindings can parse. The kflc binary
     * does this from main; the round_trip harness needs the same. */
    kflc_register_astro_builtins();

    int n_ok   = 0;
    int n_fail = 0;
    for (int i = 1; i < argc; i++) {
        int rc = run_one(argv[i]);
        if (rc == 0) {
            n_ok++;
            printf("PASS  %s\n", argv[i]);
        } else {
            n_fail++;
            printf("FAIL  %s\n", argv[i]);
        }
    }

    printf("\nround_trip: %d passed, %d failed\n", n_ok, n_fail);
    return n_fail == 0 ? 0 : 1;
}
