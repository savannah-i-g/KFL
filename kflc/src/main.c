/* kflc driver. Parses a .kfl file, emits C++, and invokes the system
 * C++ compiler to produce a native binary.
 *
 * Environment variables (override compile flags):
 *   CXX           system C++ compiler         default: c++
 *   KFLC_CFLAGS   compile flags               default: -O2 -g -std=c++11 -Wno-format-truncation
 *   KFLC_LDLIBS   link libraries              default: -lm
 *
 * The default link line is minimal; programs that consume KFL_Stack
 * libraries (libk26astro_*, libk26m3d, libk26plot, etc.) should set
 * KFLC_LDLIBS to add the matching `-l` flags.
 *
 * Handler source files are resolved relative to the .kfl file's
 * directory and prepended to the cc command's source list.
 */

#include "kflc.h"

#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

static int usage(void)
{
    fputs(
        "usage: kflc <input.kfl> [-o <output>] [-c] [--dump] [--emit] [--check]\n"
        "\n"
        "  --dump          parse and dump the AST to stdout\n"
        "  --emit          emit C++ source to stdout (no compile)\n"
        "  --check         parse + run emit silently; exit 0 on clean,\n"
        "                  1 on parse or emit errors. Nothing written to\n"
        "                  stdout. Intended for editor save-hooks, CI lint\n"
        "                  stages, and pre-commit checks.\n"
        "  -o <output>     produce a compiled binary at <output>\n"
        "  -c              also keep the emitted C++ source alongside <output>\n",
        stderr);
    return 2;
}

int main(int argc, char **argv)
{
    const char *input  = NULL;
    const char *output = NULL;
    int keep_cxx = 0;
    int dump  = 0;
    int emit  = 0;
    int check = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output = argv[++i];
        } else if (strcmp(argv[i], "-c") == 0) {
            keep_cxx = 1;
        } else if (strcmp(argv[i], "--dump") == 0) {
            dump = 1;
        } else if (strcmp(argv[i], "--emit") == 0) {
            emit = 1;
        } else if (strcmp(argv[i], "--check") == 0) {
            check = 1;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            return usage();
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "kflc: unknown option `%s`\n", argv[i]);
            return usage();
        } else if (!input) {
            input = argv[i];
        } else {
            fprintf(stderr, "kflc: unexpected argument `%s`\n", argv[i]);
            return usage();
        }
    }

    if (!input) return usage();
    /* `--check` is its own terminal mode; no implicit fall-through to
     * dump when nothing else is requested. */
    if (!output && !dump && !emit && !check) dump = 1;

    KflcArena *arena = kflc_arena_create();
    KflcDiag   diag;
    kflc_diag_init(&diag, input, stderr);

    /* Load builtin manifests from K26_KFL_BUILTIN_PATH (default
     * /usr/share/kflc/builtins/) and register the astro opaque
     * types + builtins. The hardcoded astro registrar runs after
     * as a fallback for embedders that don't ship the manifest
     * path (host build, embedded use); its collision-detection
     * makes the duplicate registrations no-op when both paths
     * populate the same names. */
    kflc_load_manifests();
    kflc_register_astro_builtins();

    KflcNode *form = kflc_parse_file(input, arena, &diag);
    if (!form || diag.errors) {
        fprintf(stderr, "kflc: parse failed (%d error%s)\n",
                diag.errors, diag.errors == 1 ? "" : "s");
        kflc_arena_release(arena);
        return 1;
    }

    /* --check: silently run emit through /dev/null so we surface
     * emit-phase diagnostics (unknown builtin, arity mismatch,
     * type errors in `compute expression`) that --dump skips. Exits
     * with the parse+emit error count as a 0/1 status. Nothing else
     * runs. */
    if (check) {
        FILE *devnull = fopen("/dev/null", "w");
        if (devnull) {
            (void)kflc_emit_cxx(devnull, form, &diag);
            fclose(devnull);
        } else {
            /* /dev/null missing is so weird we should report it
             * rather than silently skip the emit pass. */
            fprintf(stderr, "kflc: --check: cannot open /dev/null\n");
            kflc_arena_release(arena);
            return 1;
        }
        int errs = diag.errors;
        kflc_arena_release(arena);
        return errs ? 1 : 0;
    }

    if (dump) {
        kflc_dump_node(stdout, form, 0);
    }

    if (emit) {
        if (kflc_emit_cxx(stdout, form, &diag) != 0) {
            fprintf(stderr, "kflc: emit failed (%d error%s)\n",
                    diag.errors, diag.errors == 1 ? "" : "s");
            kflc_arena_release(arena);
            return 1;
        }
    }

    if (output) {
        /* Pick a path for the emitted C++. */
        char cxx_path[512];
        if (keep_cxx) {
            snprintf(cxx_path, sizeof cxx_path, "%s.kflc.cc", output);
        } else {
            snprintf(cxx_path, sizeof cxx_path, "/tmp/kflc-%d.cc", (int)getpid());
        }

        FILE *cf = fopen(cxx_path, "w");
        if (!cf) {
            fprintf(stderr, "kflc: cannot open `%s` for write\n", cxx_path);
            kflc_arena_release(arena);
            return 1;
        }
        if (kflc_emit_cxx(cf, form, &diag) != 0) {
            fprintf(stderr, "kflc: emit failed (%d error%s)\n",
                    diag.errors, diag.errors == 1 ? "" : "s");
            fclose(cf);
            if (!keep_cxx) unlink(cxx_path);
            kflc_arena_release(arena);
            return 1;
        }
        fclose(cf);

        /* Build the cc command. */
        const char *cxx    = getenv("CXX");
        if (!cxx) cxx = "c++";
        const char *cflags = getenv("KFLC_CFLAGS");
        if (!cflags) cflags = "-O2 -g -std=c++11 -Wno-format-truncation";
        const char *ldlibs = getenv("KFLC_LDLIBS");
        if (!ldlibs) ldlibs = "-lm";

        char cmd[8192];
        int  n = snprintf(cmd, sizeof cmd,
            "%s %s -o %s %s", cxx, cflags, output, cxx_path);
        n += snprintf(cmd + n, sizeof cmd - n, " %s", ldlibs);

        if ((size_t)n >= sizeof cmd) {
            fprintf(stderr, "kflc: command line too long\n");
            if (!keep_cxx) unlink(cxx_path);
            kflc_arena_release(arena);
            return 1;
        }

        fprintf(stderr, "kflc: %s\n", cmd);
        int rc = system(cmd);
        if (!keep_cxx) unlink(cxx_path);
        kflc_arena_release(arena);
        if (rc != 0) {
            fprintf(stderr, "kflc: compile failed (rc=%d)\n", rc);
            return 1;
        }
        return 0;
    }

    kflc_arena_release(arena);
    return 0;
}
