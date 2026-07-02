#define _POSIX_C_SOURCE 199309L  /* clock_gettime */

/* tests/bench_parse.c — measure in-process parse-only cost.
 *
 * Used by Docs/KFL/13-lsp-scout.md to ground the "is whole-file parse
 * fast enough for didChange" claim. Run with one input path. */

#include "kflc.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now_ms(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1000.0 + (double)t.tv_nsec / 1.0e6;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: bench_parse <file.kfl>\n");
        return 2;
    }
    /* Warmup: 100 iterations to populate page cache. */
    for (int i = 0; i < 100; i++) {
        KflcArena *a = kflc_arena_create();
        KflcDiag d; kflc_diag_init(&d, argv[1], stderr);
        kflc_parse_file(argv[1], a, &d);
        kflc_arena_release(a);
    }

    const int N = 10000;
    double t0 = now_ms();
    for (int i = 0; i < N; i++) {
        KflcArena *a = kflc_arena_create();
        KflcDiag d; kflc_diag_init(&d, argv[1], stderr);
        kflc_parse_file(argv[1], a, &d);
        kflc_arena_release(a);
    }
    double t1 = now_ms();

    double per = (t1 - t0) / (double)N;
    printf("%s: %d iterations, %.3f ms total, %.4f ms / parse\n",
           argv[1], N, t1 - t0, per);
    return 0;
}
