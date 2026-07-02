#include <math.h>
#include <stdio.h>
#include <string.h>

#include "k26geo.h"

static int fail = 0;

#define EXPECT_EQ(got, want, name)                                    \
    do {                                                              \
        if ((got) != (want)) {                                        \
            fprintf(stderr, "FAIL %s: got %zu, want %zu\n",           \
                    name, (size_t)(got), (size_t)(want));             \
            fail++;                                                   \
        }                                                             \
    } while (0)

int main(void)
{
    /* Empty + singleton + 2-point */
    {
        EXPECT_EQ(k26geo_rdp_xz(NULL, 0, 0.5), 0, "null");
        double a[] = { 1, 2 };
        EXPECT_EQ(k26geo_rdp_xz(a, 1, 0.5), 1, "singleton");
        EXPECT_EQ(k26geo_rdp_xz(a, 2, 0.5), 2, "pair");
    }

    /* Collinear points collapse to endpoints. */
    {
        double xs[] = {
            0, 0,
            1, 0,
            2, 0,
            3, 0,
            4, 0
        };
        size_t n = k26geo_rdp_xz(xs, 5, 0.5);
        EXPECT_EQ(n, 2, "collinear collapse");
        if (xs[0] != 0.0 || xs[1] != 0.0 || xs[2] != 4.0 || xs[3] != 0.0) {
            fprintf(stderr, "FAIL collinear endpoints: %g,%g,%g,%g\n",
                    xs[0], xs[1], xs[2], xs[3]);
            fail++;
        }
    }

    /* Sharp bend survives. */
    {
        double xs[] = {
            0, 0,
            2, 0,
            2, 10,
            4, 10
        };
        size_t n = k26geo_rdp_xz(xs, 4, 0.5);
        EXPECT_EQ(n, 4, "sharp corner kept");
    }

    /* Subtle bend below epsilon collapses; above doesn't. */
    {
        double xs[] = {
            0, 0,
            5, 0.1,
            10, 0
        };
        double ys[6]; memcpy(ys, xs, sizeof xs);
        size_t n1 = k26geo_rdp_xz(ys, 3, 0.5);
        EXPECT_EQ(n1, 2, "subtle bend collapses at eps=0.5");

        memcpy(ys, xs, sizeof xs);
        size_t n2 = k26geo_rdp_xz(ys, 3, 0.01);
        EXPECT_EQ(n2, 3, "subtle bend survives at eps=0.01");
    }

    /* eps = 0 (LOSSLESS) preserves everything. */
    {
        double xs[] = {
            0, 0,
            1, 0,
            2, 0
        };
        size_t n = k26geo_rdp_xz(xs, 3, K26GEO_RDP_LOSSLESS_M);
        EXPECT_EQ(n, 3, "LOSSLESS preserves count");
    }

    /* Larger zigzag stress. */
    {
        double xs[200];
        for (size_t i = 0; i < 100; i++) {
            xs[2*i]     = (double)i;
            xs[2*i + 1] = (i % 2 == 0) ? 0.0 : 0.05;   /* below 0.5 m */
        }
        size_t n = k26geo_rdp_xz(xs, 100, K26GEO_RDP_HIGH_M);
        if (n >= 100) {
            fprintf(stderr, "FAIL HIGH eps should drop noise (got %zu)\n", n);
            fail++;
        }
    }

    if (fail) {
        fprintf(stderr, "%d failure(s)\n", fail);
        return 1;
    }
    puts("test_rdp: OK");
    return 0;
}
