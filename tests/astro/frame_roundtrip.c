/* tests/astro/frame_roundtrip.c - inertial-frame round-trip gate.
 *
 * For each pair of registered inertial frames (F1, F2), pick 1000
 * random positions tagged F1, transform to F2, transform back to F1,
 * assert max residual < 1e-7 m.
 *
 * The frame layer currently supports inertial<->inertial only;
 * body-fixed pairs return K26A_FRAME_E_BODY_NOT_REGISTERED until
 * body-attitude is wired into the transform. The inertial-only subset
 * is the gate. */
#include "k26astro_core/frame.h"
#include "k26astro_core/pos.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static K26AstroPos random_pos_(unsigned *seed)
{
    K26AstroPos p;
    p.sx = (int64_t)(rand_r(seed) % 11) - 5;
    p.sy = (int64_t)(rand_r(seed) % 11) - 5;
    p.sz = (int64_t)(rand_r(seed) % 11) - 5;
    double half = (double)(1ULL << 35); /* EDGE/2 */
    p.lx = ((double)rand_r(seed) / (double)RAND_MAX) * 2.0 * half - half;
    p.ly = ((double)rand_r(seed) / (double)RAND_MAX) * 2.0 * half - half;
    p.lz = ((double)rand_r(seed) / (double)RAND_MAX) * 2.0 * half - half;
    return p;
}

static double pos_diff_(const K26AstroPos *a, const K26AstroPos *b)
{
    K26V3 d = k26astro_pos_sub(a, b);
    return sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
}

int main(void)
{
    const K26AstroFrameId inertial[] = {
        K26A_FRAME_ICRF,
        K26A_FRAME_HCRF,
        K26A_FRAME_ECI,
        K26A_FRAME_GCRS
    };
    int n = (int)(sizeof inertial / sizeof inertial[0]);

    K26AstroEpoch t = k26astro_epoch_j2000_tt();

    unsigned seed = 0xC0FFEEu;
    double max_residual = 0.0;
    int passes = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            for (int k = 0; k < 1000; k++) {
                K26AstroPos p0 = random_pos_(&seed);
                k26astro_pos_normalise(&p0);

                K26AstroPosTagged in  = { p0, inertial[i] };
                K26AstroPosTagged mid;
                K26AstroPosTagged out;
                int rc1 = k26astro_frame_transform(&mid, &in,  inertial[j], &t);
                if (rc1 != K26A_FRAME_OK) continue;
                int rc2 = k26astro_frame_transform(&out, &mid, inertial[i], &t);
                if (rc2 != K26A_FRAME_OK) continue;

                double d = pos_diff_(&p0, &out.p);
                if (d > max_residual) max_residual = d;
                passes++;
            }
        }
    }

    fprintf(stderr,
            "frame_roundtrip: %d round-trip transforms, max residual = %.3e m\n",
            passes, max_residual);
    assert(passes > 0);
    assert(max_residual < 1e-7);
    return 0;
}
