/* test_referenced_replay.c - referenced-mode determinism gate.
 *
 * Runs the same simple 3-body world twice in REFERENCED mode with
 * the same initial conditions. Each run writes to a separate log
 * file. The standalone k26astro_replay tool is invoked with
 * --check to compare the two logs byte-for-byte.
 *
 * The two runs share the same K26AstroWorld bookkeeping path, so a
 * byte-identical log is the floor: any divergence implies the
 * step-machinery is observing non-deterministic state somewhere
 * (e.g. wall clock, RNG, heap address leakage into a hashed log).
 *
 * SKIP exit 77 if the standalone tool isn't built; it lives in
 * tools/k26astro_replay and is only built by `make
 * -C tools/k26astro_replay`. */
#include "k26astro_rt/world.h"
#include "k26astro_rt/referenced.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define EXIT_SKIP 77

static int run_sim_(const char *log_path)
{
    K26AstroWorld *w = k26astro_world_create(
        K26ASTRO_MODE_REFERENCED, K26ASTRO_COORDS_Q64_64);
    if (!w) return -1;
    if (k26astro_world_set_ref_log_path(w, log_path) != 0) {
        k26astro_world_destroy(w);
        return -1;
    }

    /* Simple 3-body: Sun + 2 planets at 1 AU + 1.5 AU. */
    const double GM_SUN  = K26A_GM_SUN;
    const double R1      = K26A_AU_M;
    const double R2      = 1.5 * K26A_AU_M;
    const double V1      = sqrt(GM_SUN / R1);
    const double V2      = sqrt(GM_SUN / R2);

    K26AstroBody b;
    memset(&b, 0, sizeof b);
    b.kind = K26ASTRO_BODY_STAR;
    b.gm   = GM_SUN;
    b.mass = GM_SUN / 6.67430e-11;
    b.pos  = k26astro_pos_zero();
    b.vel  = (K26V3){ 0.0, 0.0, 0.0 };
    b.parent_body_idx = -1;
    snprintf(b.name, sizeof b.name, "%s", "sun");
    (void)k26astro_world_add_body(w, b);

    memset(&b, 0, sizeof b);
    b.kind = K26ASTRO_BODY_PLANET;
    b.gm   = 3.986004418e14;
    b.mass = 5.972e24;
    b.pos  = k26astro_pos_from_m(R1, 0.0, 0.0);
    b.vel  = (K26V3){ 0.0, V1, 0.0 };
    b.parent_body_idx = 0;
    snprintf(b.name, sizeof b.name, "%s", "earth");
    (void)k26astro_world_add_body(w, b);

    memset(&b, 0, sizeof b);
    b.kind = K26ASTRO_BODY_PLANET;
    b.gm   = 4.282837e13;
    b.mass = 6.4171e23;
    b.pos  = k26astro_pos_from_m(R2, 0.0, 0.0);
    b.vel  = (K26V3){ 0.0, V2, 0.0 };
    b.parent_body_idx = 0;
    snprintf(b.name, sizeof b.name, "%s", "mars");
    (void)k26astro_world_add_body(w, b);

    /* Step a few times. */
    double dt = 86400.0;   /* 1 day */
    for (int i = 0; i < 8; i++) {
        if (k26astro_world_step(w, dt) != K26ASTRO_RT_OK) {
            k26astro_world_destroy(w);
            return -1;
        }
    }

    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    /* Replay tool path: relative to tests/astro cwd. */
    const char *replay = "../../tools/k26astro_replay/k26astro_replay";
    struct stat st;
    if (stat(replay, &st) != 0) {
        fprintf(stderr,
                "test_referenced_replay: SKIP - replay tool not built "
                "(make -C tools/k26astro_replay)\n");
        return EXIT_SKIP;
    }

    const char *log_a = "/tmp/k26astro_referenced_a.log";
    const char *log_b = "/tmp/k26astro_referenced_b.log";
    unlink(log_a);
    unlink(log_b);

    if (run_sim_(log_a) != 0) {
        fprintf(stderr, "test_referenced_replay: run_a failed\n");
        return 1;
    }
    if (run_sim_(log_b) != 0) {
        fprintf(stderr, "test_referenced_replay: run_b failed\n");
        return 1;
    }

    /* Invoke the standalone --check tool. */
    char cmd[1024];
    snprintf(cmd, sizeof cmd, "%s --check %s %s", replay, log_a, log_b);
    int rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr,
                "test_referenced_replay: FAIL - replay --check returned "
                "%d (logs at %s and %s)\n", rc, log_a, log_b);
        return 1;
    }

    /* Sanity: also verify the logs aren't trivially empty. */
    struct stat sa, sb;
    if (stat(log_a, &sa) != 0 || stat(log_b, &sb) != 0) {
        fprintf(stderr, "test_referenced_replay: log stat failed\n");
        return 1;
    }
    if (sa.st_size < 64) {
        fprintf(stderr,
                "test_referenced_replay: FAIL - log too small (%lld bytes)\n",
                (long long)sa.st_size);
        return 1;
    }
    if (sa.st_size != sb.st_size) {
        fprintf(stderr,
                "test_referenced_replay: FAIL - log sizes differ "
                "(%lld vs %lld)\n",
                (long long)sa.st_size, (long long)sb.st_size);
        return 1;
    }

    printf("test_referenced_replay: PASS (%lld bytes, byte-identical "
           "across 2 REFERENCED-mode runs)\n", (long long)sa.st_size);
    unlink(log_a);
    unlink(log_b);
    return 0;
}
