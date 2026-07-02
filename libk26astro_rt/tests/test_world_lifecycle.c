/* test_world_lifecycle.c — create / add bodies / destroy. Verifies
 * FPU state is restored across lifecycle. */
#include "k26astro_rt/world.h"
#include "k26astro_body/body.h"

#include <stdio.h>
#include <assert.h>
#include <fenv.h>
#include <string.h>

int main(void)
{
    /* Save the caller's FPU mode (intentionally NOT pinned). */
    fesetround(FE_DOWNWARD);
    int saved = fegetround();

    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_PORTABLE,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    /* While world is alive, FPU should be pinned to FE_TONEAREST. */
    assert(fegetround() == FE_TONEAREST);

    K26AstroBody b;
    k26astro_body_init(&b);
    strncpy(b.name, "sun", sizeof b.name - 1);
    b.mass   = 1.989e30;
    b.gm     = 1.32712440018e20;
    b.radius = 6.957e8;

    int idx = k26astro_world_add_body(w, b);
    assert(idx == 0);
    assert(k26astro_world_body_count(w) == 1);

    int found = k26astro_world_find_body(w, "sun");
    assert(found == 0);

    int missing = k26astro_world_find_body(w, "nope");
    assert(missing == -1);

    k26astro_world_destroy(w);

    /* After destroy, the caller's FPU mode should be restored. */
    assert(fegetround() == saved);

    /* Restore to a sane default for any test runner after us. */
    fesetround(FE_TONEAREST);

    fprintf(stderr, "test_world_lifecycle: OK\n");
    return 0;
}
