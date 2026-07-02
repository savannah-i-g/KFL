/* tests/astro/snapshot_roundtrip.c - snapshot save/load round-trip gate.
 *
 * Build a 10-body world, save to a snapshot, load it back, then
 * step both worlds 1000 times. The two worlds must remain
 * byte-identical (in PORTABLE mode) through the replay. */
#include "k26astro_rt/world.h"
#include "k26astro_rt/snapshot.h"
#include "k26astro_body/body.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static void add_sun_(K26AstroWorld *w)
{
    K26AstroBody b; k26astro_body_init(&b);
    strncpy(b.name, "sun", sizeof b.name - 1);
    b.mass   = 1.989e30;
    b.gm     = 1.32712440018e20;
    b.radius = 6.957e8;
    (void)k26astro_world_add_body(w, b);
}

static void add_orbiter_(K26AstroWorld *w, const char *name,
                          double r_m, double m_kg, double gm)
{
    K26AstroBody b; k26astro_body_init(&b);
    strncpy(b.name, name, sizeof b.name - 1);
    b.mass = m_kg;
    b.gm   = gm;
    b.radius = 1e6;
    K26AstroPos p = k26astro_pos_zero();
    p.lx = r_m;
    b.pos = p;
    /* Circular orbit speed around the Sun (m/s): v = sqrt(GM_sun / r). */
    double v = sqrt(1.32712440018e20 / r_m);
    b.vel.x = 0.0; b.vel.y = v; b.vel.z = 0.0;
    (void)k26astro_world_add_body(w, b);
}

static int bodies_equal_(const K26AstroBody *a, const K26AstroBody *b)
{
    if (strcmp(a->name, b->name) != 0) return 0;
    if (a->pos.sx != b->pos.sx || a->pos.sy != b->pos.sy || a->pos.sz != b->pos.sz) return 0;
    if (a->pos.lx != b->pos.lx || a->pos.ly != b->pos.ly || a->pos.lz != b->pos.lz) return 0;
    if (a->vel.x != b->vel.x || a->vel.y != b->vel.y || a->vel.z != b->vel.z)       return 0;
    if (a->mass != b->mass || a->gm != b->gm || a->radius != b->radius) return 0;
    return 1;
}

int main(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_PORTABLE,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);

    add_sun_(w);
    add_orbiter_(w, "mercury", 5.79e10,  3.302e23,    2.2032e13);
    add_orbiter_(w, "venus",   1.082e11, 4.8685e24,   3.24859e14);
    add_orbiter_(w, "earth",   1.496e11, 5.972e24,    3.986004418e14);
    add_orbiter_(w, "mars",    2.279e11, 6.4171e23,   4.282837e13);
    add_orbiter_(w, "jupiter", 7.785e11, 1.898e27,    1.26686534e17);
    add_orbiter_(w, "saturn",  1.434e12, 5.683e26,    3.7931187e16);
    add_orbiter_(w, "uranus",  2.871e12, 8.681e25,    5.793939e15);
    add_orbiter_(w, "neptune", 4.495e12, 1.024e26,    6.836529e15);
    add_orbiter_(w, "pluto",   5.906e12, 1.30900e22,  8.71e11);

    assert(k26astro_world_body_count(w) == 10);

    const char *path = "/tmp/k26astro_snapshot_test.bin";
    int rc = k26astro_world_snapshot_save(w, path);
    assert(rc == K26ASTRO_RT_OK);

    K26AstroWorld *w2 = k26astro_world_snapshot_load(path);
    assert(w2);
    assert(k26astro_world_body_count(w2) == 10);

    /* Byte-equality check pre-step. */
    for (int i = 0; i < 10; i++) {
        K26AstroGravState *a = k26astro_world_grav(w);
        K26AstroGravState *b = k26astro_world_grav(w2);
        assert(bodies_equal_(&a->bodies[i], &b->bodies[i]));
    }

    /* Step both worlds 1000 times. They should remain identical. */
    for (int s = 0; s < 1000; s++) {
        rc = k26astro_world_step(w,  3600.0);  assert(rc == K26ASTRO_RT_OK);
        rc = k26astro_world_step(w2, 3600.0);  assert(rc == K26ASTRO_RT_OK);
    }

    K26AstroGravState *a = k26astro_world_grav(w);
    K26AstroGravState *b = k26astro_world_grav(w2);
    for (int i = 0; i < 10; i++) {
        if (!bodies_equal_(&a->bodies[i], &b->bodies[i])) {
            fprintf(stderr,
                "snapshot_roundtrip: body %d diverged after 1000 steps\n", i);
            assert(0);
        }
    }

    k26astro_world_destroy(w);
    k26astro_world_destroy(w2);

    fprintf(stderr, "snapshot_roundtrip: 10-body world, 1000 steps, bit-identical\n");
    return 0;
}
