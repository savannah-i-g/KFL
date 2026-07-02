/* body.c — K26AstroBody struct ops + major-body table.
 *
 * The major-body table is the authoritative source of mass / radius
 * / J2 for the canonical Solar System set. Constants:
 *
 *   - GM_sun / GM_planet      IAU 2015 Resolution B3 nominal values
 *                             (exact SI conversion factors); see
 *                             k26astro_core/consts.h K26A_GM_*.
 *   - Mass derived            mass = GM / G. G itself has CODATA
 *                             uncertainty; for gravity computations
 *                             always prefer the .gm field over .mass.
 *   - Radius                  JPL Solar System Dynamics planetary
 *                             data sheet (current as of 2024).
 *   - J2 (zonal harmonic)     IERS Conventions 2010 / planet-specific
 *                             gravity-field publications.
 *
 * Updates: when IAU revises the nominal constants (next plenary
 * resolution) or JPL refreshes the planetary data sheet, edit the
 * table entry in place. Downstream lookups (k26astro_body_load_major)
 * pick up the new value without any code change. */
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"

#include <stddef.h>
#include <string.h>
#include <strings.h>

void k26astro_body_init(K26AstroBody *b)
{
    if (!b) return;
    memset(b, 0, sizeof(*b));
    b->attitude = k26m3d_quat_identity();
    /* memset already zeroes everything else (kind, ephem_naif_id,
     * on_rails, etc.). parent_body_idx = 0 means "first body" by
     * default; callers that mean "no parent" must explicitly write
     * -1 — they typically do this via body_load_major + ephem
     * binding. */
    b->parent_body_idx = -1;
}

void k26astro_body_set_mass(K26AstroBody *b, double mass_kg)
{
    if (!b) return;
    if (mass_kg < 0.0) mass_kg = 0.0;
    b->mass = mass_kg;
    b->gm   = K26A_G * mass_kg;
}

int k26astro_body_find_by_name(const K26AstroBody *bodies, int n,
                                const char *name)
{
    if (!bodies || !name) return -1;
    for (int i = 0; i < n; i++) {
        if (strncasecmp(bodies[i].name, name, K26ASTRO_BODY_NAME_MAX) == 0) {
            return i;
        }
    }
    return -1;
}

/* ---- Major-body table -----------------------------------------
 *
 * Each entry is named so future revisions are a clean text edit. The
 * `rotation_model_name` field cross-references the IAU rotation table
 * in iau_rotations.c — the body picks up that table's coefficients at
 * world construction.
 *
 * Sources keyed in the per-entry comments. */
static const K26AstroMajorBody MAJOR_BODIES[] = {
    /* Sun.
     * GM: IAU 2015 R.B3 (exact). Radius: nominal equatorial.
     * J2: solar J2 ≈ 2e-7 (Mecheri et al. 2004); typically below
     * the threshold any K26 use cases need but recorded for
     * completeness. */
    {
        .name                 = "sun",
        .kind                 = K26ASTRO_BODY_STAR,
        .naif_id              = 10,
        .gm                   = K26A_GM_SUN,
        .mass                 = K26A_GM_SUN / K26A_G,
        .radius               = K26A_R_SUN_EQU,
        .polar_radius         = K26A_R_SUN_EQU,   /* assumed spherical */
        .j2                   = 2.0e-7,
        .rotation_model_name  = "iau2018:sun"
    },

    /* Mercury.
     * GM: DE441 (1.6601082e-7 * GM_sun) ≈ 2.2032e13.
     * Radius: 2 439.7 km equatorial.
     * J2: 5.03e-5 (MESSENGER mission, Mazarico et al. 2014). */
    {
        .name                 = "mercury",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 199,
        .gm                   = 2.2031780000e13,
        .mass                 = 2.2031780000e13 / K26A_G,
        .radius               = 2.4397e6,
        .polar_radius         = 2.4397e6,
        .j2                   = 5.03e-5,
        .rotation_model_name  = "iau2018:mercury"
    },

    /* Venus.
     * GM: DE441. Radius: 6 051.8 km. J2: ~4.4e-6 (Magellan + radar). */
    {
        .name                 = "venus",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 299,
        .gm                   = 3.2485859200e14,
        .mass                 = 3.2485859200e14 / K26A_G,
        .radius               = 6.0518e6,
        .polar_radius         = 6.0518e6,
        .j2                   = 4.458e-6,
        .rotation_model_name  = "iau2018:venus"
    },

    /* Earth.
     * GM: IAU 2015 R.B3 nominal. Radius: WGS84 equatorial.
     * Polar radius: WGS84. J2: EGM2008 leading zonal (1.0826267e-3). */
    {
        .name                 = "earth",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 399,
        .gm                   = K26A_GM_EARTH,
        .mass                 = K26A_GM_EARTH / K26A_G,
        .radius               = K26A_R_EARTH_EQU,
        .polar_radius         = K26A_R_EARTH_POL,
        .j2                   = 1.08262668e-3,
        .rotation_model_name  = "iau2018:earth"
    },

    /* Moon.
     * GM: DE441. Radius: mean 1 737.4 km.
     * J2: 2.0321e-4 (GRAIL mission, Konopliv et al. 2014). */
    {
        .name                 = "moon",
        .kind                 = K26ASTRO_BODY_MOON,
        .naif_id              = 301,
        .gm                   = 4.9028000661e12,
        .mass                 = 4.9028000661e12 / K26A_G,
        .radius               = 1.7374e6,
        .polar_radius         = 1.7374e6,
        .j2                   = 2.0321e-4,
        .rotation_model_name  = "iau2018:moon"
    },

    /* Mars.
     * GM: DE441. Radius: 3 396.2 km equatorial. J2: 1.96e-3
     * (Konopliv et al. 2016, MRO + InSight). */
    {
        .name                 = "mars",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 499,
        .gm                   = 4.2828375816e13,
        .mass                 = 4.2828375816e13 / K26A_G,
        .radius               = 3.3962e6,
        .polar_radius         = 3.3762e6,
        .j2                   = 1.96045e-3,
        .rotation_model_name  = "iau2018:mars"
    },

    /* Jupiter.
     * GM: IAU 2015 R.B3 nominal. Radius: 71 492 km equatorial.
     * J2: 1.4736e-2 (Juno + pre-mission). */
    {
        .name                 = "jupiter",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 599,
        .gm                   = K26A_GM_JUPITER,
        .mass                 = K26A_GM_JUPITER / K26A_G,
        .radius               = 7.1492e7,
        .polar_radius         = 6.6854e7,
        .j2                   = 1.4736e-2,
        .rotation_model_name  = "iau2018:jupiter"
    },

    /* Saturn.
     * GM: DE441 (~3.793e16). Radius: 60 268 km equatorial.
     * J2: 1.6298e-2 (Cassini). */
    {
        .name                 = "saturn",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 699,
        .gm                   = 3.7931187000e16,
        .mass                 = 3.7931187000e16 / K26A_G,
        .radius               = 6.0268e7,
        .polar_radius         = 5.4364e7,
        .j2                   = 1.6298e-2,
        .rotation_model_name  = "iau2018:saturn"
    },

    /* Uranus.
     * GM: DE441 (5.7940e15). Radius: 25 559 km equatorial.
     * J2: 3.3434e-3 (Voyager 2 + Earth-based). */
    {
        .name                 = "uranus",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 799,
        .gm                   = 5.7939393880e15,
        .mass                 = 5.7939393880e15 / K26A_G,
        .radius               = 2.5559e7,
        .polar_radius         = 2.4973e7,
        .j2                   = 3.3434e-3,
        .rotation_model_name  = "iau2018:uranus"
    },

    /* Neptune.
     * GM: DE441 (6.8351e15). Radius: 24 764 km. J2: 3.411e-3. */
    {
        .name                 = "neptune",
        .kind                 = K26ASTRO_BODY_PLANET,
        .naif_id              = 899,
        .gm                   = 6.8365271007e15,
        .mass                 = 6.8365271007e15 / K26A_G,
        .radius               = 2.4764e7,
        .polar_radius         = 2.4341e7,
        .j2                   = 3.411e-3,
        .rotation_model_name  = "iau2018:neptune"
    },

    /* Pluto.
     * GM: DE441. Radius: 1 188.3 km (New Horizons). J2: tiny. */
    {
        .name                 = "pluto",
        .kind                 = K26ASTRO_BODY_PLANET,   /* dwarf, but kind is informational */
        .naif_id              = 999,
        .gm                   = 8.696e11,
        .mass                 = 8.696e11 / K26A_G,
        .radius               = 1.1883e6,
        .polar_radius         = 1.1883e6,
        .j2                   = 0.0,
        .rotation_model_name  = "iau2018:pluto"
    },

    /* NULL terminator — clients iterate until name == NULL. */
    { .name = NULL }
};

const K26AstroMajorBody *k26astro_major_bodies(void)
{
    return MAJOR_BODIES;
}

int k26astro_major_body_count(void)
{
    int n = 0;
    for (const K26AstroMajorBody *m = MAJOR_BODIES; m->name; m++) n++;
    return n;
}

const K26AstroMajorBody *k26astro_major_body_find(const char *name)
{
    if (!name) return NULL;
    for (const K26AstroMajorBody *m = MAJOR_BODIES; m->name; m++) {
        if (strcasecmp(m->name, name) == 0) return m;
    }
    return NULL;
}

int k26astro_body_load_major(K26AstroBody *b, const char *name)
{
    if (!b || !name) return 1;
    const K26AstroMajorBody *m = k26astro_major_body_find(name);
    if (!m) return 2;
    k26astro_body_init(b);
    /* strncpy the lower-case name; clients writing different
     * presentation casing can override afterwards. */
    size_t nlen = strlen(m->name);
    if (nlen >= K26ASTRO_BODY_NAME_MAX) nlen = K26ASTRO_BODY_NAME_MAX - 1;
    memcpy(b->name, m->name, nlen);
    b->name[nlen] = '\0';
    b->kind          = m->kind;
    b->gm            = m->gm;
    b->mass          = m->mass;
    b->radius        = m->radius;
    b->polar_radius  = m->polar_radius;
    b->j2            = m->j2;
    b->ephem_naif_id = m->naif_id;
    b->attitude_mode = K26ASTRO_ATT_ROTATION_MODEL;
    /* rotation_model_id is resolved by libk26astro_rt at world
     * construction via rotation_model.h:k26astro_rotation_lookup. */
    b->rotation_model_id = 0;
    return 0;
}
