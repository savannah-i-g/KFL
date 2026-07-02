/* k26astro_body/body.h — K26AstroBody taxonomy + state.
 *
 * A body is a discrete simulation entity with translational state
 * (position + velocity) and optionally attitude state (quaternion +
 * angular velocity). The kind (planet, moon, spacecraft, …) drives
 * default behaviour (atmosphere flag, attitude mode, snapshot
 * serialisation) but doesn't affect the integrator — gravity sums
 * over all bodies regardless of kind.
 *
 * Bodies live in a contiguous array owned by `K26AstroWorld` (lib
 * libk26astro_rt). The body struct here exposes the data layout; the
 * world holds the lifetime. Identity ids are world-local and stable
 * across array compactions.
 *
 * The major-body lookup table at the bottom of this header gives
 * authoritative physical constants (mass, equatorial radius, J2) +
 * the NAIF id used for ephemeris queries, sourced from:
 *   - IAU 2015 Resolution B3 nominal values (GM constants)
 *   - JPL Solar System Dynamics planetary data sheet (radius, J2)
 *   - NAIF id assignments (canonical)
 */
#ifndef K26ASTRO_BODY_BODY_H
#define K26ASTRO_BODY_BODY_H

#include <stddef.h>
#include <stdint.h>

#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Taxonomy --------------------------------------------------- */

typedef enum {
    K26ASTRO_BODY_UNKNOWN      = 0,
    K26ASTRO_BODY_STAR         = 1,
    K26ASTRO_BODY_PLANET       = 2,
    K26ASTRO_BODY_MOON         = 3,
    K26ASTRO_BODY_ASTEROID     = 4,
    K26ASTRO_BODY_COMET        = 5,
    K26ASTRO_BODY_SPACECRAFT   = 6,
    K26ASTRO_BODY_DEBRIS       = 7,
    K26ASTRO_BODY_BARYCENTRE   = 8     /* abstract / virtual */
} K26AstroBodyKind;

/* ---- Attitude modes -------------------------------------------- *
 *
 * Drives how `K26AstroBody.attitude` evolves over time:
 *   FREE             — integrate ω forward each substep (used by
 *                      libk26astro_grav); ω updated by torques.
 *   ROTATION_MODEL   — attitude derived from IAU rotation model at
 *                      epoch (planets, moons with known spin).
 *   FIXED_IN_WORLD   — attitude held constant (dev / debug).
 *   FIXED_IN_PARENT  — attitude tracks the parent body's frame
 *                      (tidally-locked moons, lagrange-point objects). */
typedef enum {
    K26ASTRO_ATT_UNKNOWN          = 0,
    K26ASTRO_ATT_FREE             = 1,
    K26ASTRO_ATT_ROTATION_MODEL   = 2,
    K26ASTRO_ATT_FIXED_IN_WORLD   = 3,
    K26ASTRO_ATT_FIXED_IN_PARENT  = 4
} K26AstroAttitudeMode;

/* ---- Body struct ----------------------------------------------- */

#define K26ASTRO_BODY_NAME_MAX 32

typedef struct K26AstroBody {
    /* ---- Identity --------------------------------------------- */
    uint64_t          id;            /* world-local unique id */
    char              name[K26ASTRO_BODY_NAME_MAX];
    K26AstroBodyKind  kind;

    /* ---- Translational state ---------------------------------- */
    K26AstroPos       pos;
    K26V3             vel;             /* m/s in the world's reference frame */

    /* ---- Mass + shape ----------------------------------------- */
    double            mass;            /* kg */
    double            gm;              /* GM in m^3/s^2; redundant w/ mass but
                                          stored to avoid the G uncertainty
                                          (G has CODATA-level uncertainty; GM
                                          products are IAU-fixed exactly) */
    double            radius;          /* equatorial, m */
    double            polar_radius;    /* polar, m; 0 → assume spherical */
    double            j2;              /* second zonal harmonic (dimensionless) */
    /* future: full spherical-harmonics gravity model gated by
     * libk26astro_grav extension */

    /* ---- Attitude state (optional) --------------------------- */
    K26Quat           attitude;        /* identity if not tracked */
    K26V3             omega;           /* rad/s in body frame; zero if not tracked */
    K26AstroAttitudeMode attitude_mode;

    /* For ROTATION_MODEL: id into the IAU rotation table.
     * 0 → not assigned. See rotation_model.h. */
    int               rotation_model_id;

    /* ---- Atmosphere flag (consumed by future libk26astro_atmos) */
    uint8_t           has_atmos;

    /* ---- Ephemeris override ---------------------------------- */
    /* When > 0, the runtime treats the body as on-rails and pulls
     * its position from libk26astro_ephem using this NAIF id rather
     * than integrating it. Standard for planets in a spacecraft
     * simulation — Earth/Mars positions come from DE441; the
     * spacecraft body is integrated. */
    int               ephem_naif_id;

    /* ---- On-rails / patched-conic state (libk26astro_conics) -- */
    uint8_t           on_rails;
    int               parent_body_idx;  /* SOI parent; -1 if none */
} K26AstroBody;

/* ---- Construction + identity ----------------------------------- */

/* Zero-initialise a body with sensible defaults: identity attitude,
 * zero omega, kind = UNKNOWN. The caller fills the rest. */
void k26astro_body_init(K26AstroBody *b);

/* Update `b`'s mass and recompute `b->gm` as `K26A_G * mass_kg`.
 * Negative input clamps to 0. No-op on NULL.
 *
 * Use for bodies whose mass is being mutated at runtime, e.g.
 * spacecraft losing propellant during a burn. For celestial bodies
 * loaded from the major-body table, the table's IAU-fixed GM is
 * fundamentally more precise than `G * mass` (G carries CODATA
 * uncertainty); do not call this on celestial bodies unless their
 * mass is genuinely being mutated. */
void k26astro_body_set_mass(K26AstroBody *b, double mass_kg);

/* Set `b` from a major-body table entry (see below). Returns 0 on
 * success, non-zero if the name isn't in the table. Position + vel
 * are NOT set — the caller pulls those from ephem or computes from
 * elements. */
int  k26astro_body_load_major(K26AstroBody *b, const char *name);

/* Lookup a body by name in the body's contiguous storage. Returns
 * the index, or -1 if not found. */
int  k26astro_body_find_by_name(const K26AstroBody *bodies, int n,
                                 const char *name);

/* ---- Major-body table ----------------------------------------- *
 *
 * Authoritative physical constants for the Sun + 8 planets + Pluto +
 * Moon. Used as the default seed for new worlds — KFL programs that
 * write `astro_body earth { kind = planet; … }` without overriding
 * mass/radius pick up these values. Future updates (IAU revisions to
 * nominal constants) come through this table; downstream code never
 * embeds raw mass / radius literals. */
typedef struct {
    const char       *name;       /* lowercase; "earth", "moon", "sun", … */
    K26AstroBodyKind  kind;
    int               naif_id;    /* canonical NAIF id (see SPK reader) */
    double            gm;         /* m^3/s^2 — IAU 2015 nominal */
    double            mass;       /* kg — derived gm / G */
    double            radius;     /* m — equatorial */
    double            polar_radius;
    double            j2;
    /* IAU 2018 rotation model name — looked up at world construction
     * via rotation_model.h:k26astro_rotation_lookup. NULL means no
     * registered rotation model (used for the Sun's surface-frame
     * special case). */
    const char       *rotation_model_name;
} K26AstroMajorBody;

/* Read-only access to the table; NULL terminator. */
const K26AstroMajorBody *k26astro_major_bodies(void);

/* Total count (excluding the NULL terminator). */
int  k26astro_major_body_count(void);

/* Find a major-body entry by name (case-insensitive). NULL if not
 * found. */
const K26AstroMajorBody *k26astro_major_body_find(const char *name);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_BODY_BODY_H */
