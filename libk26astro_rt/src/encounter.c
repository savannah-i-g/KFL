/* encounter.c — MERCURIUS close-encounter primitives.
 *
 * Implements the Rein-Tamayo 2019 quintic smoothstep K(y) over the
 * transition window [y_inner, y_outer] in Hill-radius units, plus
 * the per-pair Hill radius and the per-pair encounter session
 * tracker. The orchestration that uses these primitives to actually
 * split the force between WH and IAS15 lives in orbit_step.c.
 *
 * Reference: Rein & Tamayo (2019), MNRAS 489:4632-4640,
 * "MERCURIUS: a hybrid integrator for long-term planetary
 * simulations including close encounters." */
#include "encounter_internal.h"

#include "k26astro_core/pos.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double k26astro_mercurius_K(double y, double y_inner, double y_outer)
{
    if (!(y_outer > y_inner)) return 1.0;
    if (y <= y_inner) return 1.0;
    if (y >= y_outer) return 0.0;
    double x = (y - y_inner) / (y_outer - y_inner);
    /* Quintic smoothstep S(x) = 10 x^3 - 15 x^4 + 6 x^5.
     * K(y) = 1 - S(x) so K = 1 at y_inner (full IAS15) and
     * K = 0 at y_outer (full WH). C^2 at both endpoints. */
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;
    double s  = 10.0 * x3 - 15.0 * x4 + 6.0 * x5;
    return 1.0 - s;
}

double k26astro_mercurius_hill_radius(const K26AstroBody *i,
                                       const K26AstroBody *j,
                                       double m_central)
{
    if (!i || !j) return 0.0;
    /* Semi-major axis approximated by current separation (good
     * within an order of magnitude near the transition region — the
     * exact value would require fitting an osculating conic, which
     * is too expensive per-pair per-step). Document this as a
     * heuristic in the encounter primitive. */
    K26V3 r = k26astro_pos_sub(&i->pos, &j->pos);
    double a_ij = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    if (!(a_ij > 0.0))      return 0.0;
    if (!(m_central > 0.0)) return 0.0;
    double m_sum = i->mass + j->mass;
    if (!(m_sum > 0.0))     return 0.0;
    return a_ij * cbrt(m_sum / (3.0 * m_central));
}

/* Grow the world's encounter buffer in-place. */
static int ensure_encounter_capacity_(K26AstroWorld *world, int need)
{
    if (need <= world->cap_encounters) return 0;
    int new_cap = world->cap_encounters ? world->cap_encounters * 2 : 8;
    while (new_cap < need) new_cap *= 2;
    K26AstroEncounter *p = (K26AstroEncounter *)realloc(
        world->encounters, (size_t)new_cap * sizeof(K26AstroEncounter));
    if (!p) return -1;
    world->encounters     = p;
    world->cap_encounters = new_cap;
    return 0;
}

int k26astro_mercurius_detect(K26AstroWorld *world)
{
    if (!world) return 0;
    int n = world->grav.n_bodies;
    if (n < 2) { world->n_encounters = 0; return 0; }

    /* Central body = the largest GM. In the solar system this is the
     * Sun; for moon-around-planet sub-systems the caller should set
     * up a hierarchical world (out of scope in v0.1). */
    double m_central = 0.0;
    for (int k = 0; k < n; k++) {
        double m = world->grav.bodies[k].mass;
        if (m > m_central) m_central = m;
    }
    if (!(m_central > 0.0)) { world->n_encounters = 0; return 0; }

    world->n_encounters = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            const K26AstroBody *bi = &world->grav.bodies[i];
            const K26AstroBody *bj = &world->grav.bodies[j];
            double rh = k26astro_mercurius_hill_radius(bi, bj, m_central);
            if (!(rh > 0.0)) continue;
            K26V3 r = k26astro_pos_sub(&bi->pos, &bj->pos);
            double d = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
            double y = d / rh;
            if (y >= world->mercurius_outer_factor) continue;
            /* In transition or fully near. Record with K weight. */
            double K = k26astro_mercurius_K(y,
                                              world->mercurius_hill_factor,
                                              world->mercurius_outer_factor);
            int slot = world->n_encounters;
            if (ensure_encounter_capacity_(world, slot + 1) != 0) {
                return slot;
            }
            world->encounters[slot] = (K26AstroEncounter){
                .i = i, .j = j,
                .y_last = y,
                .k_weight = K,
                .active = 1
            };
            world->n_encounters = slot + 1;
        }
    }
    return world->n_encounters;
}
