/* close_encounter.c — Hill-radius proximity detector.
 *
 * Scans the body array and returns the index of any body whose
 * distance from its parent body's center is within hill_factor × R_Hill.
 * If multiple bodies are simultaneously in close encounter, returns
 * the most-deeply-penetrating one (smallest r / r_Hill ratio).
 *
 * Hill radius approximation: r_Hill ≈ a · ∛(m_secondary / (3 · m_primary)).
 * Here we use the body's CURRENT distance from the parent (not the
 * orbital semi-major axis) — this is the "Hill sphere at this instant"
 * approximation that the MERCURIUS handoff detector wants.
 *
 * Returns -1 if no close encounter, otherwise the body index. */
#include "k26astro_grav/close_encounter.h"
#include "k26astro_body/body.h"
#include "k26astro_core/pos.h"

#include <math.h>

int k26astro_grav_close_encounter(const K26AstroGravView *view,
                                   double hill_factor)
{
    if (!view || view->n < 2 || hill_factor <= 0.0) return -1;

    int best_idx = -1;
    double best_ratio = hill_factor;   /* must be < this to register */

    for (int i = 0; i < view->n; i++) {
        const K26AstroBody *bi = &view->bodies[i];
        int p = bi->parent_body_idx;
        if (p < 0 || p >= view->n || p == i) continue;
        const K26AstroBody *parent = &view->bodies[p];
        if (parent->gm <= 0.0 || bi->gm <= 0.0) continue;

        K26V3 r = k26astro_pos_sub(&bi->pos, &parent->pos);
        double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        if (r_mag <= 0.0) continue;

        /* Hill radius using current separation as the orbit scale. */
        double r_hill = r_mag * cbrt(bi->gm / (3.0 * parent->gm));
        if (r_hill <= 0.0) continue;

        double ratio = r_mag / r_hill;
        if (ratio < best_ratio) {
            best_ratio = ratio;
            best_idx   = i;
        }
    }
    return best_idx;
}
