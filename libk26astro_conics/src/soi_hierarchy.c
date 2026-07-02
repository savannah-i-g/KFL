/* soi_hierarchy.c — parent_body_idx walker.
 *
 * For v0.1 the implementation is intentionally minimal: the
 * parent_body_idx field on K26AstroBody is the source of truth.
 * Future versions may consult the body's position to dynamically
 * re-bind (e.g. a Mars-bound probe that's outside Mars' SOI is
 * not really Mars's child) — but in v0.1, the patched-conic
 * machinery in libk26astro_rt will manage parent re-binding when
 * SOI crossings happen, and this function just reflects the
 * static state. */
#include "k26astro_conics/soi.h"

#include <stddef.h>

int k26astro_soi_parent(const K26AstroBody *bodies, int n_bodies,
                        int child_idx, double t)
{
    (void)t;   /* reserved for time-varying SOI hierarchies in v0.2+ */
    if (!bodies) return -1;
    if (child_idx < 0 || child_idx >= n_bodies) return -1;
    int p = bodies[child_idx].parent_body_idx;
    /* Validate: parent index must be in range and not self-cycling. */
    if (p < 0 || p >= n_bodies) return -1;
    if (p == child_idx) return -1;
    return p;
}
