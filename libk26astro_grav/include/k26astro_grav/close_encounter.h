/* k26astro_grav/close_encounter.h — Hill-radius proximity detector.
 *
 * The detector primitive consumed by MERCURIUS handoff orchestration.
 * Scans the body array and returns the index of any body that's
 * currently inside another body's effective Hill sphere (scaled by
 * hill_factor; typical value 3.0 for "close enough that MERCURIUS
 * should rewind WH").
 *
 * Returns -1 if no close encounter is happening. If multiple bodies
 * are simultaneously in close encounter, returns the most-deeply-
 * penetrating one (smallest r/r_Hill ratio).
 *
 * The handoff orchestration itself (rewind WH substep, switch to
 * IAS15, propagate through the close encounter, switch back) lives
 * in libk26astro_rt. */
#ifndef K26ASTRO_GRAV_CLOSE_ENCOUNTER_H
#define K26ASTRO_GRAV_CLOSE_ENCOUNTER_H

#include "k26astro_grav/forces.h"

#ifdef __cplusplus
extern "C" {
#endif

int k26astro_grav_close_encounter(const K26AstroGravView *view,
                                   double hill_factor);

#ifdef __cplusplus
}
#endif

#endif
