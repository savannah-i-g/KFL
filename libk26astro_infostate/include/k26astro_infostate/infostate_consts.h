/* k26astro_infostate/infostate_consts.h — compile-time tunables.
 *
 * Speed-of-light value is delegated to libk26astro_core's K26A_C so
 * a single canonical constant is used across the astro tier. The
 * fixed-point tolerance and iteration cap are chosen to converge
 * cleanly for engagement-bubble light-times (microseconds to
 * minutes); both can be tightened by consumers in REFERENCED-mode
 * builds where bit-identical replay is required.
 *
 * K26ASTRO_INFOSTATE_MAX_TARGETS bounds the per-infostate target
 * table. The default 64 comfortably covers a two-faction scenario
 * (10 craft per faction plus a handful of celestial bodies).
 * Larger scenarios can rebuild with a higher cap. */
#ifndef K26ASTRO_INFOSTATE_CONSTS_H
#define K26ASTRO_INFOSTATE_CONSTS_H

#define K26ASTRO_INFOSTATE_MAX_TARGETS              64

#define K26ASTRO_INFOSTATE_DEFAULT_HISTORY_CAPACITY 1024
#define K26ASTRO_INFOSTATE_MIN_HISTORY_CAPACITY     2

#define K26ASTRO_INFOSTATE_FIXED_POINT_TOL_S        1.0e-12
#define K26ASTRO_INFOSTATE_FIXED_POINT_MAX_ITERS    16

#endif /* K26ASTRO_INFOSTATE_CONSTS_H */
