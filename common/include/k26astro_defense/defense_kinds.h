/* defense_kinds.h — payload kind-tag registry for the KFL_Defense
 * libraries.
 *
 * Every defense subsystem binds into libk26astro_vehicle's generic
 * PAYLOAD composition slots (k26astro_vehicle/payload.h) and
 * identifies itself through the open kind-tag namespace carried in
 * the K26AstroPayload base. This header is the single registry for
 * the defense stack's tags so co-linked defense libraries never
 * collide.
 *
 * Tag scheme: high half 0x4446 ("DF"); low half enumerates. 0 is
 * reserved by payload.h as "none". */
#ifndef K26ASTRO_DEFENSE_KINDS_H
#define K26ASTRO_DEFENSE_KINDS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
    K26ASTRO_DEFENSE_KIND_IMPACTOR      = 0x44460001, /* kinetic impactor */
    K26ASTRO_DEFENSE_KIND_LASER         = 0x44460002,
    /* 0x44460003 is reserved and must not be reassigned. */
    K26ASTRO_DEFENSE_KIND_DECOY         = 0x44460004,
    K26ASTRO_DEFENSE_KIND_JAMMER        = 0x44460005,
    K26ASTRO_DEFENSE_KIND_DAZZLER       = 0x44460006,
    K26ASTRO_DEFENSE_KIND_DETECT_SENSOR = 0x44460007,
    K26ASTRO_DEFENSE_KIND_INFOSTATE     = 0x44460008
};

/* Effector-class tags occupy [IMPACTOR, DAZZLER]; sensing/state tags
 * follow. A dispatch loop that needs "is this payload one of ours
 * and an effector" checks the range below. */
#define K26ASTRO_DEFENSE_KIND_IS_EFFECTOR(k) \
    ((k) >= K26ASTRO_DEFENSE_KIND_IMPACTOR && (k) <= K26ASTRO_DEFENSE_KIND_DAZZLER)

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DEFENSE_KINDS_H */
