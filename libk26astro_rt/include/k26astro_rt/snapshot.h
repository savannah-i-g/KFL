/* libk26astro_rt — versioned binary snapshot format.
 *
 * On-disk layout (little-endian; reader byte-swaps if endian probe
 * mismatches). See the runtime reference manual §snapshot format.
 *
 *   struct K26AstroSnapshotHeader {
 *       char     magic[8];      "K26ASNP\0"
 *       uint32_t version;       K26ASTRO_SNAPSHOT_VERSION_CURRENT
 *       uint32_t endian_probe;  0x01020304
 *       uint32_t flags;         bit 0 = Q64_64; bit 1 = FAST mode
 *       uint32_t body_count;
 *       K26AstroEpoch epoch;    24 bytes
 *       char     integrator[16];
 *       uint8_t  reserved[8];   header total = 80 bytes
 *   };
 *
 * Each body record:
 *   uint16_t type_tag;    0 = base K26AstroBody, future = extensions
 *   uint16_t flags;       reserved
 *   <K26AstroBody fields serialised one at a time>
 *
 * Version 1 is the only version v0.1 reads/writes. A v2 reader would
 * dispatch by `version` to a migration routine. */
#ifndef K26ASTRO_RT_SNAPSHOT_H
#define K26ASTRO_RT_SNAPSHOT_H

#include "k26astro_rt/world.h"

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_SNAPSHOT_MAGIC          "K26ASNP\0"
#define K26ASTRO_SNAPSHOT_VERSION_CURRENT 1u
#define K26ASTRO_SNAPSHOT_ENDIAN_PROBE   0x01020304u
#define K26ASTRO_SNAPSHOT_HEADER_BYTES   80

/* Snapshot flags. */
#define K26ASTRO_SNAPSHOT_FLAG_Q64_64    (1u << 0)
#define K26ASTRO_SNAPSHOT_FLAG_FAST_MODE (1u << 1)

/* Write the world's state to `path`. Returns 0 on success or a
 * negative K26ASTRO_RT_E_* code. In PORTABLE mode the output is
 * byte-identical for the same world state — this is Test 9's gate. */
int k26astro_world_snapshot_save(const K26AstroWorld *world,
                                  const char *path);

/* Load a world from `path`. The returned world is a fresh allocation
 * (the caller owns it; destroy with k26astro_world_destroy). On
 * version mismatch returns NULL — v2 worlds need their own loader. */
K26AstroWorld *k26astro_world_snapshot_load(const char *path);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_SNAPSHOT_H */
