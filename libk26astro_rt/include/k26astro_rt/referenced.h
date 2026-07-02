/* libk26astro_rt — REFERENCED determinism mode op-log.
 *
 * The third determinism mode (alongside FAST and PORTABLE). When a
 * world is created with K26ASTRO_MODE_REFERENCED and a log path is
 * set via k26astro_world_set_ref_log_path(), every world-stepping
 * event writes a packed-binary record to the log. Two runs from the
 * same initial conditions produce byte-identical log files — the
 * cross-libc / cross-platform / cross-build determinism regression
 * net.
 *
 * Why packed binary instead of textual `%a` records:
 *   - musl strtod / glibc strtod round-trip differently in the last
 *     ULP on some bit patterns; text-based comparison adds a
 *     spurious divergence source.
 *   - Hex-encoded IEEE-754 doubles ARE round-trip-stable, but writing
 *     them as ASCII costs ~3x storage + downstream readers need a
 *     custom parser. Binary is the simpler ground truth.
 *
 * Log structure:
 *   - K26ASTRO_REF_MAGIC ("K26REF\0\0") + uint32 version + uint32
 *     reserved at the file head.
 *   - One IC block (K26AstroRefIC) — body count + per-body
 *     {mass, gm, pos{x,y,z}, vel{x,y,z}, name[K26ASTRO_BODY_NAME_MAX]}.
 *   - A sequence of K26AstroOpRecord packed-binary records.
 *
 * Replay tool: tools/k26astro_replay/.  `--check` compares two logs
 * byte-for-byte; `--dump --sample N` produces a human-readable
 * sampled dump for triage.
 */
#ifndef K26ASTRO_RT_REFERENCED_H
#define K26ASTRO_RT_REFERENCED_H

#include "k26astro_rt/world.h"

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 8-byte magic at offset 0. Use a tail null so casual `cat` shows a
 * recognisable prefix. */
#define K26ASTRO_REF_MAGIC      "K26REF\0\0"
#define K26ASTRO_REF_MAGIC_LEN  8
#define K26ASTRO_REF_VERSION    1u

/* Op identifiers. Stable across builds — appending only, never
 * renumbering. Used as packed `uint16_t op_id` in record headers. */
typedef enum {
    K26ASTRO_OP_IC_BLOCK     = 0,   /* initial-conditions header */
    K26ASTRO_OP_STEP_BEGIN   = 1,   /* dt requested */
    K26ASTRO_OP_STEP_END     = 2,   /* post-step world.time */
    K26ASTRO_OP_ENCOUNTER    = 3    /* close-encounter pair detect */
} K26AstroRefOpId;

/* Per-op argument arity (number of double args following the
 * header). Static-const at compile time — `op_arity[op_id]` gives
 * the args count. Replay tools dispatch via this table. */
extern const uint8_t k26astro_ref_op_arity[];
/* Number of valid op ids (size of the arity table). */
extern const uint8_t k26astro_ref_op_count;

/* Packed-binary op record header. The on-disk layout is exactly:
 *   op_id (u16) | body_count (u16) | step_idx (u32)
 *   body_count × int32 body_indices
 *   k26astro_ref_op_arity[op_id] × uint64 (hex-stable IEEE-754
 *     double bit patterns; writer uses memcpy of a union, no float
 *     arithmetic in serialization). */
typedef struct {
    uint16_t op_id;
    uint16_t body_count;
    uint32_t step_idx;
} K26AstroRefOpHeader;

/* ---- World-side API --------------------------------------------- */

/* Set the log file path. Caller-owned string. The world opens the
 * file for write+append on the next step and keeps it open until
 * destroy. Set the path to NULL to disable logging. Only valid
 * when the world's mode is K26ASTRO_MODE_REFERENCED — returns
 * K26ASTRO_RT_E_BAD_ARG otherwise.
 *
 * The IC block is written on the FIRST call to k26astro_world_step
 * after the path is set, so callers must add all bodies before
 * stepping. */
int k26astro_world_set_ref_log_path(K26AstroWorld *world, const char *path);

/* Force-flush the ring buffer to disk + close the fd. Normally
 * called from world_destroy; exposed for test/debug use. After
 * close, the path is detached — set_ref_log_path() can be called
 * again to start a fresh log. */
int k26astro_world_ref_log_close(K26AstroWorld *world);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_REFERENCED_H */
