/* referenced.c — REFERENCED determinism mode op-log writer.
 *
 * Ring-buffered append. Each writer call (k26astro_ref_emit_*)
 * marshals the op record into the buffer; the ring flushes to the
 * fd at world-step boundaries OR when the next write wouldn't fit.
 *
 * Determinism contract: all doubles in the log are written via a
 * union (no float arithmetic in serialization) so the bit pattern
 * reaches disk exactly as it lives in memory. Two runs from the
 * same initial conditions produce byte-identical log files.
 *
 * Layout note: the file is opened with O_TRUNC | O_WRONLY — the
 * v0.1 surface is "write a complete log per simulation"; mid-run
 * appending across processes isn't supported. */
#define _GNU_SOURCE
#include "k26astro_rt/referenced.h"
#include "world_internal.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Op arity table — number of double args following the header for
 * each op id. */
const uint8_t k26astro_ref_op_arity[] = {
    [K26ASTRO_OP_IC_BLOCK]   = 0,    /* IC block carries its own variable payload */
    [K26ASTRO_OP_STEP_BEGIN] = 1,    /* dt */
    [K26ASTRO_OP_STEP_END]   = 1,    /* world.time after step (TDB seconds) */
    [K26ASTRO_OP_ENCOUNTER]  = 1     /* k_weight */
};
const uint8_t k26astro_ref_op_count =
    (uint8_t)(sizeof(k26astro_ref_op_arity) / sizeof(k26astro_ref_op_arity[0]));

#define K26ASTRO_REF_RING_SIZE  (16 * 1024)

/* Internal ring buffer + fd. Allocated lazily on first emit. */
typedef struct {
    FILE       *fp;            /* stream over the log file */
    uint8_t    *buf;           /* ring storage; 16 KB */
    size_t      used;          /* bytes pending flush */
    uint32_t    step_idx;      /* monotonic counter, ++ on STEP_BEGIN */
    int         header_written;
    int         ic_written;
    char       *path;          /* owned copy of caller's path */
} K26AstroRefCtx;

/* Helper: append `len` raw bytes to the ring; flush if it would
 * overflow. */
static int ring_append_(K26AstroRefCtx *ctx, const void *data, size_t len)
{
    if (!ctx || !ctx->fp || !ctx->buf) return -1;
    if (ctx->used + len > K26ASTRO_REF_RING_SIZE) {
        /* Flush the existing tail then continue. */
        if (ctx->used > 0) {
            if (fwrite(ctx->buf, 1, ctx->used, ctx->fp) != ctx->used)
                return -1;
            ctx->used = 0;
        }
        if (len > K26ASTRO_REF_RING_SIZE) {
            /* Oversize write: stream directly. */
            if (fwrite(data, 1, len, ctx->fp) != len) return -1;
            return 0;
        }
    }
    memcpy(ctx->buf + ctx->used, data, len);
    ctx->used += len;
    return 0;
}

/* Helper: pack a double into its IEEE-754 bit pattern as uint64 and
 * append. No float arithmetic — bit-stable across CPUs. */
static int ring_append_double_(K26AstroRefCtx *ctx, double d)
{
    union { double d; uint64_t u; } cvt;
    cvt.d = d;
    return ring_append_(ctx, &cvt.u, sizeof(uint64_t));
}

/* Helper: pack an integer (host-endian, fixed-width) and append. */
static int ring_append_u32_(K26AstroRefCtx *ctx, uint32_t v)
{ return ring_append_(ctx, &v, sizeof v); }
static int ring_append_i32_(K26AstroRefCtx *ctx, int32_t v)
{ return ring_append_(ctx, &v, sizeof v); }

/* Allocate the ctx, open the file, write the file header.
 * Returns 0 on success, -1 on error. */
static int ref_ctx_open_(K26AstroWorld *world)
{
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || !ctx->path) return -1;
    ctx->fp = fopen(ctx->path, "wb");
    if (!ctx->fp) {
        fprintf(stderr, "k26astro_ref: cannot open %s: %s\n",
                ctx->path, strerror(errno));
        return -1;
    }
    if (!ctx->buf) {
        ctx->buf = (uint8_t *)malloc(K26ASTRO_REF_RING_SIZE);
        if (!ctx->buf) { fclose(ctx->fp); ctx->fp = NULL; return -1; }
    }
    ctx->used = 0;
    /* File header. */
    if (ring_append_(ctx, K26ASTRO_REF_MAGIC, K26ASTRO_REF_MAGIC_LEN) != 0)
        return -1;
    uint32_t version  = K26ASTRO_REF_VERSION;
    uint32_t reserved = 0;
    if (ring_append_u32_(ctx, version)  != 0) return -1;
    if (ring_append_u32_(ctx, reserved) != 0) return -1;
    ctx->header_written = 1;
    return 0;
}

/* Write the IC block: body_count + per-body packed record. Called
 * lazily on the first step after a log path is set. Writes through
 * the ring like everything else. */
static int ref_write_ic_(K26AstroWorld *world)
{
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || ctx->ic_written) return 0;

    /* Op header: op_id=IC, body_count, step_idx=0. body_count is
     * really the IC's payload, not a body_indices count; here we
     * encode the body count directly in body_count. */
    int n = world->grav.n_bodies;
    if (n < 0 || n > 65535) return -1;
    K26AstroRefOpHeader hdr = {
        .op_id     = (uint16_t)K26ASTRO_OP_IC_BLOCK,
        .body_count = (uint16_t)n,
        .step_idx  = 0
    };
    if (ring_append_(ctx, &hdr, sizeof hdr) != 0) return -1;

    /* Per-body payload: name[K26ASTRO_BODY_NAME_MAX] + mass + gm +
     * pos.{x,y,z} as doubles (m) + vel.{x,y,z}. */
    for (int i = 0; i < n; i++) {
        K26AstroBody *b = &world->grav.bodies[i];
        if (ring_append_(ctx, b->name, K26ASTRO_BODY_NAME_MAX) != 0) return -1;
        if (ring_append_double_(ctx, b->mass) != 0) return -1;
        if (ring_append_double_(ctx, b->gm)   != 0) return -1;
        K26V3 r = k26astro_pos_to_m_approx(&b->pos);
        if (ring_append_double_(ctx, r.x) != 0) return -1;
        if (ring_append_double_(ctx, r.y) != 0) return -1;
        if (ring_append_double_(ctx, r.z) != 0) return -1;
        if (ring_append_double_(ctx, b->vel.x) != 0) return -1;
        if (ring_append_double_(ctx, b->vel.y) != 0) return -1;
        if (ring_append_double_(ctx, b->vel.z) != 0) return -1;
    }
    ctx->ic_written = 1;
    return 0;
}

/* Public API ----------------------------------------------------- */

int k26astro_world_set_ref_log_path(K26AstroWorld *world, const char *path)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    if (world->mode != K26ASTRO_MODE_REFERENCED) return -K26ASTRO_RT_E_BAD_ARG;
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx) {
        ctx = (K26AstroRefCtx *)calloc(1, sizeof(*ctx));
        if (!ctx) return -K26ASTRO_RT_E_OOM;
        world->ref_ctx = ctx;
    }
    /* Close any existing log on path change. */
    if (ctx->fp) { fflush(ctx->fp); fclose(ctx->fp); ctx->fp = NULL; }
    free(ctx->path);
    ctx->path = path ? strdup(path) : NULL;
    ctx->header_written = 0;
    ctx->ic_written     = 0;
    ctx->step_idx       = 0;
    ctx->used           = 0;
    return K26ASTRO_RT_OK;
}

int k26astro_world_ref_log_close(K26AstroWorld *world)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx) return K26ASTRO_RT_OK;
    if (ctx->fp) {
        if (ctx->used > 0) {
            (void)fwrite(ctx->buf, 1, ctx->used, ctx->fp);
            ctx->used = 0;
        }
        fflush(ctx->fp);
        fclose(ctx->fp);
        ctx->fp = NULL;
    }
    free(ctx->buf); ctx->buf = NULL;
    free(ctx->path); ctx->path = NULL;
    free(ctx);
    world->ref_ctx = NULL;
    return K26ASTRO_RT_OK;
}

/* ---- Internal emit helpers (called from world.c, orbit_step.c) ---- */

/* Lazy-init: open the file + write header + IC on first emit. */
static int ref_lazy_init_(K26AstroWorld *world)
{
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || !ctx->path) return -1;
    if (!ctx->fp) {
        if (ref_ctx_open_(world) != 0) return -1;
    }
    if (!ctx->ic_written) {
        if (ref_write_ic_(world) != 0) return -1;
    }
    return 0;
}

/* Step-begin emit (called from world_step before scheduler advance). */
void k26astro_rt_ref_emit_step_begin(K26AstroWorld *world, double dt)
{
    if (!world || world->mode != K26ASTRO_MODE_REFERENCED) return;
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || !ctx->path) return;
    if (ref_lazy_init_(world) != 0) return;
    ctx->step_idx++;
    K26AstroRefOpHeader hdr = {
        .op_id      = (uint16_t)K26ASTRO_OP_STEP_BEGIN,
        .body_count = 0,
        .step_idx   = ctx->step_idx
    };
    (void)ring_append_(ctx, &hdr, sizeof hdr);
    (void)ring_append_double_(ctx, dt);
}

/* Step-end emit (called from world_step after scheduler advance). */
void k26astro_rt_ref_emit_step_end(K26AstroWorld *world, double world_time_s)
{
    if (!world || world->mode != K26ASTRO_MODE_REFERENCED) return;
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || !ctx->path || !ctx->fp) return;
    K26AstroRefOpHeader hdr = {
        .op_id      = (uint16_t)K26ASTRO_OP_STEP_END,
        .body_count = 0,
        .step_idx   = ctx->step_idx
    };
    (void)ring_append_(ctx, &hdr, sizeof hdr);
    (void)ring_append_double_(ctx, world_time_s);
    /* Flush at step boundaries — guarantees mid-run readability. */
    if (ctx->used > 0) {
        (void)fwrite(ctx->buf, 1, ctx->used, ctx->fp);
        ctx->used = 0;
    }
}

/* Encounter emit (called from encounter.c on pair detect). */
void k26astro_rt_ref_emit_encounter(K26AstroWorld *world,
                                     int i, int j, double k_weight)
{
    if (!world || world->mode != K26ASTRO_MODE_REFERENCED) return;
    K26AstroRefCtx *ctx = (K26AstroRefCtx *)world->ref_ctx;
    if (!ctx || !ctx->path) return;
    if (ref_lazy_init_(world) != 0) return;
    K26AstroRefOpHeader hdr = {
        .op_id      = (uint16_t)K26ASTRO_OP_ENCOUNTER,
        .body_count = 2,
        .step_idx   = ctx->step_idx
    };
    (void)ring_append_(ctx, &hdr, sizeof hdr);
    (void)ring_append_i32_(ctx, (int32_t)i);
    (void)ring_append_i32_(ctx, (int32_t)j);
    (void)ring_append_double_(ctx, k_weight);
}
