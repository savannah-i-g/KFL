/* k26astro_replay — REFERENCED-mode log inspector + comparator.
 *
 * The op-log written by libk26astro_rt's referenced mode is packed
 * binary; this tool provides two functions:
 *
 *   --check log1 log2      byte-for-byte comparison; reports first
 *                          divergent record + offset
 *   --dump log [--sample N] human-readable dump of every Nth record
 *
 * Determinism contract: two runs of the same simulation with
 * identical initial conditions must produce byte-identical logs.
 * `--check` is the regression net; `--dump` is for triage when
 * --check finds a divergence.
 *
 * Host-only tool. Links against libk26astro_rt for the schema
 * (op enum, arity table, magic constants). */
#define _GNU_SOURCE
#include "k26astro_rt/referenced.h"
#include "k26astro_body/body.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *op_name_(uint16_t op_id)
{
    switch (op_id) {
        case K26ASTRO_OP_IC_BLOCK:   return "IC";
        case K26ASTRO_OP_STEP_BEGIN: return "STEP_BEGIN";
        case K26ASTRO_OP_STEP_END:   return "STEP_END";
        case K26ASTRO_OP_ENCOUNTER:  return "ENCOUNTER";
        default:                     return "?";
    }
}

/* Slurp a whole file into memory. Returns 0/-1; fills *out_buf, *out_len. */
static int slurp_(const char *path, uint8_t **out_buf, size_t *out_len)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "open %s: %s\n", path, strerror(errno));
        return -1;
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    if (sz < 0) { fclose(f); return -1; }
    fseek(f, 0, SEEK_SET);
    uint8_t *buf = (uint8_t *)malloc((size_t)sz);
    if (!buf) { fclose(f); return -1; }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        fclose(f); free(buf); return -1;
    }
    fclose(f);
    *out_buf = buf;
    *out_len = (size_t)sz;
    return 0;
}

static int verify_header_(const uint8_t *buf, size_t len, size_t *out_off)
{
    if (len < K26ASTRO_REF_MAGIC_LEN + 8) {
        fprintf(stderr, "log truncated (header)\n");
        return -1;
    }
    if (memcmp(buf, K26ASTRO_REF_MAGIC, K26ASTRO_REF_MAGIC_LEN) != 0) {
        fprintf(stderr, "magic mismatch (not a K26REF log)\n");
        return -1;
    }
    uint32_t version;
    memcpy(&version, buf + K26ASTRO_REF_MAGIC_LEN, sizeof version);
    if (version != K26ASTRO_REF_VERSION) {
        fprintf(stderr, "version %u != expected %u\n",
                version, (unsigned)K26ASTRO_REF_VERSION);
        return -1;
    }
    *out_off = K26ASTRO_REF_MAGIC_LEN + 8;
    return 0;
}

static int cmd_check_(const char *p1, const char *p2)
{
    uint8_t *a = NULL, *b = NULL;
    size_t la = 0, lb = 0;
    if (slurp_(p1, &a, &la) != 0) return 2;
    if (slurp_(p2, &b, &lb) != 0) { free(a); return 2; }
    if (la != lb) {
        fprintf(stderr,
                "DIVERGENT: log sizes differ — %s=%zu bytes vs %s=%zu bytes\n",
                p1, la, p2, lb);
        free(a); free(b);
        return 1;
    }
    /* Find first divergent byte. */
    for (size_t i = 0; i < la; i++) {
        if (a[i] != b[i]) {
            fprintf(stderr,
                "DIVERGENT: first byte differs at offset %zu — %s=0x%02x vs %s=0x%02x\n",
                i, p1, a[i], p2, b[i]);
            free(a); free(b);
            return 1;
        }
    }
    printf("identical (%zu bytes)\n", la);
    free(a); free(b);
    return 0;
}

static int cmd_dump_(const char *path, int sample)
{
    uint8_t *buf = NULL;
    size_t len = 0;
    if (slurp_(path, &buf, &len) != 0) return 2;
    size_t off = 0;
    if (verify_header_(buf, len, &off) != 0) { free(buf); return 1; }
    printf("# K26REF v%u log: %s (%zu bytes)\n",
           (unsigned)K26ASTRO_REF_VERSION, path, len);

    int rec_idx = 0;
    while (off < len) {
        if (len - off < sizeof(K26AstroRefOpHeader)) {
            fprintf(stderr,
                    "truncated record at offset %zu (need %zu, have %zu)\n",
                    off, sizeof(K26AstroRefOpHeader), len - off);
            free(buf); return 1;
        }
        K26AstroRefOpHeader hdr;
        memcpy(&hdr, buf + off, sizeof hdr);
        off += sizeof hdr;

        size_t payload_bytes;
        if (hdr.op_id == K26ASTRO_OP_IC_BLOCK) {
            /* IC: body_count bodies, each name[32] + 8 doubles. */
            payload_bytes = (size_t)hdr.body_count
                          * (K26ASTRO_BODY_NAME_MAX + 8 * sizeof(uint64_t));
        } else {
            payload_bytes = (size_t)hdr.body_count * sizeof(int32_t);
            if (hdr.op_id < k26astro_ref_op_count) {
                payload_bytes += (size_t)k26astro_ref_op_arity[hdr.op_id]
                               * sizeof(uint64_t);
            }
        }
        if (off + payload_bytes > len) {
            fprintf(stderr, "truncated payload at offset %zu\n", off);
            free(buf); return 1;
        }

        if (sample <= 1 || (rec_idx % sample) == 0) {
            printf("[%d] off=%zu op=%-12s body_count=%u step=%u",
                   rec_idx, off - sizeof hdr, op_name_(hdr.op_id),
                   (unsigned)hdr.body_count, (unsigned)hdr.step_idx);
            if (hdr.op_id == K26ASTRO_OP_IC_BLOCK) {
                printf("  (IC payload %zu bytes)\n", payload_bytes);
            } else {
                /* Print body indices + double args. */
                size_t cur = off;
                for (int i = 0; i < hdr.body_count; i++) {
                    int32_t v;
                    memcpy(&v, buf + cur, sizeof v);
                    cur += sizeof v;
                    printf(" b[%d]=%d", i, v);
                }
                int arity = (hdr.op_id < k26astro_ref_op_count)
                              ? k26astro_ref_op_arity[hdr.op_id] : 0;
                for (int i = 0; i < arity; i++) {
                    union { uint64_t u; double d; } cvt;
                    memcpy(&cvt.u, buf + cur, sizeof cvt.u);
                    cur += sizeof cvt.u;
                    printf(" a[%d]=%a", i, cvt.d);
                }
                printf("\n");
            }
        }
        off += payload_bytes;
        rec_idx++;
    }
    printf("# %d records total\n", rec_idx);
    free(buf);
    return 0;
}

static void usage_(void)
{
    fprintf(stderr,
        "usage: k26astro_replay --check log1 log2\n"
        "       k26astro_replay --dump  log  [--sample N]\n");
}

int main(int argc, char **argv)
{
    if (argc < 2) { usage_(); return 2; }
    if (strcmp(argv[1], "--check") == 0) {
        if (argc != 4) { usage_(); return 2; }
        return cmd_check_(argv[2], argv[3]);
    }
    if (strcmp(argv[1], "--dump") == 0) {
        if (argc < 3) { usage_(); return 2; }
        int sample = 1;
        for (int i = 3; i < argc; i++) {
            if (strcmp(argv[i], "--sample") == 0 && i + 1 < argc) {
                sample = atoi(argv[i + 1]);
                if (sample < 1) sample = 1;
                i++;
            }
        }
        return cmd_dump_(argv[2], sample);
    }
    usage_();
    return 2;
}
