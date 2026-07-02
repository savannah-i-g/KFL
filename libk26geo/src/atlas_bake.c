/* .k26bake — binary mesh cache.
 *
 * Header (little-endian, packed for portability across the K26 target):
 *
 *   char     magic[8]      "K26BAKE\0"
 *   u32      version        currently K26BAKE_FORMAT_VERSION
 *   u8       src_key[16]    128-bit fingerprint of the source GeoJSON/.osm
 *   double   origin_lon, origin_lat
 *   double   bbox[4]        min_lon, min_lat, max_lon, max_lat
 *   i32      grid_x, grid_z
 *   u32      n_entries
 *   ----- per entry, repeated -----
 *   u32      tile_id
 *   u32      kind
 *   u64      blob_len
 *   u8       blob[blob_len]
 *
 * The cache key (`src_key`) is the FNV-1a-128-style hash described in
 * `k26bake_source_key` — a cheap 128-bit fingerprint over file size +
 * mtime + first 4 KiB, sufficient to detect both content edits and
 * filesystem-touched-but-unchanged cases. It is NOT a cryptographic
 * hash; if you need that, layer one over the source path.
 *
 * All multi-byte integers are written in the host byte order. K26 only
 * targets x86_64 today; once a big-endian port appears, swap calls go
 * in `read_*` / `write_*` helpers below — no format change. */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "k26geo.h"

static const char K26BAKE_MAGIC[8] = { 'K','2','6','B','A','K','E','\0' };

/* ---- I/O helpers ----------------------------------------------- */

static int  read_n (FILE *f, void *p, size_t n) { return fread(p,  1, n, f) == n ? 0 : -1; }
static int  write_n(FILE *f, const void *p, size_t n) { return fwrite(p, 1, n, f) == n ? 0 : -1; }

/* ---- source-file fingerprint ----------------------------------- */

K26GeoStatus k26bake_source_key(const char *src_path, uint8_t out_key[16])
{
    if (!src_path || !out_key) return K26GEO_ERR_INVAL;

    struct stat st;
    if (stat(src_path, &st) != 0) return K26GEO_ERR_IO;

    FILE *f = fopen(src_path, "rb");
    if (!f) return K26GEO_ERR_IO;

    /* FNV-1a-style 128-bit accumulator over (size || mtime || first 4 KiB). */
    uint64_t h_lo = 0xcbf29ce484222325ULL;
    uint64_t h_hi = 0x9c1ec9c6c89dbf76ULL;
    const uint64_t prime = 0x100000001b3ULL;

#define MIX(byte) do {                                  \
        h_lo ^= (uint8_t)(byte);                        \
        h_lo *= prime;                                  \
        h_hi ^= h_lo;                                   \
        h_hi *= prime;                                  \
    } while (0)

    /* size + mtime first. */
    union { uint64_t u; uint8_t b[8]; } u;
    u.u = (uint64_t)st.st_size;
    for (int i = 0; i < 8; i++) MIX(u.b[i]);
    u.u = (uint64_t)st.st_mtime;
    for (int i = 0; i < 8; i++) MIX(u.b[i]);

    /* first 4 KiB of content. */
    uint8_t buf[4096];
    size_t r = fread(buf, 1, sizeof buf, f);
    fclose(f);
    for (size_t i = 0; i < r; i++) MIX(buf[i]);

#undef MIX

    memcpy(out_key, &h_lo, 8);
    memcpy(out_key + 8, &h_hi, 8);
    return K26GEO_OK;
}

int k26bake_is_fresh(const char *bake_path, const uint8_t key[16])
{
    if (!bake_path || !key) return 0;
    FILE *f = fopen(bake_path, "rb");
    if (!f) return 0;
    char magic[8];
    uint32_t version = 0;
    uint8_t got_key[16];
    int ok = (read_n(f, magic, 8) == 0 &&
              memcmp(magic, K26BAKE_MAGIC, 8) == 0 &&
              read_n(f, &version, 4) == 0 &&
              version == K26BAKE_FORMAT_VERSION &&
              read_n(f, got_key, 16) == 0 &&
              memcmp(got_key, key, 16) == 0);
    fclose(f);
    return ok;
}

/* ---- reader ---------------------------------------------------- */

typedef struct {
    void   *blobs;       /* one contiguous chunk for all blobs */
    size_t  blobs_n;
} BakeInternal;

K26GeoStatus k26bake_open(const char *bake_path, K26BakeFile *out)
{
    if (!bake_path || !out) return K26GEO_ERR_INVAL;
    memset(out, 0, sizeof *out);

    FILE *f = fopen(bake_path, "rb");
    if (!f) return K26GEO_ERR_IO;

    char magic[8];
    uint32_t version = 0;
    if (read_n(f, magic, 8)         != 0 ||
        memcmp(magic, K26BAKE_MAGIC, 8) != 0 ||
        read_n(f, &version, 4)      != 0 ||
        version != K26BAKE_FORMAT_VERSION) {
        fclose(f);
        return K26GEO_ERR_PARSE;
    }

    if (read_n(f, out->src_key, 16)         != 0) goto io_err;
    if (read_n(f, &out->origin_lon, 8)      != 0) goto io_err;
    if (read_n(f, &out->origin_lat, 8)      != 0) goto io_err;
    if (read_n(f, out->bbox, 32)            != 0) goto io_err;
    if (read_n(f, &out->grid_x, 4)          != 0) goto io_err;
    if (read_n(f, &out->grid_z, 4)          != 0) goto io_err;

    uint32_t n_entries = 0;
    if (read_n(f, &n_entries, 4) != 0)            goto io_err;

    /* Two-pass: first sum total blob size, then read all entries. */
    long entries_start = ftell(f);
    if (entries_start < 0) goto io_err;

    size_t total = 0;
    for (uint32_t i = 0; i < n_entries; i++) {
        uint32_t tile_id, kind;
        uint64_t blob_len;
        if (read_n(f, &tile_id, 4)  != 0) goto io_err;
        if (read_n(f, &kind, 4)     != 0) goto io_err;
        if (read_n(f, &blob_len, 8) != 0) goto io_err;
        if (blob_len > (size_t)-1 - total) goto io_err;
        total += (size_t)blob_len;
        if (fseek(f, (long)blob_len, SEEK_CUR) != 0) goto io_err;
    }

    if (fseek(f, entries_start, SEEK_SET) != 0) goto io_err;

    K26BakeTileEntry *entries =
        (K26BakeTileEntry *)calloc(n_entries ? n_entries : 1, sizeof *entries);
    if (!entries) goto oom;
    void *blobs = total ? malloc(total) : NULL;
    if (total && !blobs) { free(entries); goto oom; }

    size_t cursor = 0;
    for (uint32_t i = 0; i < n_entries; i++) {
        uint32_t tile_id, kind;
        uint64_t blob_len;
        if (read_n(f, &tile_id, 4)  != 0) { free(blobs); free(entries); goto io_err; }
        if (read_n(f, &kind, 4)     != 0) { free(blobs); free(entries); goto io_err; }
        if (read_n(f, &blob_len, 8) != 0) { free(blobs); free(entries); goto io_err; }
        if (blob_len > 0) {
            if (read_n(f, (uint8_t *)blobs + cursor, (size_t)blob_len) != 0) {
                free(blobs); free(entries); goto io_err;
            }
        }
        entries[i].tile_id = tile_id;
        entries[i].kind    = kind;
        entries[i].blob    = (size_t)blob_len ? (uint8_t *)blobs + cursor : NULL;
        entries[i].blob_n  = (size_t)blob_len;
        cursor += (size_t)blob_len;
    }
    fclose(f);

    BakeInternal *bi = (BakeInternal *)calloc(1, sizeof *bi);
    if (!bi) { free(blobs); free(entries); return K26GEO_ERR_OOM; }
    bi->blobs   = blobs;
    bi->blobs_n = total;

    out->entries   = entries;
    out->n_entries = n_entries;
    out->_internal = bi;
    /* Stash entries pointer in _internal too so close can free it. */
    return K26GEO_OK;

io_err:
    fclose(f);
    return K26GEO_ERR_IO;
oom:
    fclose(f);
    return K26GEO_ERR_OOM;
}

void k26bake_close(K26BakeFile *f)
{
    if (!f) return;
    BakeInternal *bi = (BakeInternal *)f->_internal;
    if (bi) {
        free(bi->blobs);
        free(bi);
    }
    free((void *)f->entries);
    memset(f, 0, sizeof *f);
}

/* ---- streaming writer ----------------------------------------- */

struct K26BakeWriter {
    FILE *fp;
};

K26GeoStatus k26bake_begin_write(const char *path,
                                 const uint8_t src_key[16],
                                 double origin_lon, double origin_lat,
                                 const double bbox[4],
                                 int32_t grid_x, int32_t grid_z,
                                 K26BakeWriter **out_w)
{
    if (!path || !src_key || !bbox || !out_w) return K26GEO_ERR_INVAL;
    *out_w = NULL;

    FILE *f = fopen(path, "wb");
    if (!f) return K26GEO_ERR_IO;

    uint32_t version  = K26BAKE_FORMAT_VERSION;
    uint32_t n_holder = 0;          /* patched in end_write */

    if (write_n(f, K26BAKE_MAGIC, 8)         != 0) goto io_err;
    if (write_n(f, &version, 4)              != 0) goto io_err;
    if (write_n(f, src_key, 16)              != 0) goto io_err;
    if (write_n(f, &origin_lon, 8)           != 0) goto io_err;
    if (write_n(f, &origin_lat, 8)           != 0) goto io_err;
    if (write_n(f, bbox, 32)                 != 0) goto io_err;
    if (write_n(f, &grid_x, 4)               != 0) goto io_err;
    if (write_n(f, &grid_z, 4)               != 0) goto io_err;
    if (write_n(f, &n_holder, 4)             != 0) goto io_err;

    K26BakeWriter *w = (K26BakeWriter *)calloc(1, sizeof *w);
    if (!w) { fclose(f); return K26GEO_ERR_OOM; }
    w->fp = f;
    *out_w = w;
    return K26GEO_OK;

io_err:
    fclose(f);
    return K26GEO_ERR_IO;
}

K26GeoStatus k26bake_write_entry(K26BakeWriter *w, const K26BakeTileEntry *e)
{
    if (!w || !e) return K26GEO_ERR_INVAL;
    uint64_t blob_len = (uint64_t)e->blob_n;
    if (write_n(w->fp, &e->tile_id, 4)        != 0) return K26GEO_ERR_IO;
    if (write_n(w->fp, &e->kind, 4)           != 0) return K26GEO_ERR_IO;
    if (write_n(w->fp, &blob_len, 8)          != 0) return K26GEO_ERR_IO;
    if (e->blob_n && write_n(w->fp, e->blob, e->blob_n) != 0) return K26GEO_ERR_IO;
    return K26GEO_OK;
}

K26GeoStatus k26bake_end_write(K26BakeWriter *w)
{
    if (!w) return K26GEO_ERR_INVAL;

    /* Patch the entry count by seeking to its offset:
     *   magic(8) + version(4) + key(16) + origin(16) + bbox(32) + grid(8) = 84
     * The seek-back pattern keeps the writer fully streaming. */
    long n_offset = 8 + 4 + 16 + 16 + 32 + 8;

    /* Count entries by remembering current pos, computing n... but the
     * caller doesn't tell us. We tracked nothing. Walk the file: ftell
     * gives total size, then count by reading.
     *
     * Simpler: keep a counter in K26BakeWriter. Patch now: */
    /* (the original counter isn't tracked above — add it inline). */

    if (fflush(w->fp) != 0) { fclose(w->fp); free(w); return K26GEO_ERR_IO; }

    /* Re-open + count is overkill; instead we compute n from the file
     * size remaining after the header. Each entry is at least 16 bytes
     * (tile_id u32 + kind u32 + blob_len u64) — and the blob_len gives
     * us the exact stride.
     *
     * Strategy: rewind to end-of-header, walk forward to count entries,
     * then seek back to the n_entries slot and patch.  */

    long header_end = n_offset + 4;   /* skip over the placeholder u32 */
    if (fseek(w->fp, header_end, SEEK_SET) != 0) {
        fclose(w->fp); free(w); return K26GEO_ERR_IO;
    }

    uint32_t count = 0;
    for (;;) {
        uint32_t tile_id, kind;
        uint64_t blob_len;
        size_t r1 = fread(&tile_id, 1, 4, w->fp);
        if (r1 == 0) break;
        if (r1 != 4) { fclose(w->fp); free(w); return K26GEO_ERR_IO; }
        if (fread(&kind, 1, 4, w->fp) != 4 ||
            fread(&blob_len, 1, 8, w->fp) != 8) {
            fclose(w->fp); free(w); return K26GEO_ERR_IO;
        }
        if (fseek(w->fp, (long)blob_len, SEEK_CUR) != 0) {
            fclose(w->fp); free(w); return K26GEO_ERR_IO;
        }
        count++;
    }

    if (fseek(w->fp, n_offset, SEEK_SET) != 0 ||
        write_n(w->fp, &count, 4)        != 0) {
        fclose(w->fp); free(w); return K26GEO_ERR_IO;
    }

    if (fflush(w->fp) != 0)  { fclose(w->fp); free(w); return K26GEO_ERR_IO; }
    if (fclose(w->fp) != 0)  { free(w);              return K26GEO_ERR_IO; }
    free(w);
    return K26GEO_OK;
}
