/* spk.c — minimal DAF/SPK Type 2 reader + matching synthetic writer.
 *
 * Format reference: NAIF Required Reading on the SPICE Kernel format
 * (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html and
 * daf.html). The salient pieces for us:
 *
 *   - File record (bytes 0..1023): "DAF/SPK " magic, ND/NI counts,
 *     FWARD/BWARD record numbers (1-indexed) of the first/last
 *     summary records, FREE pointer, LOCFMT="LTL-IEEE".
 *
 *   - Summary records form a doubly-linked list (NEXT, PREV in the
 *     first three doubles of each record, NSUM in the third). Each
 *     summary is (ND + ceil(NI/2)) doubles wide: for SPK that's
 *     2 + 3 = 5 doubles = 40 bytes. The NI ints are packed as
 *     int32s in the trailing doubles.
 *
 *   - For SPK: ND=2 (start_et, end_et in TDB seconds past J2000),
 *     NI=6 (target, centre, frame, type, start_addr, end_addr).
 *     Addresses are 1-indexed double-precision word offsets from
 *     the file start.
 *
 *   - Element records (the segment data) sit at the addresses named
 *     in the summaries. Type 2 segment data: N records of RSIZE
 *     doubles, followed by 4 trailing doubles (INIT, INTLEN, RSIZE,
 *     N). Each record holds (MID, RADIUS, x coeffs..., y coeffs...,
 *     z coeffs...). Position at epoch t:
 *       record  = floor((t - INIT) / INTLEN)
 *       s       = (t - record.MID) / record.RADIUS    ∈ [-1, 1]
 *       (x, y, z) = Chebyshev(coeffs, s)
 */
#include "k26astro_ephem/spk.h"
#include "k26astro_ephem/cheby.h"

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define DAF_RECORD_BYTES 1024
#define DAF_DOUBLES_PER_RECORD (DAF_RECORD_BYTES / 8)
#define SPK_ND 2
#define SPK_NI 6
/* Doubles per summary = ND + ceil(NI/2). For SPK: 2 + 3 = 5. */
#define SPK_SUMMARY_DOUBLES (SPK_ND + (SPK_NI + 1) / 2)
#define SPK_SUMMARY_BYTES   (SPK_SUMMARY_DOUBLES * 8)

/* Max summaries that fit in one summary record: 3 header doubles
 * (NEXT, PREV, NSUM) + N * SPK_SUMMARY_DOUBLES doubles.
 * (128 - 3) / 5 = 25 summaries. */
#define SPK_SUMMARIES_PER_RECORD ((DAF_DOUBLES_PER_RECORD - 3) / SPK_SUMMARY_DOUBLES)

#define SPK_MAX_SEGMENTS 128

/* ---- Opaque struct ---------------------------------------------- */
struct K26AstroSpk {
    void                *map;
    size_t               size;
    int                  n_segments;
    K26AstroSpkSegment   segments[SPK_MAX_SEGMENTS];
};

/* ---- Header validation ------------------------------------------ */

static int validate_header_(const uint8_t *base, size_t size,
                            int *out_fward, int *out_bward)
{
    if (size < DAF_RECORD_BYTES) return 1;

    if (memcmp(base, "DAF/SPK ", 8) != 0
     && memcmp(base, "NAIF/DAF", 8) != 0
     && memcmp(base, "DAF/PCK ", 8) != 0) {
        /* Accept the DAF/SPK and the older DAF/NAIF prefixes; reject
         * other DAF subtypes (PCK is included only for diagnostic
         * tolerance; segments will all fail the type check anyway). */
        if (memcmp(base, "DAF/PCK ", 8) != 0) return 2;
    }

    int32_t nd, ni;
    memcpy(&nd, base + 8,  4);
    memcpy(&ni, base + 12, 4);
    if (nd != SPK_ND || ni != SPK_NI) return 3;

    /* LOCFMT at byte 88, 8 chars. Must be "LTL-IEEE". */
    if (memcmp(base + 88, "LTL-IEEE", 8) != 0) return 4;

    int32_t fward, bward;
    memcpy(&fward, base + 76, 4);
    memcpy(&bward, base + 80, 4);
    if (fward < 1 || bward < 1) return 5;
    if (out_fward) *out_fward = fward;
    if (out_bward) *out_bward = bward;
    return 0;
}

/* ---- Summary chain walk ----------------------------------------- */

static const double *record_doubles_(const uint8_t *base, int record_number)
{
    /* 1-indexed record → byte position. */
    return (const double *)(base + (size_t)(record_number - 1) * DAF_RECORD_BYTES);
}

static int parse_segment_meta_(K26AstroSpkSegment *seg,
                                const uint8_t *base,
                                size_t          size,
                                const double   *sum_doubles)
{
    /* sum_doubles points at the 5-double summary. Layout:
     *   [0] start_et (double)
     *   [1] end_et   (double)
     *   [2..4] three doubles holding 6 int32s. */
    seg->start_et = sum_doubles[0];
    seg->end_et   = sum_doubles[1];

    int32_t ints[6];
    memcpy(ints, &sum_doubles[2], 24);
    seg->target_body = ints[0];
    seg->center_body = ints[1];
    seg->frame_id    = ints[2];
    int  type        = ints[3];
    int  start_addr  = ints[4];
    int  end_addr    = ints[5];

    if (type != 2) return 1;             /* not Type 2 — caller skips */
    if (seg->frame_id != 1) return 2;     /* not J2000/ICRF — caller skips */

    /* Trailing 4 doubles of the segment: INIT, INTLEN, RSIZE, N. */
    size_t end_byte = (size_t)(end_addr) * 8;
    if (end_byte > size) return 3;
    const double *seg_doubles = (const double *)(base + (size_t)(start_addr - 1) * 8);
    const double *trailing    = (const double *)(base + end_byte - 32);

    seg->init_et             = trailing[0];
    seg->interval_seconds    = trailing[1];
    seg->record_size_doubles = (int)trailing[2];
    seg->n_records           = (int)trailing[3];
    if (seg->record_size_doubles < 4 || seg->n_records < 1) return 4;
    if ((seg->record_size_doubles - 2) % 3 != 0) return 5;
    seg->coeffs_per_axis = (seg->record_size_doubles - 2) / 3;

    /* Sanity: data range covers what the summary claims. */
    size_t expected_bytes = (size_t)seg->n_records
                          * (size_t)seg->record_size_doubles * 8
                          + 32;   /* + 4 trailing doubles */
    size_t actual_bytes = (size_t)(end_addr - start_addr + 1) * 8;
    if (actual_bytes < expected_bytes) return 6;

    seg->records = seg_doubles;
    return 0;
}

K26AstroSpk *k26astro_spk_open(const char *path)
{
    if (!path) return NULL;
    int fd = open(path, O_RDONLY);
    if (fd < 0) return NULL;
    struct stat st;
    if (fstat(fd, &st) != 0) { close(fd); return NULL; }
    if (st.st_size < DAF_RECORD_BYTES) { close(fd); return NULL; }

    void *map = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);   /* mmap retains the reference */
    if (map == MAP_FAILED) return NULL;

    int fward, bward;
    if (validate_header_((const uint8_t *)map, (size_t)st.st_size,
                         &fward, &bward) != 0) {
        munmap(map, (size_t)st.st_size);
        return NULL;
    }

    K26AstroSpk *spk = (K26AstroSpk *)calloc(1, sizeof(*spk));
    if (!spk) { munmap(map, (size_t)st.st_size); return NULL; }
    spk->map  = map;
    spk->size = (size_t)st.st_size;

    /* Walk summary chain FWARD → ... → BWARD. */
    const uint8_t *base = (const uint8_t *)map;
    int rec = fward;
    while (rec > 0 && spk->n_segments < SPK_MAX_SEGMENTS) {
        const double *rd = record_doubles_(base, rec);
        int next  = (int)rd[0];
        int nsum  = (int)rd[2];
        (void)rd;   /* PREV at rd[1] unused */

        if (nsum < 0 || nsum > SPK_SUMMARIES_PER_RECORD) break;

        for (int i = 0; i < nsum && spk->n_segments < SPK_MAX_SEGMENTS; i++) {
            const double *sum = rd + 3 + (size_t)i * SPK_SUMMARY_DOUBLES;
            K26AstroSpkSegment seg = {0};
            int rc = parse_segment_meta_(&seg, base, spk->size, sum);
            if (rc == 0) {
                spk->segments[spk->n_segments++] = seg;
            }
            /* rc != 0 → silently skip; types other than 2 and frames
             * other than J2000 are not errors, they're just outside
             * v0.1's scope. */
        }

        if (next == 0 || next == rec) break;
        rec = next;
    }

    return spk;
}

void k26astro_spk_close(K26AstroSpk *spk)
{
    if (!spk) return;
    if (spk->map) munmap(spk->map, spk->size);
    free(spk);
}

int k26astro_spk_n_segments(const K26AstroSpk *spk)
{
    return spk ? spk->n_segments : 0;
}

const K26AstroSpkSegment *k26astro_spk_segment(const K26AstroSpk *spk, int idx)
{
    if (!spk || idx < 0 || idx >= spk->n_segments) return NULL;
    return &spk->segments[idx];
}

const K26AstroSpkSegment *
k26astro_spk_find_segment(const K26AstroSpk *spk, int target_body, double et)
{
    if (!spk) return NULL;
    for (int i = 0; i < spk->n_segments; i++) {
        const K26AstroSpkSegment *s = &spk->segments[i];
        if (s->target_body != target_body) continue;
        if (et < s->start_et || et > s->end_et) continue;
        return s;
    }
    return NULL;
}

/* ---- Evaluation ------------------------------------------------- */

static int eval_segment_(const K26AstroSpkSegment *s, double et,
                         double out_xyz[3], double out_xyz_dot[3])
{
    int idx = (int)((et - s->init_et) / s->interval_seconds);
    if (idx < 0) idx = 0;
    if (idx >= s->n_records) idx = s->n_records - 1;
    const double *rec = s->records + (size_t)idx * (size_t)s->record_size_doubles;
    double mid    = rec[0];
    double radius = rec[1];
    if (radius == 0.0) return 1;
    double s_norm = (et - mid) / radius;
    /* Clamp ε past boundary in case of integer-record-index rounding
     * at the very edge of the interval; the Chebyshev domain is
     * mathematically [-1, 1] and going slightly outside damages
     * numerical stability. */
    if (s_norm >  1.0) s_norm =  1.0;
    if (s_norm < -1.0) s_norm = -1.0;

    int K = s->coeffs_per_axis;
    const double *cx = rec + 2;
    const double *cy = cx + K;
    const double *cz = cy + K;

    if (out_xyz_dot) {
        double vx, dx, vy, dy, vz, dz;
        k26astro_cheby_eval_both(cx, K, s_norm, &vx, &dx);
        k26astro_cheby_eval_both(cy, K, s_norm, &vy, &dy);
        k26astro_cheby_eval_both(cz, K, s_norm, &vz, &dz);
        out_xyz[0] = vx; out_xyz[1] = vy; out_xyz[2] = vz;
        /* Chain ds/dt = 1/radius. */
        double inv_r = 1.0 / radius;
        out_xyz_dot[0] = dx * inv_r;
        out_xyz_dot[1] = dy * inv_r;
        out_xyz_dot[2] = dz * inv_r;
    } else {
        out_xyz[0] = k26astro_cheby_eval(cx, K, s_norm);
        out_xyz[1] = k26astro_cheby_eval(cy, K, s_norm);
        out_xyz[2] = k26astro_cheby_eval(cz, K, s_norm);
    }
    return 0;
}

int k26astro_spk_pos(const K26AstroSpk *spk,
                     int target_body, double et,
                     double out_xyz[3])
{
    const K26AstroSpkSegment *s = k26astro_spk_find_segment(spk, target_body, et);
    if (!s) return 1;
    return eval_segment_(s, et, out_xyz, NULL);
}

int k26astro_spk_pos_vel(const K26AstroSpk *spk,
                         int target_body, double et,
                         double out_xyz[3], double out_xyz_dot[3])
{
    const K26AstroSpkSegment *s = k26astro_spk_find_segment(spk, target_body, et);
    if (!s) return 1;
    return eval_segment_(s, et, out_xyz, out_xyz_dot);
}

/* ---- Synthetic writer ------------------------------------------- *
 *
 * Produces a valid DAF/SPK file with the requested Type-2 segments.
 * Layout:
 *   record 1     : file record (header)
 *   record 2     : single summary record (NEXT=0, PREV=0)
 *   record 3     : name record (40 bytes per summary, zero-padded
 *                  to fill the 1024-byte record)
 *   record 4..N  : element data — segments packed sequentially, each
 *                  ending with the 4-double trailer.
 *
 * Segments don't share records; each segment's data starts at the
 * next 1024-byte record boundary. This wastes some space (up to
 * 1023 bytes per segment) but keeps the writer trivial and the
 * resulting addresses pleasingly aligned.
 */
int k26astro_spk_write_synthetic(const char *path,
                                  const K26AstroSpkWriteSegment *segs,
                                  int n_segs)
{
    if (!path || !segs || n_segs <= 0) return 1;
    if (n_segs > (int)SPK_SUMMARIES_PER_RECORD) return 2;

    FILE *f = fopen(path, "wb");
    if (!f) return 3;

    /* Build the file record. */
    uint8_t file_rec[DAF_RECORD_BYTES] = {0};
    memcpy(file_rec, "DAF/SPK ", 8);
    int32_t nd = SPK_ND, ni = SPK_NI;
    memcpy(file_rec + 8,  &nd, 4);
    memcpy(file_rec + 12, &ni, 4);
    /* Internal filename — fill with spaces. */
    memset(file_rec + 16, ' ', 60);
    int32_t fward = 2, bward = 2;
    memcpy(file_rec + 76, &fward, 4);
    memcpy(file_rec + 80, &bward, 4);
    /* FREE filled in after we know all addresses. Leave as 0 for now. */
    memcpy(file_rec + 88, "LTL-IEEE", 8);

    /* Build the summary record. Layout: 3 header doubles + summaries. */
    double sum_rec[DAF_DOUBLES_PER_RECORD] = {0};
    sum_rec[0] = 0.0;            /* NEXT — none */
    sum_rec[1] = 0.0;            /* PREV — none */
    sum_rec[2] = (double)n_segs; /* NSUM */

    /* We'll compute start_addr for each segment after writing the
     * file/summary/name records, since each segment starts at the
     * next record boundary. */
    int32_t cur_record = 4;   /* first segment-data record number */
    int32_t seg_start_addr[SPK_SUMMARIES_PER_RECORD] = {0};
    int32_t seg_end_addr  [SPK_SUMMARIES_PER_RECORD] = {0};

    for (int i = 0; i < n_segs; i++) {
        if (segs[i].coeffs_per_axis < 1 || segs[i].n_records < 1) {
            fclose(f); return 4;
        }
        int32_t rsize = 2 + 3 * segs[i].coeffs_per_axis;
        int32_t total_doubles = segs[i].n_records * rsize + 4;   /* + INIT/INTLEN/RSIZE/N */

        int32_t start_addr = (cur_record - 1) * DAF_DOUBLES_PER_RECORD + 1;
        int32_t end_addr   = start_addr + total_doubles - 1;
        seg_start_addr[i] = start_addr;
        seg_end_addr  [i] = end_addr;

        /* Advance cur_record past this segment's data, rounded up to
         * the next record boundary. */
        int32_t total_records = (total_doubles + DAF_DOUBLES_PER_RECORD - 1)
                              / DAF_DOUBLES_PER_RECORD;
        cur_record += total_records;
    }

    /* Now fill the summary doubles. */
    for (int i = 0; i < n_segs; i++) {
        const K26AstroSpkWriteSegment *s = &segs[i];
        double *sm = &sum_rec[3 + (size_t)i * SPK_SUMMARY_DOUBLES];
        sm[0] = s->start_et;
        sm[1] = s->end_et;
        int32_t ints[6] = {
            s->target_body, s->center_body, 1, 2,
            seg_start_addr[i], seg_end_addr[i]
        };
        memcpy(&sm[2], ints, 24);
    }

    /* FREE = first byte past all segment data. */
    int32_t free_addr = (cur_record - 1) * DAF_DOUBLES_PER_RECORD + 1;
    memcpy(file_rec + 84, &free_addr, 4);

    /* Write file record. */
    if (fwrite(file_rec, 1, DAF_RECORD_BYTES, f) != DAF_RECORD_BYTES) {
        fclose(f); return 5;
    }
    /* Write summary record. */
    if (fwrite(sum_rec, 8, DAF_DOUBLES_PER_RECORD, f) != DAF_DOUBLES_PER_RECORD) {
        fclose(f); return 5;
    }
    /* Write name record (zeros). */
    uint8_t name_rec[DAF_RECORD_BYTES] = {0};
    if (fwrite(name_rec, 1, DAF_RECORD_BYTES, f) != DAF_RECORD_BYTES) {
        fclose(f); return 5;
    }

    /* Now write each segment's data, padded to a record boundary. */
    for (int i = 0; i < n_segs; i++) {
        const K26AstroSpkWriteSegment *s = &segs[i];
        int32_t rsize = 2 + 3 * s->coeffs_per_axis;
        size_t records_bytes = (size_t)s->n_records * (size_t)rsize * 8;

        /* Sanity check that the file position matches the address we
         * promised in the summary. */
        long now = ftell(f);
        long expected = (long)(seg_start_addr[i] - 1) * 8;
        if (now != expected) {
            /* Pad if needed. */
            if (now < expected) {
                uint8_t zero = 0;
                for (long k = now; k < expected; k++) fputc(zero, f);
            } else {
                fclose(f); return 6;
            }
        }

        /* Records. */
        if (fwrite(s->records, 8, (size_t)s->n_records * (size_t)rsize, f)
            != (size_t)s->n_records * (size_t)rsize) {
            fclose(f); return 5;
        }
        /* Trailer. */
        double trailer[4] = {
            s->start_et, s->interval_seconds, (double)rsize, (double)s->n_records
        };
        if (fwrite(trailer, 8, 4, f) != 4) { fclose(f); return 5; }

        /* Pad to record boundary. */
        long cur = ftell(f);
        long boundary = ((cur + DAF_RECORD_BYTES - 1) / DAF_RECORD_BYTES) * DAF_RECORD_BYTES;
        for (long k = cur; k < boundary; k++) fputc(0, f);

        (void)records_bytes;
    }

    fclose(f);
    return 0;
}
