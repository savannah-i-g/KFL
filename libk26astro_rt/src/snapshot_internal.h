/* snapshot_internal.h — private I/O helpers. */
#ifndef K26ASTRO_RT_SNAPSHOT_INTERNAL_H
#define K26ASTRO_RT_SNAPSHOT_INTERNAL_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Field-by-field LE-encoded write helpers. Writer always emits
 * little-endian; reader byte-swaps if the header's endian probe
 * disagrees. */
int k26astro_snap_write_u16 (FILE *f, uint16_t v);
int k26astro_snap_write_u32 (FILE *f, uint32_t v);
int k26astro_snap_write_u64 (FILE *f, uint64_t v);
int k26astro_snap_write_i64 (FILE *f, int64_t  v);
int k26astro_snap_write_dbl (FILE *f, double   v);
int k26astro_snap_write_buf (FILE *f, const void *p, size_t n);

int k26astro_snap_read_u16  (FILE *f, uint16_t *out, int swap);
int k26astro_snap_read_u32  (FILE *f, uint32_t *out, int swap);
int k26astro_snap_read_u64  (FILE *f, uint64_t *out, int swap);
int k26astro_snap_read_i64  (FILE *f, int64_t  *out, int swap);
int k26astro_snap_read_dbl  (FILE *f, double   *out, int swap);
int k26astro_snap_read_buf  (FILE *f, void *p, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_SNAPSHOT_INTERNAL_H */
