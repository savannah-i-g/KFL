#ifndef K26_ATOMIC_IO_H
#define K26_ATOMIC_IO_H

#include <stddef.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Atomic write: stage at "<path>.new", fsync, rename over the target.
 * Leaves no half-written file if power cuts mid-save. Returns 0 on
 * success, -1 on failure (errno set by the underlying syscall).
 *
 * Use for small-to-medium config files; the staging file path is
 * bounded by an internal buffer so very long `path` arguments (>511
 * bytes including the ".new" suffix) are rejected. */
int k26_atomic_write(const char *path, const void *data, size_t len, mode_t mode);

/* Slurp a regular file into a fresh malloc'd buffer. *out points at a
 * NUL-terminated heap string on success (caller frees); *outlen is the
 * byte count NOT counting the trailing NUL. Returns 0 on success, -1 on
 * failure (errno set; *out left NULL).
 *
 * No atomicity guarantee — read is naturally atomic against partial
 * writes only when the writer used k26_atomic_write (rename-over). The
 * symmetric naming is a convention, not a promise. */
int k26_atomic_read(const char *path, char **out, size_t *outlen);

#ifdef __cplusplus
}
#endif

#endif /* K26_ATOMIC_IO_H */
