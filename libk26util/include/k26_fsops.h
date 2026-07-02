#ifndef K26_FSOPS_H
#define K26_FSOPS_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* libk26util: filesystem mutation primitives.
 *
 * Pure POSIX (musl-friendly); no glibc-isms like nftw. Usable from
 * any consumer. No FLTK or X dependency, so CLI tools and headless
 * backup utilities link this the same way GUI file managers do.
 *
 * Trash policy: when the host boots from a live medium with persistent
 * overlays, deletes can be routed through a recoverable trash
 * directory; otherwise they fall back to immediate removal.
 * k26_fs_trash_available() returns 1 when ~/.local lives on a
 * non-tmpfs mount (statfs check) and 0 otherwise. Callers should
 * either route deletes through trash_move (which writes a sidecar
 * .trashinfo file recording the original path, time, and mode for
 * subsequent undo) or fall back to confirm + delete_tree when trash
 * is unavailable. */

typedef int (*k26_progress_cb)(const char *current, uint64_t bytes,
                               uint64_t total, void *ud);

/* Recursively delete a path. Symlinks are unlink'd, never followed.
 * Continues on per-entry failures, accumulating a count and recording
 * the first failure path so the caller can surface a single status
 * line. Returns 0 on full success, -1 on any failure (errno reflects
 * the LAST failure; first_fail is more useful for diagnostics).
 *
 * failed_n / first_fail may be NULL if the caller doesn't care.
 * first_fail buffer is left empty if no failures occurred. */
int  k26_fs_delete_tree(const char *path, int *failed_n,
                        char *first_fail, size_t cap);

/* Returns 1 when the trash directory can be created/used (i.e. ~/.local
 * is on a non-tmpfs mount), 0 otherwise. The check is statfs(home) +
 * f_type compare against TMPFS_MAGIC (0x01021994). Caller decides
 * whether to use trash_move or delete_tree based on this. */
int  k26_fs_trash_available(void);

/* Returns the canonical trash root path (~/.local/share/k26/trash) into
 * out, creating the directory hierarchy if it doesn't yet exist. Returns
 * 0 on success, -1 on failure (errno set). out_cap is the buffer size;
 * the path is written with snprintf bounds. */
int  k26_fs_trash_root(char *out, size_t out_cap);

/* Move a path into the trash directory, choosing a unique destination
 * name and writing a sidecar .trashinfo file recording the original
 * path, deletion time, and original mode (consumed by an undo path).
 * The chosen trashed path is written to out_trashed (may be NULL).
 *
 * Returns 0 on success, -1 on failure. Cross-filesystem trash uses
 * rename(2) when possible, falls back to copy+delete on EXDEV (the
 * trash and source must straddle a mount boundary in that case). */
int  k26_fs_trash_move(const char *src, const char *trash_root,
                       char *out_trashed, size_t cap);

/* Walk a tree to compute total size + item count. Used as a pre-pass
 * before copy/move so a progress dialog has totals to render. cancel
 * (may be NULL) is polled between entries; setting *cancel = 1 from
 * another thread aborts the walk. Returns 0 on success, -1 on opendir
 * failure of the root, -2 on cancel. */
int  k26_fs_walk_size(const char *root, off_t *bytes, int *items,
                      volatile int *cancel);

/* Format byte count as a human-readable string ("1.2 MiB"). Returns
 * the number of characters written (excluding NUL), or -1 if cap is
 * too small. cap of 16 is sufficient for any off_t value. */
int  k26_fs_human_size(off_t bytes, char *out, size_t cap);

/* Generate a non-colliding filename inside `dir` based on `base`.
 * If "<dir>/<base>" doesn't exist, writes that. Otherwise tries
 * "<dir>/<base> (copy)", "<dir>/<base> (copy 2)", ... up to a small
 * cap. Returns 0 on success, -1 on cap-exceeded or output overflow. */
int  k26_fs_unique_name(const char *dir, const char *base,
                        char *out, size_t cap);

/* Recursive copy. Symlinks copied as symlinks (readlink + symlink),
 * regular files as O_RDONLY -> O_WRONLY|O_CREAT|O_EXCL with fchmod +
 * utimensat to preserve mode + mtime. Refuses to overwrite existing
 * dst (caller must pre-resolve via k26_fs_unique_name). progress_cb
 * may be NULL. Returns 0 on success, -1 on first failure (first_fail
 * captures the path that failed; errno reflects the failure cause). */
int  k26_fs_copy_tree(const char *src, const char *dst,
                      k26_progress_cb cb, void *ud,
                      char *first_fail, size_t cap);

/* Move src to dst. Tries rename(2) first; on EXDEV falls through to
 * copy_tree + delete_tree. Same return convention as copy_tree. */
int  k26_fs_move(const char *src, const char *dst,
                 k26_progress_cb cb, void *ud,
                 char *first_fail, size_t cap);

#ifdef __cplusplus
}
#endif

#endif /* K26_FSOPS_H */
