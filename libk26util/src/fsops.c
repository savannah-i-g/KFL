/* libk26util: filesystem mutation primitives.
 *
 * Pure POSIX. opendir/readdir/lstat recursion, never nftw. lstat
 * everywhere so we never follow symlinks into a target tree. */

#include "k26_fsops.h"

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/statfs.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/* musl + glibc both define these constants but only inside <linux/magic.h>
 * which we don't want to drag in. The value is stable kernel ABI. */
#ifndef TMPFS_MAGIC
#define TMPFS_MAGIC 0x01021994
#endif

/* ---- helpers ----------------------------------------------------- */

static const char *fsops_home(void)
{
    const char *h = getenv("HOME");
    if (h && *h) return h;
    struct passwd *pw = getpwuid(getuid());
    return pw ? pw->pw_dir : "/root";
}

static int mkdir_p(const char *path)
{
    if (!path || !*path) return -1;
    char buf[1024];
    size_t len = strlen(path);
    if (len >= sizeof(buf)) { errno = ENAMETOOLONG; return -1; }
    memcpy(buf, path, len + 1);
    for (size_t i = 1; i <= len; i++) {
        if (buf[i] == '/' || buf[i] == '\0') {
            char saved = buf[i];
            buf[i] = '\0';
            if (mkdir(buf, 0755) != 0 && errno != EEXIST) return -1;
            buf[i] = saved;
        }
    }
    return 0;
}

/* ---- delete_tree ------------------------------------------------- */

static void delete_tree_rec(const char *path, int *failed_n,
                            char *first_fail, size_t cap)
{
    struct stat st;
    if (lstat(path, &st) != 0) {
        if (failed_n) (*failed_n)++;
        if (first_fail && cap && first_fail[0] == '\0')
            snprintf(first_fail, cap, "%s", path);
        return;
    }
    if (S_ISDIR(st.st_mode)) {
        DIR *d = opendir(path);
        if (!d) {
            if (failed_n) (*failed_n)++;
            if (first_fail && cap && first_fail[0] == '\0')
                snprintf(first_fail, cap, "%s", path);
            return;
        }
        struct dirent *de;
        while ((de = readdir(d)) != NULL) {
            if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, "..")) continue;
            char child[2048];
            snprintf(child, sizeof(child), "%s/%s", path, de->d_name);
            delete_tree_rec(child, failed_n, first_fail, cap);
        }
        closedir(d);
        if (rmdir(path) != 0) {
            if (failed_n) (*failed_n)++;
            if (first_fail && cap && first_fail[0] == '\0')
                snprintf(first_fail, cap, "%s", path);
        }
    } else {
        /* Regular file, symlink, special — unlink unconditionally
         * (symlinks: never follow into the target). */
        if (unlink(path) != 0) {
            if (failed_n) (*failed_n)++;
            if (first_fail && cap && first_fail[0] == '\0')
                snprintf(first_fail, cap, "%s", path);
        }
    }
}

int k26_fs_delete_tree(const char *path, int *failed_n,
                       char *first_fail, size_t cap)
{
    if (!path || !*path) return -1;
    int local_failed = 0;
    int *fp = failed_n ? failed_n : &local_failed;
    *fp = 0;
    if (first_fail && cap) first_fail[0] = '\0';

    delete_tree_rec(path, fp, first_fail, cap);
    return *fp == 0 ? 0 : -1;
}

/* ---- trash ------------------------------------------------------- */

int k26_fs_trash_available(void)
{
    const char *home = fsops_home();
    char local[1024];
    snprintf(local, sizeof(local), "%s/.local", home);

    /* If ~/.local doesn't exist yet, statfs the home itself — it's the
     * mount we'd put trash under regardless. */
    struct statfs sf;
    const char *probe = local;
    struct stat st;
    if (stat(local, &st) != 0) probe = home;
    if (statfs(probe, &sf) != 0) return 0;
    return sf.f_type == TMPFS_MAGIC ? 0 : 1;
}

int k26_fs_trash_root(char *out, size_t out_cap)
{
    if (!out || out_cap == 0) return -1;
    const char *home = fsops_home();
    int n = snprintf(out, out_cap, "%s/.local/share/k26/trash", home);
    if (n < 0 || (size_t)n >= out_cap) return -1;
    if (mkdir_p(out) != 0) return -1;
    return 0;
}

int k26_fs_trash_move(const char *src, const char *trash_root,
                      char *out_trashed, size_t cap)
{
    if (!src || !*src || !trash_root || !*trash_root) return -1;

    char root[1024];
    if (trash_root && *trash_root) {
        snprintf(root, sizeof(root), "%s", trash_root);
    } else {
        if (k26_fs_trash_root(root, sizeof(root)) != 0) return -1;
    }
    if (mkdir_p(root) != 0) return -1;

    /* Choose unique destination basename = source basename, suffixed
     * with a counter on collision. */
    const char *slash = strrchr(src, '/');
    const char *base = slash ? slash + 1 : src;
    char dst[2048];
    if (k26_fs_unique_name(root, base, dst, sizeof(dst)) != 0) return -1;

    /* Write sidecar .trashinfo before moving — if the move fails we
     * leave a stale sidecar, but losing the sidecar is cheaper than
     * losing the file. */
    struct stat st;
    if (lstat(src, &st) != 0) return -1;

    char info_path[2048];
    snprintf(info_path, sizeof(info_path), "%s.trashinfo", dst);
    FILE *fp = fopen(info_path, "w");
    if (fp) {
        time_t now = time(NULL);
        struct tm tm;
        localtime_r(&now, &tm);
        char ts[32];
        snprintf(ts, sizeof(ts), "%04d-%02d-%02dT%02d:%02d:%02d",
                 1900 + tm.tm_year, 1 + tm.tm_mon, tm.tm_mday,
                 tm.tm_hour, tm.tm_min, tm.tm_sec);
        fprintf(fp, "[Trash]\n");
        fprintf(fp, "Path=%s\n", src);
        fprintf(fp, "DeletionDate=%s\n", ts);
        fprintf(fp, "OriginalMode=%o\n", (unsigned)st.st_mode);
        fclose(fp);
    }

    /* rename(2) first; on EXDEV a copy + delete would be required.
     * EXDEV is reported as failure here so the caller can fall back
     * to confirm + delete_tree. */
    if (rename(src, dst) != 0) {
        unlink(info_path);
        return -1;
    }
    if (out_trashed && cap) snprintf(out_trashed, cap, "%s", dst);
    return 0;
}

/* ---- walk_size --------------------------------------------------- */

static int walk_size_rec(const char *path, off_t *bytes, int *items,
                         volatile int *cancel)
{
    if (cancel && *cancel) return -2;
    struct stat st;
    if (lstat(path, &st) != 0) return -1;

    if (items) (*items)++;
    if (bytes) (*bytes) += st.st_size;

    if (S_ISDIR(st.st_mode)) {
        DIR *d = opendir(path);
        if (!d) return -1;
        struct dirent *de;
        while ((de = readdir(d)) != NULL) {
            if (cancel && *cancel) { closedir(d); return -2; }
            if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, "..")) continue;
            char child[2048];
            snprintf(child, sizeof(child), "%s/%s", path, de->d_name);
            int r = walk_size_rec(child, bytes, items, cancel);
            if (r == -2) { closedir(d); return -2; }
        }
        closedir(d);
    }
    return 0;
}

int k26_fs_walk_size(const char *root, off_t *bytes, int *items,
                     volatile int *cancel)
{
    if (!root || !*root) return -1;
    if (bytes) *bytes = 0;
    if (items) *items = 0;
    return walk_size_rec(root, bytes, items, cancel);
}

/* ---- human_size -------------------------------------------------- */

int k26_fs_human_size(off_t bytes, char *out, size_t cap)
{
    if (!out || cap == 0) return -1;
    int n;
    if (bytes < 1024)
        n = snprintf(out, cap, "%lld B", (long long)bytes);
    else if (bytes < (off_t)1024 * 1024)
        n = snprintf(out, cap, "%.1f KiB", (double)bytes / 1024.0);
    else if (bytes < (off_t)1024 * 1024 * 1024)
        n = snprintf(out, cap, "%.1f MiB", (double)bytes / (1024.0 * 1024.0));
    else
        n = snprintf(out, cap, "%.2f GiB",
                     (double)bytes / (1024.0 * 1024.0 * 1024.0));
    return (n < 0 || (size_t)n >= cap) ? -1 : n;
}

/* ---- unique_name ------------------------------------------------- */

int k26_fs_unique_name(const char *dir, const char *base,
                       char *out, size_t cap)
{
    if (!dir || !base || !out || cap == 0) return -1;

    char candidate[2048];
    int n = snprintf(candidate, sizeof(candidate), "%s/%s", dir, base);
    if (n < 0 || (size_t)n >= sizeof(candidate)) return -1;

    struct stat st;
    if (lstat(candidate, &st) != 0 && errno == ENOENT) {
        if ((size_t)n >= cap) return -1;
        memcpy(out, candidate, (size_t)n + 1);
        return 0;
    }

    /* Split base into stem + ext to keep the suffix readable.
     * "report.tar.gz" → stem="report", ext=".tar.gz" — but we only
     * separate on the LAST dot to keep the rule simple; if the user
     * has multi-extension files the suffix lands before the final
     * extension, e.g. "report.tar (copy).gz". Single-dot separation
     * is the pragmatic choice for this surface. */
    const char *dot = strrchr(base, '.');
    char stem[512], ext[64];
    if (dot && dot != base) {
        size_t sl = (size_t)(dot - base);
        if (sl >= sizeof(stem)) sl = sizeof(stem) - 1;
        memcpy(stem, base, sl);
        stem[sl] = '\0';
        snprintf(ext, sizeof(ext), "%s", dot);
    } else {
        snprintf(stem, sizeof(stem), "%s", base);
        ext[0] = '\0';
    }

    for (int i = 1; i <= 999; i++) {
        if (i == 1)
            n = snprintf(candidate, sizeof(candidate),
                         "%s/%s (copy)%s", dir, stem, ext);
        else
            n = snprintf(candidate, sizeof(candidate),
                         "%s/%s (copy %d)%s", dir, stem, i, ext);
        if (n < 0 || (size_t)n >= sizeof(candidate)) return -1;
        if (lstat(candidate, &st) != 0 && errno == ENOENT) {
            if ((size_t)n >= cap) return -1;
            memcpy(out, candidate, (size_t)n + 1);
            return 0;
        }
    }
    return -1;
}

/* ---- copy_tree --------------------------------------------------- */

#define K26_FS_COPY_BUF_SZ (64 * 1024)

static void copy_record_failure_(const char *path, char *first_fail, size_t cap)
{
    if (first_fail && cap && first_fail[0] == '\0')
        snprintf(first_fail, cap, "%s", path);
}

static int copy_regular_(const char *src, const char *dst, mode_t mode,
                         time_t mtime,
                         char *first_fail, size_t cap)
{
    int sfd = open(src, O_RDONLY);
    if (sfd < 0) { copy_record_failure_(src, first_fail, cap); return -1; }

    int dfd = open(dst, O_WRONLY | O_CREAT | O_EXCL, mode & 0777);
    if (dfd < 0) {
        copy_record_failure_(dst, first_fail, cap);
        close(sfd);
        return -1;
    }

    static char buf[K26_FS_COPY_BUF_SZ];
    for (;;) {
        ssize_t r = read(sfd, buf, sizeof(buf));
        if (r == 0) break;
        if (r < 0) {
            if (errno == EINTR) continue;
            copy_record_failure_(src, first_fail, cap);
            close(sfd); close(dfd); unlink(dst);
            return -1;
        }
        ssize_t off = 0;
        while (off < r) {
            ssize_t w = write(dfd, buf + off, (size_t)(r - off));
            if (w < 0) {
                if (errno == EINTR) continue;
                copy_record_failure_(dst, first_fail, cap);
                close(sfd); close(dfd); unlink(dst);
                return -1;
            }
            off += w;
        }
    }

    /* Preserve mode (best effort — already created with mode&0777) and
     * mtime via utimensat. fchmod ensures setuid/setgid/sticky come
     * through if the caller had them. */
    fchmod(dfd, mode & 07777);
    close(sfd);
    close(dfd);

    struct timespec ts[2];
    ts[0].tv_sec  = mtime;
    ts[0].tv_nsec = 0;
    ts[1].tv_sec  = mtime;
    ts[1].tv_nsec = 0;
    utimensat(AT_FDCWD, dst, ts, 0);
    return 0;
}

static int copy_tree_rec_(const char *src, const char *dst,
                          k26_progress_cb cb, void *ud,
                          uint64_t *done_bytes, uint64_t total_bytes,
                          char *first_fail, size_t cap)
{
    struct stat st;
    if (lstat(src, &st) != 0) {
        copy_record_failure_(src, first_fail, cap);
        return -1;
    }

    if (cb) {
        if (cb(src, *done_bytes, total_bytes, ud) != 0) {
            errno = ECANCELED;
            return -1;
        }
    }

    if (S_ISLNK(st.st_mode)) {
        char link_target[2048];
        ssize_t ln = readlink(src, link_target, sizeof(link_target) - 1);
        if (ln < 0) { copy_record_failure_(src, first_fail, cap); return -1; }
        link_target[ln] = '\0';
        if (symlink(link_target, dst) != 0) {
            copy_record_failure_(dst, first_fail, cap);
            return -1;
        }
        return 0;
    }

    if (S_ISDIR(st.st_mode)) {
        if (mkdir(dst, st.st_mode & 0777) != 0 && errno != EEXIST) {
            copy_record_failure_(dst, first_fail, cap);
            return -1;
        }
        DIR *d = opendir(src);
        if (!d) { copy_record_failure_(src, first_fail, cap); return -1; }
        struct dirent *de;
        int rc = 0;
        while ((de = readdir(d)) != NULL) {
            if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, "..")) continue;
            char child_src[2048], child_dst[2048];
            snprintf(child_src, sizeof(child_src), "%s/%s", src, de->d_name);
            snprintf(child_dst, sizeof(child_dst), "%s/%s", dst, de->d_name);
            if (copy_tree_rec_(child_src, child_dst, cb, ud,
                               done_bytes, total_bytes,
                               first_fail, cap) != 0) {
                rc = -1;
                break;
            }
        }
        closedir(d);
        /* Preserve directory mtime (set last so children's writes don't
         * bump it). Mode already set via mkdir + umask; fix via chmod. */
        chmod(dst, st.st_mode & 07777);
        struct timespec ts[2];
        ts[0].tv_sec = st.st_mtime; ts[0].tv_nsec = 0;
        ts[1].tv_sec = st.st_mtime; ts[1].tv_nsec = 0;
        utimensat(AT_FDCWD, dst, ts, 0);
        return rc;
    }

    if (S_ISREG(st.st_mode)) {
        int rc = copy_regular_(src, dst, st.st_mode, st.st_mtime,
                               first_fail, cap);
        if (rc == 0) *done_bytes += (uint64_t)st.st_size;
        return rc;
    }

    /* Special files (char/block/fifo/sock): not supported. Skip
     * silently rather than fail — matches cp -r without --preserve. */
    return 0;
}

int k26_fs_copy_tree(const char *src, const char *dst,
                     k26_progress_cb cb, void *ud,
                     char *first_fail, size_t cap)
{
    if (!src || !*src || !dst || !*dst) return -1;
    if (first_fail && cap) first_fail[0] = '\0';

    /* Pre-walk for total bytes (used by progress cb). cancel ptr is
     * NULL — caller cancels via the progress callback returning != 0. */
    off_t total_bytes_off = 0;
    int   total_items     = 0;
    if (k26_fs_walk_size(src, &total_bytes_off, &total_items, NULL) != 0) {
        copy_record_failure_(src, first_fail, cap);
        return -1;
    }
    uint64_t total_bytes = (uint64_t)total_bytes_off;
    uint64_t done_bytes  = 0;

    return copy_tree_rec_(src, dst, cb, ud,
                          &done_bytes, total_bytes,
                          first_fail, cap);
}

/* ---- move -------------------------------------------------------- */

int k26_fs_move(const char *src, const char *dst,
                k26_progress_cb cb, void *ud,
                char *first_fail, size_t cap)
{
    if (!src || !*src || !dst || !*dst) return -1;
    if (first_fail && cap) first_fail[0] = '\0';

    if (rename(src, dst) == 0) return 0;
    if (errno != EXDEV) {
        copy_record_failure_(src, first_fail, cap);
        return -1;
    }

    /* Cross-device: copy the tree, then delete the source. If the
     * delete partially fails we don't unwind the copy — better to have
     * stale source files than lose data. */
    if (k26_fs_copy_tree(src, dst, cb, ud, first_fail, cap) != 0) {
        return -1;
    }
    int failed = 0;
    char delete_fail[1024]; delete_fail[0] = '\0';
    if (k26_fs_delete_tree(src, &failed, delete_fail, sizeof(delete_fail)) != 0) {
        if (first_fail && cap && first_fail[0] == '\0' && delete_fail[0])
            snprintf(first_fail, cap, "%s", delete_fail);
        return -1;
    }
    return 0;
}
