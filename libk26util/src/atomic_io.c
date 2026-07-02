#include "k26_atomic_io.h"

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

int k26_atomic_write(const char *path, const void *data, size_t len, mode_t mode)
{
    if (!path || (!data && len)) { errno = EINVAL; return -1; }

    char tmp[512];
    if ((size_t)snprintf(tmp, sizeof tmp, "%s.new", path) >= sizeof tmp) {
        errno = ENAMETOOLONG;
        return -1;
    }

    int fd = open(tmp, O_WRONLY | O_CREAT | O_TRUNC, mode);
    if (fd < 0) return -1;

    const char *p = (const char *)data;
    size_t left = len;
    while (left > 0) {
        ssize_t w = write(fd, p, left);
        if (w < 0) {
            if (errno == EINTR) continue;
            close(fd);
            unlink(tmp);
            return -1;
        }
        p    += w;
        left -= (size_t)w;
    }

    if (fsync(fd) != 0) {
        int saved = errno;
        close(fd);
        unlink(tmp);
        errno = saved;
        return -1;
    }
    close(fd);

    if (rename(tmp, path) != 0) {
        int saved = errno;
        unlink(tmp);
        errno = saved;
        return -1;
    }
    return 0;
}

int k26_atomic_read(const char *path, char **out, size_t *outlen)
{
    if (!path || !out) { errno = EINVAL; return -1; }
    *out = NULL;
    if (outlen) *outlen = 0;

    int fd = open(path, O_RDONLY);
    if (fd < 0) return -1;

    struct stat st;
    if (fstat(fd, &st) != 0) { int saved = errno; close(fd); errno = saved; return -1; }
    if (!S_ISREG(st.st_mode)) { close(fd); errno = EINVAL; return -1; }

    size_t cap = (size_t)st.st_size;
    char *buf = malloc(cap + 1);
    if (!buf) { close(fd); errno = ENOMEM; return -1; }

    size_t got = 0;
    while (got < cap) {
        ssize_t r = read(fd, buf + got, cap - got);
        if (r < 0) {
            if (errno == EINTR) continue;
            int saved = errno;
            free(buf);
            close(fd);
            errno = saved;
            return -1;
        }
        if (r == 0) break;            /* file shrunk under us */
        got += (size_t)r;
    }
    close(fd);

    buf[got] = '\0';
    *out = buf;
    if (outlen) *outlen = got;
    return 0;
}
