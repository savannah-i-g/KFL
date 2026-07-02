/* libk26util: string helpers exposed to KFL's expression layer.
 * See k26_str.h for the contract. */

#include "k26_str.h"

#include <string.h>

size_t k26_str_len(const char *s)
{
    return s ? strlen(s) : 0u;
}

int k26_str_eq(const char *a, const char *b)
{
    if (!a || !b) return 0;
    return strcmp(a, b) == 0;
}

int k26_str_starts_with(const char *s, const char *prefix)
{
    if (!s || !prefix) return 0;
    size_t pl = strlen(prefix);
    return strncmp(s, prefix, pl) == 0;
}

int k26_str_ends_with(const char *s, const char *suffix)
{
    if (!s || !suffix) return 0;
    size_t sl = strlen(s);
    size_t tl = strlen(suffix);
    if (tl > sl) return 0;
    return strcmp(s + (sl - tl), suffix) == 0;
}

size_t k26_str_concat(char *out, size_t cap, const char *a, const char *b)
{
    if (!out || cap == 0) return 0;
    out[0] = '\0';
    size_t off = 0;
    if (a) {
        size_t la = strlen(a);
        if (la > cap - 1) la = cap - 1;
        memcpy(out, a, la);
        off = la;
        out[off] = '\0';
    }
    if (b && off + 1 < cap) {
        size_t lb = strlen(b);
        size_t room = cap - off - 1;
        if (lb > room) lb = room;
        memcpy(out + off, b, lb);
        off += lb;
        out[off] = '\0';
    }
    return off;
}
