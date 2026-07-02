/* runtime.c: KFL memory runtime.
 *
 * See include/kflc_runtime.h for the API surface and rationale.
 * This file is small by design: ~150 LOC of bump allocator + move
 * primitive. The arena uses a single malloc-backed buffer with an
 * 8-byte stride alignment (alignof(max_align_t) on every supported
 * target). Allocation is bump-and-check; out-of-capacity returns
 * NULL and the caller's emit-time check for NULL is the only error
 * surface.
 */

#include "kflc_runtime.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

struct K26KflArena {
    unsigned char *buf;
    size_t         capacity;
    size_t         used;
};

#define K26KFL_ALIGN  (sizeof(long double))  /* >= max_align_t on x86_64 */

static size_t align_up_(size_t n)
{
    return (n + (K26KFL_ALIGN - 1)) & ~(K26KFL_ALIGN - 1);
}

K26KflArena *k26kfl_arena_create(size_t capacity_bytes)
{
    if (capacity_bytes == 0) return NULL;
    K26KflArena *a = (K26KflArena *)malloc(sizeof(*a));
    if (!a) return NULL;
    a->buf = (unsigned char *)malloc(capacity_bytes);
    if (!a->buf) { free(a); return NULL; }
    a->capacity = capacity_bytes;
    a->used     = 0;
    return a;
}

void k26kfl_arena_destroy(K26KflArena *a)
{
    if (!a) return;
    free(a->buf);
    free(a);
}

void *k26kfl_arena_alloc(K26KflArena *a, size_t bytes)
{
    if (!a || bytes == 0) return NULL;
    size_t need = align_up_(bytes);
    if (a->used + need > a->capacity) return NULL;
    void *p = a->buf + a->used;
    a->used += need;
    return p;
}

void k26kfl_arena_reset(K26KflArena *a)
{
    if (!a) return;
    a->used = 0;
}

size_t k26kfl_arena_used(const K26KflArena *a)
{
    return a ? a->used : 0;
}

size_t k26kfl_arena_capacity(const K26KflArena *a)
{
    return a ? a->capacity : 0;
}

void k26kfl_move(void **dst, void **src)
{
    if (!dst || !src) return;
    *dst = *src;
    *src = NULL;
}
