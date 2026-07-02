/* kflc — arena allocator + diagnostics. */

#include "kflc.h"
#include "internal.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

typedef struct ArenaChunk {
    struct ArenaChunk *next;
    size_t             cap;
    size_t             used;
    char               data[];
} ArenaChunk;

struct KflcArena {
    ArenaChunk *head;
};

#define ARENA_DEFAULT_CHUNK (64u * 1024u)
#define ARENA_ALIGN          16u

static size_t align_up(size_t n, size_t a)
{
    return (n + (a - 1)) & ~(a - 1);
}

KflcArena *kflc_arena_create(void)
{
    KflcArena *a = (KflcArena *)calloc(1, sizeof(*a));
    return a;
}

void kflc_arena_release(KflcArena *a)
{
    if (!a) return;
    ArenaChunk *c = a->head;
    while (c) {
        ArenaChunk *next = c->next;
        free(c);
        c = next;
    }
    free(a);
}

void *kflc_arena_alloc(KflcArena *a, size_t n)
{
    if (!a || n == 0) return NULL;
    size_t need = align_up(n, ARENA_ALIGN);
    ArenaChunk *c = a->head;
    if (!c || c->used + need > c->cap) {
        size_t cap = ARENA_DEFAULT_CHUNK;
        if (cap < need) cap = need;
        ArenaChunk *fresh = (ArenaChunk *)malloc(sizeof(*fresh) + cap);
        if (!fresh) return NULL;
        fresh->next = a->head;
        fresh->cap  = cap;
        fresh->used = 0;
        a->head = fresh;
        c = fresh;
    }
    void *p = c->data + c->used;
    c->used += need;
    memset(p, 0, n);
    return p;
}

char *kflc_arena_strdup(KflcArena *a, const char *s)
{
    if (!s) return NULL;
    return kflc_arena_strndup(a, s, strlen(s));
}

char *kflc_arena_strndup(KflcArena *a, const char *s, size_t n)
{
    if (!s) return NULL;
    char *out = (char *)kflc_arena_alloc(a, n + 1);
    if (!out) return NULL;
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

/* ---- Diagnostics -------------------------------------------------- */

void kflc_diag_init(KflcDiag *d, const char *path, FILE *errstream)
{
    d->path      = path ? path : "<input>";
    d->errors    = 0;
    d->warnings  = 0;
    d->errstream = errstream ? errstream : stderr;
}

/* Shared printer. `col <= 0` suppresses the column field so editors
 * that only understand `<path>:<line>:` still parse the location. */
static void diag_print_v_(KflcDiag *d, const char *severity,
                          int line, int col,
                          const char *fmt, va_list ap)
{
    if (col > 0) {
        fprintf(d->errstream, "%s:%d:%d: %s: ",
                d->path, line, col, severity);
    } else {
        fprintf(d->errstream, "%s:%d: %s: ",
                d->path, line, severity);
    }
    vfprintf(d->errstream, fmt, ap);
    fputc('\n', d->errstream);
}

void kflc_diag_errorf(KflcDiag *d, int line, const char *fmt, ...)
{
    if (!d) return;
    d->errors++;
    va_list ap;
    va_start(ap, fmt);
    diag_print_v_(d, "error", line, 0, fmt, ap);
    va_end(ap);
}

void kflc_diag_errorf_col(KflcDiag *d, int line, int col,
                          const char *fmt, ...)
{
    if (!d) return;
    d->errors++;
    va_list ap;
    va_start(ap, fmt);
    diag_print_v_(d, "error", line, col, fmt, ap);
    va_end(ap);
}

void kflc_diag_warnf(KflcDiag *d, int line, const char *fmt, ...)
{
    if (!d) return;
    d->warnings++;
    va_list ap;
    va_start(ap, fmt);
    diag_print_v_(d, "warning", line, 0, fmt, ap);
    va_end(ap);
}

void kflc_diag_warnf_col(KflcDiag *d, int line, int col,
                         const char *fmt, ...)
{
    if (!d) return;
    d->warnings++;
    va_list ap;
    va_start(ap, fmt);
    diag_print_v_(d, "warning", line, col, fmt, ap);
    va_end(ap);
}

/* ---- Reserved-word table ----------------------------------------- */

int kfl_is_reserved_future(const char *name)
{
    if (!name) return 0;
    static const char *const RESERVED[] = {
        "thread", "join", "mutex", "lock", "unlock",
        "parallel_for", "bind", "import",
        NULL
    };
    for (int i = 0; RESERVED[i]; i++) {
        if (strcmp(name, RESERVED[i]) == 0) return 1;
    }
    return 0;
}

/* ---- AST structural equality ------------------------------------- */

static int str_equal_(const char *a, const char *b)
{
    if (a == b) return 1;
    if (!a || !b) return 0;
    return strcmp(a, b) == 0;
}

static int value_equal_(const KflcValue *a, const KflcValue *b)
{
    if (a->kind != b->kind) return 0;
    switch (a->kind) {
    case KFLV_NONE:
        return 1;
    case KFLV_INT:
        return a->u.i == b->u.i;
    case KFLV_FLOAT:
        /* Bitwise compare. Documented behaviour: if a round-trip
         * surfaces a mismatch here, fix the serializer precision
         * (use %.17g) rather than weakening this check. */
        return a->u.f == b->u.f;
    case KFLV_STR:
    case KFLV_IDENT:
    case KFLV_SHORTCUT:
        return str_equal_(a->u.s, b->u.s);
    case KFLV_COLOR:
        return a->u.rgb == b->u.rgb;
    case KFLV_SIZE:
        return a->u.size.w == b->u.size.w &&
               a->u.size.h == b->u.size.h;
    case KFLV_ALIGN:
        return a->u.align == b->u.align;
    }
    return 0;
}

static int expr_equal_(const KflcExpr *a, const KflcExpr *b)
{
    if (a == b) return 1;
    if (!a || !b) return 0;
    if (a->kind != b->kind) return 0;
    switch (a->kind) {
    case KFLE_INT_LIT:
        return a->u.i == b->u.i;
    case KFLE_FLOAT_LIT:
        return a->u.f == b->u.f;
    case KFLE_IDENT:
        return str_equal_(a->u.ident, b->u.ident);
    case KFLE_UNARY:
        return a->u.un.op == b->u.un.op &&
               expr_equal_(a->u.un.operand, b->u.un.operand);
    case KFLE_BINARY:
        return a->u.bin.op == b->u.bin.op &&
               expr_equal_(a->u.bin.lhs, b->u.bin.lhs) &&
               expr_equal_(a->u.bin.rhs, b->u.bin.rhs);
    case KFLE_CALL: {
        if (!str_equal_(a->u.call.name, b->u.call.name)) return 0;
        if (a->u.call.n_args != b->u.call.n_args) return 0;
        for (int i = 0; i < a->u.call.n_args; i++) {
            if (!expr_equal_(a->u.call.args[i],
                             b->u.call.args[i])) return 0;
        }
        return 1;
    }
    case KFLE_VEC_LIT: {
        if (a->u.vec.n_elems != b->u.vec.n_elems) return 0;
        for (int i = 0; i < a->u.vec.n_elems; i++) {
            if (!expr_equal_(a->u.vec.elems[i],
                             b->u.vec.elems[i])) return 0;
        }
        return 1;
    }
    case KFLE_INDEX:
        return expr_equal_(a->u.index.base, b->u.index.base) &&
               expr_equal_(a->u.index.idx,  b->u.index.idx);
    }
    return 0;
}

static int attr_list_equal_(const KflcAttr *a, const KflcAttr *b)
{
    while (a && b) {
        if (!str_equal_(a->name, b->name)) return 0;
        if (!value_equal_(&a->value, &b->value)) return 0;
        if (!expr_equal_(a->expr, b->expr)) return 0;
        a = a->next;
        b = b->next;
    }
    return a == NULL && b == NULL;
}

static int node_equal_(const KflcNode *a, const KflcNode *b);

static int node_chain_equal_(const KflcNode *a, const KflcNode *b)
{
    while (a && b) {
        if (!node_equal_(a, b)) return 0;
        a = a->next;
        b = b->next;
    }
    return a == NULL && b == NULL;
}

static int node_equal_(const KflcNode *a, const KflcNode *b)
{
    if (a == b) return 1;
    if (!a || !b) return 0;
    if (a->kind      != b->kind)      return 0;
    if (a->container != b->container) return 0;
    if (a->widget    != b->widget)    return 0;
    if (a->type      != b->type)      return 0;
    if (!str_equal_(a->type_subtype, b->type_subtype)) return 0;
    if (a->flags     != b->flags)     return 0;
    if (!str_equal_(a->name, b->name)) return 0;
    if (!value_equal_(&a->position, &b->position)) return 0;
    if (!attr_list_equal_(a->attrs, b->attrs)) return 0;
    if (!node_chain_equal_(a->children,      b->children))      return 0;
    if (!node_chain_equal_(a->else_children, b->else_children)) return 0;
    if (!expr_equal_(a->expr,  b->expr))  return 0;
    if (!expr_equal_(a->expr2, b->expr2)) return 0;
    /* line numbers intentionally skipped — serializer may renumber */
    return 1;
}

int kflc_ast_equal(const KflcNode *a, const KflcNode *b)
{
    return node_chain_equal_(a, b);
}
