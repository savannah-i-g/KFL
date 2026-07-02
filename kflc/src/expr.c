/* kflc — expression sub-language.
 *
 * A self-contained lexer + parser + emitter for KFL expressions.
 * Operates on a NUL-terminated source string (typically the contents
 * of an `expression "..."` attribute, or the right-hand side of an
 * inline `fn` statement).
 *
 * Grammar (precedence climbing, low to high):
 *   logical_or  = logical_and  ( "||" logical_and )*
 *   logical_and = comparison   ( "&&" comparison )*
 *   comparison  = additive     ( ("<"|"<="|">"|">="|"=="|"!=") additive )*
 *   additive    = multiplicative ( ("+"|"-") multiplicative )*
 *   multiplicative = unary     ( ("*"|"/"|"%") unary )*
 *   unary       = ("-"|"+"|"!") unary | primary
 *   primary     = NUMBER | IDENT [ "(" arglist ")" ] | "(" logical_or ")"
 *   arglist     = ( logical_or ( "," logical_or )* )?
 *
 * Literals, unary +/-, the five arithmetic operators, parenthesised
 * sub-expressions, and function calls against a libm / libk26compute
 * whitelist resolved by the emitter are supported. The lexer accepts
 * comparison + logical operators; the emitter only emits them in
 * predicate contexts (`if` / `while`).
 */

#include "kflc.h"
#include "internal.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Whitelist of identifiers callable from an expression -------- */

/* Each entry maps a KFL function name to a C++ expression. The emitter
 * looks up the KFL name and emits the C++ form; this lets us rename
 * (e.g. KFL `abs` → C++ `fabs`) and gives a single point of truth for
 * the set of available calls. `arity` is checked at emit time; -1
 * means variadic and is currently unused. */
typedef struct {
    const char *kfl_name;
    const char *cxx_name;
    int         arity;
} BuiltinFn;

static const BuiltinFn BUILTINS[] = {
    /* libm scalar — one argument unless noted */
    { "sin",   "sin",     1 },
    { "cos",   "cos",     1 },
    { "tan",   "tan",     1 },
    { "asin",  "asin",    1 },
    { "acos",  "acos",    1 },
    { "atan",  "atan",    1 },
    { "atan2", "atan2",   2 },
    { "exp",   "exp",     1 },
    { "log",   "log",     1 },
    { "log10", "log10",   1 },
    { "sqrt",  "sqrt",    1 },
    { "pow",   "pow",     2 },
    { "abs",   "fabs",    1 },
    { "floor", "floor",   1 },
    { "ceil",  "ceil",    1 },
    { "round", "round",   1 },
    { "min",   "fmin",    2 },
    { "max",   "fmax",    2 },
    { "fmod",  "fmod",    2 },
    /* String helpers from libk26util. */
    { "strlen",          "k26_str_len",                    1 },
    { "streq",           "k26_str_eq",                     2 },
    { "starts_with",     "k26_str_starts_with",            2 },
    { "ends_with",       "k26_str_ends_with",              2 },
    { "concat",          "k26_str_concat",                 2 },
    { NULL,    NULL,      0 }
};

/* Dynamic registry, populated by kflc_register_builtin. External
 * libraries (or the future .kflbi manifest loader) call the public API
 * before parsing starts; lookup_builtin walks the static table first,
 * then this one. Cap is generous — astro alone projects ~100 entries
 * across its 6+1 libs. The strings are caller-owned; the loader passes
 * arena-allocated copies. */
#define KFLC_DYNAMIC_BUILTIN_MAX 256
static BuiltinFn g_dynamic_builtins[KFLC_DYNAMIC_BUILTIN_MAX];
static int       g_dynamic_builtin_count = 0;

int kflc_register_builtin(const char *kfl_name, const char *cxx_name, int arity)
{
    if (!kfl_name || !cxx_name) return 1;
    /* Refuse a duplicate KFL name regardless of which table it's in —
     * collisions usually mean a manifest conflict that the caller
     * should resolve, not a behaviour we should paper over. */
    for (const BuiltinFn *b = BUILTINS; b->kfl_name; b++) {
        if (strcmp(b->kfl_name, kfl_name) == 0) return 2;
    }
    for (int i = 0; i < g_dynamic_builtin_count; i++) {
        if (strcmp(g_dynamic_builtins[i].kfl_name, kfl_name) == 0) return 2;
    }
    if (g_dynamic_builtin_count >= KFLC_DYNAMIC_BUILTIN_MAX) return 3;
    g_dynamic_builtins[g_dynamic_builtin_count].kfl_name = kfl_name;
    g_dynamic_builtins[g_dynamic_builtin_count].cxx_name = cxx_name;
    g_dynamic_builtins[g_dynamic_builtin_count].arity    = arity;
    g_dynamic_builtin_count++;
    return 0;
}

void kflc_clear_builtins(void)
{
    g_dynamic_builtin_count = 0;
}

/* Returns 1 if `name` is one of the string-flavoured builtins that
 * needs custom emit (return-type ≠ double, or scratch-buffer trick).
 * Centralised here so the KFLE_CALL switch can fast-path them. */
static int is_string_builtin_(const char *name)
{
    if (!name) return 0;
    return (strcmp(name, "strlen")       == 0 ||
            strcmp(name, "streq")        == 0 ||
            strcmp(name, "starts_with")  == 0 ||
            strcmp(name, "ends_with")    == 0 ||
            strcmp(name, "concat")       == 0);
}

static const BuiltinFn *lookup_builtin(const char *name)
{
    if (!name) return NULL;
    for (const BuiltinFn *b = BUILTINS; b->kfl_name; b++) {
        if (strcmp(b->kfl_name, name) == 0) return b;
    }
    for (int i = 0; i < g_dynamic_builtin_count; i++) {
        if (strcmp(g_dynamic_builtins[i].kfl_name, name) == 0)
            return &g_dynamic_builtins[i];
    }
    return NULL;
}

/* ---- Expression lexer -------------------------------------------- */

typedef enum {
    ET_EOF = 0,
    ET_NUMBER,
    ET_IDENT,
    ET_LPAREN,
    ET_RPAREN,
    ET_LBRACK,
    ET_RBRACK,
    ET_COMMA,
    ET_PLUS,
    ET_MINUS,
    ET_STAR,
    ET_SLASH,
    ET_PERCENT,
    ET_BANG,
    ET_LT,
    ET_LE,
    ET_GT,
    ET_GE,
    ET_EQ,
    ET_NE,
    ET_ANDAND,
    ET_OROR,
    ET_ERROR
} ETKind;

typedef struct {
    ETKind  kind;
    /* For NUMBER: parsed value + is_float flag. */
    double  fval;
    long    ival;
    int     is_float;
    /* For IDENT: arena-owned NUL-terminated string. */
    char   *str;
} EToken;

typedef struct {
    const char *src;
    size_t      pos;
    size_t      len;
    int         line;
    KflcArena  *arena;
    KflcDiag   *diag;
    int         had_error;
} EL;

static void el_skip_ws(EL *L)
{
    while (L->pos < L->len) {
        int c = (unsigned char)L->src[L->pos];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            L->pos++;
            continue;
        }
        break;
    }
}

static int el_is_ident_start(int c) { return c == '_' || isalpha(c); }
static int el_is_ident_cont (int c) { return c == '_' || isalnum(c); }

static EToken el_next(EL *L)
{
    EToken t;
    memset(&t, 0, sizeof(t));
    el_skip_ws(L);
    if (L->pos >= L->len) { t.kind = ET_EOF; return t; }
    int c = (unsigned char)L->src[L->pos];

    if (isdigit(c) || (c == '.' && L->pos + 1 < L->len && isdigit((unsigned char)L->src[L->pos + 1]))) {
        /* Number: a run of digits with optional fractional part. */
        size_t start = L->pos;
        int is_float = 0;
        while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
        if (L->pos < L->len && L->src[L->pos] == '.') {
            is_float = 1;
            L->pos++;
            while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
        }
        /* Optional exponent. */
        if (L->pos < L->len && (L->src[L->pos] == 'e' || L->src[L->pos] == 'E')) {
            is_float = 1;
            L->pos++;
            if (L->pos < L->len && (L->src[L->pos] == '+' || L->src[L->pos] == '-')) L->pos++;
            while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
        }
        char buf[64];
        size_t n = L->pos - start;
        if (n >= sizeof(buf)) {
            kflc_diag_errorf(L->diag, L->line, "expr: number literal too long");
            L->had_error = 1;
            t.kind = ET_ERROR;
            return t;
        }
        memcpy(buf, L->src + start, n);
        buf[n] = '\0';
        t.kind = ET_NUMBER;
        t.is_float = is_float;
        if (is_float) t.fval = strtod(buf, NULL);
        else          t.ival = strtol(buf, NULL, 10);
        return t;
    }

    if (el_is_ident_start(c)) {
        size_t start = L->pos;
        while (L->pos < L->len && el_is_ident_cont((unsigned char)L->src[L->pos])) L->pos++;
        size_t n = L->pos - start;
        t.kind = ET_IDENT;
        t.str  = kflc_arena_strndup(L->arena, L->src + start, n);
        return t;
    }

    L->pos++;
    switch (c) {
    case '(': t.kind = ET_LPAREN;  return t;
    case ')': t.kind = ET_RPAREN;  return t;
    case '[': t.kind = ET_LBRACK;  return t;
    case ']': t.kind = ET_RBRACK;  return t;
    case ',': t.kind = ET_COMMA;   return t;
    case '+': t.kind = ET_PLUS;    return t;
    case '-': t.kind = ET_MINUS;   return t;
    case '*': t.kind = ET_STAR;    return t;
    case '/': t.kind = ET_SLASH;   return t;
    case '%': t.kind = ET_PERCENT; return t;
    case '<':
        if (L->pos < L->len && L->src[L->pos] == '=') { L->pos++; t.kind = ET_LE; return t; }
        t.kind = ET_LT; return t;
    case '>':
        if (L->pos < L->len && L->src[L->pos] == '=') { L->pos++; t.kind = ET_GE; return t; }
        t.kind = ET_GT; return t;
    case '=':
        if (L->pos < L->len && L->src[L->pos] == '=') { L->pos++; t.kind = ET_EQ; return t; }
        kflc_diag_errorf(L->diag, L->line, "expr: bare `=` is not yet an operator (use `==`)");
        L->had_error = 1;
        t.kind = ET_ERROR; return t;
    case '!':
        if (L->pos < L->len && L->src[L->pos] == '=') { L->pos++; t.kind = ET_NE; return t; }
        t.kind = ET_BANG; return t;
    case '&':
        if (L->pos < L->len && L->src[L->pos] == '&') { L->pos++; t.kind = ET_ANDAND; return t; }
        break;
    case '|':
        if (L->pos < L->len && L->src[L->pos] == '|') { L->pos++; t.kind = ET_OROR; return t; }
        break;
    }
    kflc_diag_errorf(L->diag, L->line,
                     "expr: unexpected character `%c` (0x%02x)",
                     (c >= 32 && c < 127) ? c : '?', c);
    L->had_error = 1;
    t.kind = ET_ERROR;
    return t;
}

/* ---- Expression parser ------------------------------------------- */

typedef struct {
    EL      L;
    EToken  cur;
    KflcArena *arena;
    KflcDiag  *diag;
    int        had_error;
    int        line;
} EP;

static void ep_advance(EP *P)
{
    P->cur = el_next(&P->L);
    if (P->cur.kind == ET_ERROR) P->had_error = 1;
}

static KflcExpr *new_expr(EP *P, KflcExprKind k)
{
    KflcExpr *e = (KflcExpr *)kflc_arena_alloc(P->arena, sizeof(*e));
    e->kind = k;
    e->line = P->line;
    return e;
}

/* Forward decl: top-level parse. */
static KflcExpr *parse_or(EP *P);

static KflcExpr *parse_primary(EP *P)
{
    if (P->cur.kind == ET_NUMBER) {
        KflcExpr *e;
        if (P->cur.is_float) {
            e = new_expr(P, KFLE_FLOAT_LIT);
            e->u.f = P->cur.fval;
        } else {
            e = new_expr(P, KFLE_INT_LIT);
            e->u.i = P->cur.ival;
        }
        ep_advance(P);
        return e;
    }
    if (P->cur.kind == ET_IDENT) {
        char *name = P->cur.str;
        ep_advance(P);
        /* Postfix `[...]` chain. Single-level (`v[i]`) is the vector
         * form; two-level (`m[i][j]`) is matrix cell access. Build
         * left-associatively so the outermost INDEX node holds the
         * rightmost subscript; the emitter walks the chain inside-out
         * for the deref. The binding-type check (vector vs matrix vs
         * scalar) happens at emit time so a single error message can
         * cite the actual ident + type. */
        if (P->cur.kind == ET_LBRACK) {
            KflcExpr *base = new_expr(P, KFLE_IDENT);
            base->u.ident = name;
            while (P->cur.kind == ET_LBRACK) {
                ep_advance(P);
                KflcExpr *idx = parse_or(P);
                if (!idx) return NULL;
                if (P->cur.kind != ET_RBRACK) {
                    kflc_diag_errorf(P->diag, P->line,
                                     "expr: expected `]` after index");
                    P->had_error = 1;
                    return NULL;
                }
                ep_advance(P);
                KflcExpr *e = new_expr(P, KFLE_INDEX);
                e->u.index.base = base;
                e->u.index.idx  = idx;
                base = e;
            }
            return base;
        }
        if (P->cur.kind == ET_LPAREN) {
            ep_advance(P);
            KflcExpr *call = new_expr(P, KFLE_CALL);
            call->u.call.name   = name;
            call->u.call.args   = NULL;
            call->u.call.n_args = 0;
            int  cap = 0;
            if (P->cur.kind != ET_RPAREN) {
                for (;;) {
                    KflcExpr *a = parse_or(P);
                    if (!a) return NULL;
                    if (call->u.call.n_args == cap) {
                        int nc = cap ? cap * 2 : 4;
                        KflcExpr **na = (KflcExpr **)kflc_arena_alloc(
                            P->arena, sizeof(KflcExpr *) * (size_t)nc);
                        if (call->u.call.n_args > 0) {
                            memcpy(na, call->u.call.args,
                                   sizeof(KflcExpr *) * (size_t)call->u.call.n_args);
                        }
                        call->u.call.args = na;
                        cap = nc;
                    }
                    call->u.call.args[call->u.call.n_args++] = a;
                    if (P->cur.kind == ET_COMMA) { ep_advance(P); continue; }
                    break;
                }
            }
            if (P->cur.kind != ET_RPAREN) {
                kflc_diag_errorf(P->diag, P->line,
                                 "expr: expected `)` in call to `%s`", name);
                P->had_error = 1;
                return NULL;
            }
            ep_advance(P);
            return call;
        }
        KflcExpr *e = new_expr(P, KFLE_IDENT);
        e->u.ident = name;
        return e;
    }
    if (P->cur.kind == ET_LPAREN) {
        ep_advance(P);
        KflcExpr *e = parse_or(P);
        if (!e) return NULL;
        if (P->cur.kind != ET_RPAREN) {
            kflc_diag_errorf(P->diag, P->line, "expr: expected `)`");
            P->had_error = 1;
            return NULL;
        }
        ep_advance(P);
        return e;
    }
    if (P->cur.kind == ET_LBRACK) {
        /* Vector or matrix literal. Vector elements are scalar
         * expressions; matrix elements are themselves vector literals
         * (rows). The parser builds a flat KFLE_VEC_LIT either way; the
         * downstream emitter checks if every element is itself a
         * KFLE_VEC_LIT and switches to matrix codegen if so. */
        ep_advance(P);
        KflcExpr *v = new_expr(P, KFLE_VEC_LIT);
        v->u.vec.elems   = NULL;
        v->u.vec.n_elems = 0;
        int cap = 0;
        if (P->cur.kind != ET_RBRACK) {
            for (;;) {
                KflcExpr *e = parse_or(P);
                if (!e) return NULL;
                if (v->u.vec.n_elems == cap) {
                    int nc = cap ? cap * 2 : 4;
                    KflcExpr **na = (KflcExpr **)kflc_arena_alloc(
                        P->arena, sizeof(KflcExpr *) * (size_t)nc);
                    if (v->u.vec.n_elems > 0) {
                        memcpy(na, v->u.vec.elems,
                               sizeof(KflcExpr *) * (size_t)v->u.vec.n_elems);
                    }
                    v->u.vec.elems = na;
                    cap = nc;
                }
                v->u.vec.elems[v->u.vec.n_elems++] = e;
                if (P->cur.kind == ET_COMMA) { ep_advance(P); continue; }
                break;
            }
        }
        if (P->cur.kind != ET_RBRACK) {
            kflc_diag_errorf(P->diag, P->line,
                             "expr: expected `]` to close vector/matrix literal");
            P->had_error = 1;
            return NULL;
        }
        ep_advance(P);
        return v;
    }
    kflc_diag_errorf(P->diag, P->line,
                     "expr: unexpected token (kind %d) in primary expression",
                     P->cur.kind);
    P->had_error = 1;
    return NULL;
}

static KflcExpr *parse_unary(EP *P)
{
    if (P->cur.kind == ET_MINUS || P->cur.kind == ET_PLUS || P->cur.kind == ET_BANG) {
        ETKind k = P->cur.kind;
        ep_advance(P);
        KflcExpr *operand = parse_unary(P);
        if (!operand) return NULL;
        KflcExpr *u = new_expr(P, KFLE_UNARY);
        u->u.un.op = (k == ET_MINUS) ? KFLOP_NEG
                  : (k == ET_PLUS)  ? KFLOP_POS
                  :                   KFLOP_NOT;
        u->u.un.operand = operand;
        return u;
    }
    return parse_primary(P);
}

static KflcExpr *parse_mul(EP *P)
{
    KflcExpr *lhs = parse_unary(P);
    if (!lhs) return NULL;
    while (P->cur.kind == ET_STAR || P->cur.kind == ET_SLASH || P->cur.kind == ET_PERCENT) {
        ETKind k = P->cur.kind;
        ep_advance(P);
        KflcExpr *rhs = parse_unary(P);
        if (!rhs) return NULL;
        KflcExpr *b = new_expr(P, KFLE_BINARY);
        b->u.bin.op = (k == ET_STAR)  ? KFLOP_MUL
                   : (k == ET_SLASH) ? KFLOP_DIV
                   :                   KFLOP_MOD;
        b->u.bin.lhs = lhs;
        b->u.bin.rhs = rhs;
        lhs = b;
    }
    return lhs;
}

static KflcExpr *parse_add(EP *P)
{
    KflcExpr *lhs = parse_mul(P);
    if (!lhs) return NULL;
    while (P->cur.kind == ET_PLUS || P->cur.kind == ET_MINUS) {
        ETKind k = P->cur.kind;
        ep_advance(P);
        KflcExpr *rhs = parse_mul(P);
        if (!rhs) return NULL;
        KflcExpr *b = new_expr(P, KFLE_BINARY);
        b->u.bin.op = (k == ET_PLUS) ? KFLOP_ADD : KFLOP_SUB;
        b->u.bin.lhs = lhs;
        b->u.bin.rhs = rhs;
        lhs = b;
    }
    return lhs;
}

static KflcExpr *parse_cmp(EP *P)
{
    KflcExpr *lhs = parse_add(P);
    if (!lhs) return NULL;
    while (P->cur.kind == ET_LT || P->cur.kind == ET_LE ||
           P->cur.kind == ET_GT || P->cur.kind == ET_GE ||
           P->cur.kind == ET_EQ || P->cur.kind == ET_NE)
    {
        ETKind k = P->cur.kind;
        ep_advance(P);
        KflcExpr *rhs = parse_add(P);
        if (!rhs) return NULL;
        KflcExpr *b = new_expr(P, KFLE_BINARY);
        switch (k) {
        case ET_LT: b->u.bin.op = KFLOP_LT; break;
        case ET_LE: b->u.bin.op = KFLOP_LE; break;
        case ET_GT: b->u.bin.op = KFLOP_GT; break;
        case ET_GE: b->u.bin.op = KFLOP_GE; break;
        case ET_EQ: b->u.bin.op = KFLOP_EQ; break;
        case ET_NE: b->u.bin.op = KFLOP_NE; break;
        default: break;
        }
        b->u.bin.lhs = lhs;
        b->u.bin.rhs = rhs;
        lhs = b;
    }
    return lhs;
}

static KflcExpr *parse_and(EP *P)
{
    KflcExpr *lhs = parse_cmp(P);
    if (!lhs) return NULL;
    while (P->cur.kind == ET_ANDAND) {
        ep_advance(P);
        KflcExpr *rhs = parse_cmp(P);
        if (!rhs) return NULL;
        KflcExpr *b = new_expr(P, KFLE_BINARY);
        b->u.bin.op = KFLOP_AND;
        b->u.bin.lhs = lhs;
        b->u.bin.rhs = rhs;
        lhs = b;
    }
    return lhs;
}

static KflcExpr *parse_or(EP *P)
{
    KflcExpr *lhs = parse_and(P);
    if (!lhs) return NULL;
    while (P->cur.kind == ET_OROR) {
        ep_advance(P);
        KflcExpr *rhs = parse_and(P);
        if (!rhs) return NULL;
        KflcExpr *b = new_expr(P, KFLE_BINARY);
        b->u.bin.op = KFLOP_OR;
        b->u.bin.lhs = lhs;
        b->u.bin.rhs = rhs;
        lhs = b;
    }
    return lhs;
}

KflcExpr *kflc_parse_expr(const char *src, KflcArena *arena,
                          KflcDiag *diag, int line_base)
{
    if (!src) return NULL;
    EP P;
    memset(&P, 0, sizeof(P));
    P.arena     = arena;
    P.diag      = diag;
    P.line      = line_base > 0 ? line_base : 1;
    P.L.src     = src;
    P.L.pos     = 0;
    P.L.len     = strlen(src);
    P.L.line    = P.line;
    P.L.arena   = arena;
    P.L.diag    = diag;
    P.L.had_error = 0;
    ep_advance(&P);
    KflcExpr *e = parse_or(&P);
    if (P.had_error || P.L.had_error) return NULL;
    if (P.cur.kind != ET_EOF) {
        kflc_diag_errorf(diag, P.line,
                         "expr: trailing tokens after expression");
        return NULL;
    }
    return e;
}

/* ---- Expression emitter ----------------------------------------- */

static const KflcExprBinding *lookup_binding(const char *name,
                                              const KflcExprCtx *ctx)
{
    if (!ctx || !ctx->bindings) return NULL;
    for (int i = 0; i < ctx->n_bindings; i++) {
        if (ctx->bindings[i].name && strcmp(ctx->bindings[i].name, name) == 0) {
            return &ctx->bindings[i];
        }
    }
    return NULL;
}

static const KflcExprFn *lookup_user_fn(const char *name, const KflcExprCtx *ctx)
{
    if (!ctx || !ctx->fns) return NULL;
    for (int i = 0; i < ctx->n_fns; i++) {
        if (ctx->fns[i].name && strcmp(ctx->fns[i].name, name) == 0) {
            return &ctx->fns[i];
        }
    }
    return NULL;
}

/* Vector-arg compute builtins. Each entry's `cxx_template` is split
 * by `$1` / `$2` etc. into a printf-style emission. The table is
 * small and hardcoded; growing the surface is a matter of adding
 * entries. The expected arg kind matters: V = a vector binding
 * (struct lvalue), S = a scalar expression, M = a matrix binding. */
typedef struct {
    const char *kfl_name;
    const char *cxx_template;     /* simple replacement: $1, $2 = arg indices */
    const char *arg_kinds;        /* string of 'V' 'S' 'M' chars; length = arity */
} ComputeFn;

static const ComputeFn COMPUTE_FNS[] = {
    /* Vector statistics (scalar-returning). */
    { "mean",      "k26c_stats_mean($1.data, $1.n)",                   "V" },
    { "var",       "k26c_stats_var($1.data, $1.n, 1)",                 "V" },
    { "std",       "k26c_stats_std($1.data, $1.n, 1)",                 "V" },
    { "min_v",     "k26c_stats_min($1.data, $1.n)",                    "V" },
    { "max_v",     "k26c_stats_max($1.data, $1.n)",                    "V" },
    { "pearson",   "k26c_stats_pearson($1.data, $2.data, $1.n)",       "VV" },
    /* Vector linear algebra. */
    { "dot",       "k26c_vec_dot(&($1), &($2))",                       "VV" },
    { "norm",      "k26c_vec_norm(&($1))",                             "V" },
    /* Matrix scalar-returning. mat_det is wrapped via a helper emitted
     * once in emit.c when needed; we just call that helper here. */
    { "mat_det",   "_kfl_mat_det(&($1))",                              "M" },
    /* Matrix out-parameter ops. Used as expression-statements:
     *   mat_mul(out, a, b)
     *   mat_vec(out_vec, a, x)
     *   mat_inv(out, in)
     *   transpose(out, in)
     *   solve(out_x, a, b)
     * The "return type" is K26CStatus but the surrounding (void) cast
     * in the expression-statement emit makes that irrelevant. */
    { "mat_mul",   "k26c_mat_mul(&($1), &($2), &($3))",                "MMM" },
    { "mat_vec",   "k26c_mat_vec(&($1), &($2), &($3))",                "VMV" },
    { "mat_inv",   "k26c_mat_inverse(&($1), &($2))",                   "MM" },
    { "transpose", "k26c_mat_transpose(&($1), &($2))",                 "MM" },
    { "solve",     "k26c_mat_solve(&($1), &($2), &($3))",              "VMV" },
    { "vec_add",   "k26c_vec_add(&($1), &($2), &($3))",                "VVV" },
    { "vec_sub",   "k26c_vec_sub(&($1), &($2), &($3))",                "VVV" },
    { "vec_scale", "k26c_vec_scale(&($1), ($2))",                      "VS" },
    { NULL, NULL, NULL }
};

static const ComputeFn *lookup_compute_fn(const char *name)
{
    if (!name) return NULL;
    for (const ComputeFn *c = COMPUTE_FNS; c->kfl_name; c++) {
        if (strcmp(c->kfl_name, name) == 0) return c;
    }
    return NULL;
}

static int emit_expr_rec(FILE *out, const KflcExpr *e,
                          const KflcExprCtx *ctx, KflcDiag *diag);

/* Emit a compute fn call. The arg expressions must conform to the
 * arg_kinds string: 'V' / 'M' args require an in-scope vector /
 * matrix binding (by identifier; no anonymous temporaries). Scalar
 * args (S) are emitted via the recursive expression emitter. */
static int emit_compute_call(FILE *out, const ComputeFn *cf,
                              const KflcExpr *call,
                              const KflcExprCtx *ctx, KflcDiag *diag)
{
    int n_args = (int)strlen(cf->arg_kinds);
    if (call->u.call.n_args != n_args) {
        kflc_diag_errorf(diag, call->line,
                         "expr: `%s` expects %d arg(s), got %d",
                         cf->kfl_name, n_args, call->u.call.n_args);
        return 1;
    }
    /* Walk the template, replacing $N with the rendered Nth argument. */
    for (const char *p = cf->cxx_template; *p; p++) {
        if (*p == '$' && p[1] >= '1' && p[1] <= '9') {
            int idx = p[1] - '1';
            p++;
            if (idx >= n_args) {
                kflc_diag_errorf(diag, call->line,
                                 "expr: template ref $%d out of range", idx + 1);
                return 1;
            }
            char want = cf->arg_kinds[idx];
            const KflcExpr *a = call->u.call.args[idx];
            if (want == 'V' || want == 'M') {
                /* Must be an in-scope ident of matching type. */
                if (a->kind != KFLE_IDENT) {
                    kflc_diag_errorf(diag, call->line,
                                     "expr: `%s` arg %d must be a %s identifier",
                                     cf->kfl_name, idx + 1,
                                     want == 'V' ? "vector" : "matrix");
                    return 1;
                }
                const KflcExprBinding *b = lookup_binding(a->u.ident, ctx);
                KflcType expected = (want == 'V') ? KFLT_VECTOR : KFLT_MATRIX;
                if (!b || b->type != expected) {
                    kflc_diag_errorf(diag, call->line,
                                     "expr: `%s` arg %d: `%s` is not a %s in scope",
                                     cf->kfl_name, idx + 1, a->u.ident,
                                     want == 'V' ? "vector" : "matrix");
                    return 1;
                }
                fputs(a->u.ident, out);
            } else {
                /* Scalar: emit as a parenthesised double expression. */
                fputc('(', out);
                if (emit_expr_rec(out, a, ctx, diag)) return 1;
                fputc(')', out);
            }
        } else {
            fputc(*p, out);
        }
    }
    return 0;
}

static int emit_expr_rec(FILE *out, const KflcExpr *e,
                          const KflcExprCtx *ctx, KflcDiag *diag)
{
    if (!e) return 1;
    switch (e->kind) {
    case KFLE_INT_LIT:
        fprintf(out, "(double)(%ld)", e->u.i);
        return 0;
    case KFLE_FLOAT_LIT:
        /* %.17g gives full double precision round-trip. */
        fprintf(out, "(%.17g)", e->u.f);
        return 0;
    case KFLE_VEC_LIT:
        kflc_diag_errorf(diag, e->line,
                         "expr: vector/matrix literal not allowed in scalar context");
        return 1;
    case KFLE_INDEX: {
        /* Two supported shapes:
         *   v[i]      vector element read (KFLC_VEC_R wraps it for
         *             opt-in bounds checking)
         *   m[i][j]   matrix cell read (KFLC_MAT_R wraps it for
         *             opt-in bounds checking; the macro does the
         *             two-compare row/col guard).
         * Anything else is an error. */
        const KflcExpr *base = e->u.index.base;
        if (!base) {
            kflc_diag_errorf(diag, e->line, "expr: index missing base");
            return 1;
        }
        if (base->kind == KFLE_IDENT) {
            const KflcExprBinding *b = lookup_binding(base->u.ident, ctx);
            if (!b) {
                kflc_diag_errorf(diag, e->line,
                                 "expr: unknown identifier `%s`",
                                 base->u.ident);
                return 1;
            }
            if (b->type == KFLT_MATRIX) {
                kflc_diag_errorf(diag, e->line,
                    "expr: matrix `%s` requires `[row][col]` not just `[row]`",
                    base->u.ident);
                return 1;
            }
            if (b->type != KFLT_VECTOR) {
                kflc_diag_errorf(diag, e->line,
                                 "expr: `%s` is not a vector in scope",
                                 base->u.ident);
                return 1;
            }
            fprintf(out, "KFLC_VEC_R(%s, ", base->u.ident);
            if (emit_expr_rec(out, e->u.index.idx, ctx, diag)) return 1;
            fputs(")", out);
            return 0;
        }
        /* Matrix m[i][j]: base is itself a KFLE_INDEX whose own
         * base is the matrix identifier. */
        if (base->kind == KFLE_INDEX &&
            base->u.index.base &&
            base->u.index.base->kind == KFLE_IDENT)
        {
            const char *mname = base->u.index.base->u.ident;
            const KflcExprBinding *b = lookup_binding(mname, ctx);
            if (!b || b->type != KFLT_MATRIX) {
                kflc_diag_errorf(diag, e->line,
                    "expr: `%s[i][j]` requires `%s` to be a matrix in scope",
                    mname, mname);
                return 1;
            }
            fprintf(out, "KFLC_MAT_R(%s, ", mname);
            if (emit_expr_rec(out, base->u.index.idx, ctx, diag)) return 1;
            fputs(", ", out);
            if (emit_expr_rec(out, e->u.index.idx, ctx, diag)) return 1;
            fputs(")", out);
            return 0;
        }
        kflc_diag_errorf(diag, e->line,
            "expr: index chains deeper than `m[i][j]` are not supported");
        return 1;
    }
    case KFLE_IDENT: {
        const char *nm = e->u.ident;
        const KflcExprBinding *b = lookup_binding(nm, ctx);
        if (b) {
            /* Use-after-move check (intra-block).
             *
             * If this binding was previously moved-from via move(<name>)
             * or own-to-own assignment, any subsequent read is a compile
             * error. Cross-block flow analysis is not performed;
             * intra-block coverage catches the common case.
             *
             * The check fires BEFORE the type dispatch so all paths
             * (opaque, scalar, vector) are protected. */
            if (b->moved_from) {
                kflc_diag_errorf(diag, e->line,
                    "expr: use after move of `%s` (binding was moved earlier "
                    "in this scope)",
                    nm);
                return 1;
            }
            /* Borrow source-alive check. A borrow binding's source
             * must not be moved-from at borrow read time; otherwise
             * the borrow points at NULL / freed memory. */
            if (b->lifetime_qualifier == KFL_LQ_BORROW
                && b->borrow_source_idx >= 0
                && b->borrow_source_idx < ctx->n_bindings
                && ctx->bindings[b->borrow_source_idx].moved_from) {
                kflc_diag_errorf(diag, e->line,
                    "expr: borrow `%s` read after its source `%s` was "
                    "moved (the borrow now points at NULL)",
                    nm, ctx->bindings[b->borrow_source_idx].name);
                return 1;
            }
            if (b->type == KFLT_VECTOR || b->type == KFLT_MATRIX) {
                /* Bare lvalue — caller (typically a compute-fn template)
                 * uses the struct fields directly. Emitting a vector
                 * inside a scalar expression context is a user error
                 * the compute-call emitter catches before we get here. */
                fputs(nm, out);
                return 0;
            }
            const char *pfx = b->is_form_arg ? "kfl_arg_" : "";
            if (b->type == KFLT_STRING) {
                /* String bindings (form-args only) emit as the bare
                 * C global. Consumed by builtins like dir_count(path)
                 * that take const char *. */
                fprintf(out, "%s%s", pfx, nm);
                return 0;
            }
            if (b->type == KFLT_BOOL || b->type == KFLT_INT) {
                /* Bool / int bindings cast to double for scalar
                 * context. Form-args use `kfl_arg_<name>`; local
                 * lets use the bare identifier. */
                fprintf(out, "((double)(%s%s))", pfx, nm);
                return 0;
            }
            if (b->type == KFLT_DOUBLE) {
                /* Double bindings — form-arg or local-let — read
                 * directly without a cast. */
                fprintf(out, "(%s%s)", pfx, nm);
                return 0;
            }
            if (b->type == KFLT_OPAQUE) {
                /* Opaque handles emit as the bare identifier so
                 * dynamic-builtin call sites receive the pointer
                 * value at full width. The receiving builtin
                 * signature carries the correct pointer type via
                 * the .kflbi manifest. */
                fprintf(out, "%s%s", pfx, nm);
                return 0;
            }
            kflc_diag_errorf(diag, e->line,
                             "expr: identifier `%s` has unsupported KFL type %d",
                             nm, (int)b->type);
            return 1;
        }
        kflc_diag_errorf(diag, e->line,
                         "expr: unknown identifier `%s`", nm);
        return 1;
    }
    case KFLE_UNARY: {
        const char *op = "?";
        switch (e->u.un.op) {
        case KFLOP_NEG: op = "-"; break;
        case KFLOP_POS: op = "+"; break;
        case KFLOP_NOT: op = "!"; break;
        default: break;
        }
        fputs("(", out);
        fputs(op, out);
        if (emit_expr_rec(out, e->u.un.operand, ctx, diag)) return 1;
        fputs(")", out);
        return 0;
    }
    case KFLE_BINARY: {
        const char *op = "?";
        switch (e->u.bin.op) {
        case KFLOP_ADD: op = "+"; break;
        case KFLOP_SUB: op = "-"; break;
        case KFLOP_MUL: op = "*"; break;
        case KFLOP_DIV: op = "/"; break;
        case KFLOP_MOD: op = ""; break;   /* handled below */
        case KFLOP_LT:  op = "<";  break;
        case KFLOP_LE:  op = "<="; break;
        case KFLOP_GT:  op = ">";  break;
        case KFLOP_GE:  op = ">="; break;
        case KFLOP_EQ:  op = "=="; break;
        case KFLOP_NE:  op = "!="; break;
        case KFLOP_AND: op = "&&"; break;
        case KFLOP_OR:  op = "||"; break;
        default: break;
        }
        if (e->u.bin.op == KFLOP_MOD) {
            /* C's `%` only works on integers; for doubles emit fmod(). */
            fputs("fmod(", out);
            if (emit_expr_rec(out, e->u.bin.lhs, ctx, diag)) return 1;
            fputs(", ", out);
            if (emit_expr_rec(out, e->u.bin.rhs, ctx, diag)) return 1;
            fputs(")", out);
            return 0;
        }
        fputs("(", out);
        if (emit_expr_rec(out, e->u.bin.lhs, ctx, diag)) return 1;
        fprintf(out, " %s ", op);
        if (emit_expr_rec(out, e->u.bin.rhs, ctx, diag)) return 1;
        fputs(")", out);
        return 0;
    }
    case KFLE_CALL: {
        /* `move(<ident>)` builtin.
         *
         * Explicit ownership transfer: takes a single own-qualified
         * identifier, marks its binding as moved_from = 1 (so
         * subsequent reads in the same scope are compile errors),
         * and emits a GCC statement-expression that yields the
         * pointer value and nulls the source.
         *
         * Required for own-to-own transfers across fn calls
         * (foo(move(x))) and return statements (return move(x)).
         * Auto-move on own-to-own let assignment is also handled in
         * stmt.c by detecting bare own-binding RHS, but the explicit
         * form is the recommended user-facing pattern. */
        if (e->u.call.name && strcmp(e->u.call.name, "move") == 0) {
            if (e->u.call.n_args != 1) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move() requires exactly one identifier argument");
                return 1;
            }
            const KflcExpr *arg = e->u.call.args[0];
            if (arg->kind != KFLE_IDENT) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move() argument must be a bare identifier "
                    "(an own-qualified binding name)");
                return 1;
            }
            const KflcExprBinding *src = lookup_binding(arg->u.ident, ctx);
            if (!src) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move(`%s`): unknown identifier",
                    arg->u.ident);
                return 1;
            }
            if (src->moved_from) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move(`%s`): binding was already moved-from "
                    "earlier in this scope",
                    arg->u.ident);
                return 1;
            }
            if (src->lifetime_qualifier != KFL_LQ_OWN &&
                src->lifetime_qualifier != KFL_LQ_PTR) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move(`%s`): binding must be `own` or `ptr` "
                    "qualified (got %s)",
                    arg->u.ident,
                    kflc_lifetime_qualifier_kfl_str(src->lifetime_qualifier));
                return 1;
            }
            if (src->type != KFLT_OPAQUE) {
                kflc_diag_errorf(diag, e->line,
                    "expr: move(`%s`): only opaque-typed bindings can "
                    "be moved (scalars are value-copy)",
                    arg->u.ident);
                return 1;
            }
            /* Mark source as moved. Cast away const — the bindings
             * array is logically mutable; the const-on-pointer is a
             * "don't realloc" hint, not a "don't update fields"
             * promise (see KflcExprCtx struct comment). */
            ((KflcExprBinding *)src)->moved_from = 1;
            const char *cxx_ty = kflc_opaque_cxx(src->type_subtype);
            if (!cxx_ty) cxx_ty = "void *";
            const char *pfx = src->is_form_arg ? "kfl_arg_" : "";
            /* GCC statement-expr: capture the pointer value, null the
             * source, return the captured pointer. K26 targets gcc
             * everywhere; portability is fine. Note: kflc_opaque_cxx
             * returns the full pointer type (e.g. "K26AstroBody *"),
             * so no extra `*` here. */
            fprintf(out,
                "({ %s _kfl_mv = %s%s; %s%s = NULL; _kfl_mv; })",
                cxx_ty, pfx, arg->u.ident, pfx, arg->u.ident);
            return 0;
        }

        const BuiltinFn *b = lookup_builtin(e->u.call.name);
        if (b) {
            if (b->arity >= 0 && b->arity != e->u.call.n_args) {
                kflc_diag_errorf(diag, e->line,
                                 "expr: `%s` expects %d arg(s), got %d",
                                 e->u.call.name, b->arity, e->u.call.n_args);
                return 1;
            }
            /* String helpers: emit-path specials. */
            if (is_string_builtin_(b->kfl_name)) {
                if (strcmp(b->kfl_name, "concat") == 0) {
                    /* Per-callsite static char[256] via GCC statement
                     * expression. Each textual occurrence gets its
                     * own block-scope static so nested / sibling
                     * concats don't stomp each other. */
                    fputs("({ static char _kfl_cbuf[256]; "
                          "k26_str_concat(_kfl_cbuf, sizeof _kfl_cbuf, ",
                          out);
                    if (emit_expr_rec(out, e->u.call.args[0], ctx, diag)) return 1;
                    fputs(", ", out);
                    if (emit_expr_rec(out, e->u.call.args[1], ctx, diag)) return 1;
                    fputs("); (const char *)_kfl_cbuf; })", out);
                    return 0;
                }
                /* Scalar-returning string builtins — wrap in
                 * (double) so they slot into a scalar expression
                 * context. */
                fprintf(out, "((double)%s(", b->cxx_name);
                for (int i = 0; i < e->u.call.n_args; i++) {
                    if (i > 0) fputs(", ", out);
                    if (emit_expr_rec(out, e->u.call.args[i], ctx, diag)) return 1;
                }
                fputs("))", out);
                return 0;
            }
            fprintf(out, "%s(", b->cxx_name);
            for (int i = 0; i < e->u.call.n_args; i++) {
                if (i > 0) fputs(", ", out);
                if (emit_expr_rec(out, e->u.call.args[i], ctx, diag)) return 1;
            }
            fputs(")", out);
            return 0;
        }
        /* libk26compute vector/matrix builtins. */
        const ComputeFn *cf = lookup_compute_fn(e->u.call.name);
        if (cf) return emit_compute_call(out, cf, e, ctx, diag);
        /* Not a builtin — try user-declared fns. */
        const KflcExprFn *uf = lookup_user_fn(e->u.call.name, ctx);
        if (uf) {
            if (uf->arity >= 0 && uf->arity != e->u.call.n_args) {
                kflc_diag_errorf(diag, e->line,
                                 "expr: `%s` expects %d arg(s), got %d",
                                 e->u.call.name, uf->arity, e->u.call.n_args);
                return 1;
            }
            fprintf(out, "%s(", e->u.call.name);
            for (int i = 0; i < e->u.call.n_args; i++) {
                if (i > 0) fputs(", ", out);
                if (emit_expr_rec(out, e->u.call.args[i], ctx, diag)) return 1;
            }
            fputs(")", out);
            return 0;
        }
        kflc_diag_errorf(diag, e->line,
                         "expr: unknown function `%s`", e->u.call.name);
        return 1;
    }
    }
    return 1;
}

int kflc_emit_expr(FILE *out, const KflcExpr *e,
                   const KflcExprCtx *ctx, KflcDiag *diag)
{
    if (!out || !e) {
        if (diag) kflc_diag_errorf(diag, 0, "emit_expr: NULL input");
        return 1;
    }
    return emit_expr_rec(out, e, ctx, diag);
}

/* Parse-time-known body name check. */
int kflc_body_idx_known(const KflcExprCtx *ctx, const char *name)
{
    if (!ctx || !name) return 0;
    if (!ctx->known_body_names || ctx->n_known_body_names <= 0) return 0;
    for (int i = 0; i < ctx->n_known_body_names; i++) {
        const char *kn = ctx->known_body_names[i];
        if (kn && strcmp(kn, name) == 0) return 1;
    }
    return 0;
}

/* Emit a KflcExpr in lvalue position. See kflc.h for the contract.
 * Each shape mirrors the corresponding rvalue emit (KFLE_IDENT and
 * KFLE_INDEX cases above) but produces an assignable expression:
 *  - scalar idents skip the (double) cast emit_expr_rec applies
 *  - vector index uses KFLC_VEC_W (the lvalue-returning lambda) so
 *    -DKFLC_BOUNDS_CHECK still guards writes
 *  - matrix double-index uses KFLC_MAT_W/R so bounds-check builds
 *    catch OOB on both row and col */
int kflc_emit_lvalue(FILE *out, const KflcExpr *e,
                     const KflcExprCtx *ctx, KflcDiag *diag,
                     KflcLvMode mode)
{
    if (!out || !e) {
        if (diag) kflc_diag_errorf(diag, 0, "emit_lvalue: NULL input");
        return 1;
    }
    /* The walker dispatches on mode for vector/matrix access.
     * Scalar IDENT path is mode-independent. WRITE callers (assign
     * cases in stmt.c) get KFLC_VEC_W / KFLC_MAT_W (returning
     * `double&`); READ callers route bounds-checked rvalue access
     * through KFLC_VEC_R / KFLC_MAT_R. The macros collapse to plain
     * `.data[...]` in non-bounds-check builds (emit.c preamble). */
    const char *vec_macro = (mode == KFLC_LV_WRITE) ? "KFLC_VEC_W" : "KFLC_VEC_R";
    const char *mat_macro = (mode == KFLC_LV_WRITE) ? "KFLC_MAT_W" : "KFLC_MAT_R";
    if (e->kind == KFLE_IDENT) {
        /* Scalar lvalue. Form-arg writes prefix `kfl_arg_`; locals
         * are bare. No double-cast (we're on the LHS of `=`). */
        const KflcExprBinding *b = lookup_binding(e->u.ident, ctx);
        if (b && b->type == KFLT_VECTOR) {
            if (diag) kflc_diag_errorf(diag, e->line,
                "lvalue: cannot assign to whole vector `%s` "
                "(use `%s[i] = ...` for element write)",
                e->u.ident, e->u.ident);
            return 1;
        }
        if (b && b->type == KFLT_MATRIX) {
            if (diag) kflc_diag_errorf(diag, e->line,
                "lvalue: cannot assign to whole matrix `%s` "
                "(use `%s[r][c] = ...` for cell write)",
                e->u.ident, e->u.ident);
            return 1;
        }
        const char *pfx = (b && b->is_form_arg) ? "kfl_arg_" : "";
        fprintf(out, "%s%s", pfx, e->u.ident);
        return 0;
    }
    if (e->kind == KFLE_INDEX) {
        const KflcExpr *base = e->u.index.base;
        if (!base) {
            if (diag) kflc_diag_errorf(diag, e->line,
                "lvalue: index missing base");
            return 1;
        }
        if (base->kind == KFLE_IDENT) {
            /* v[i] vector write. Look up to confirm it really is a
             * vector; give a typed error if the user wrote `m[i]` on
             * a matrix (must use `m[i][j]`). */
            const KflcExprBinding *b = lookup_binding(base->u.ident, ctx);
            if (b && b->type == KFLT_MATRIX) {
                if (diag) kflc_diag_errorf(diag, e->line,
                    "lvalue: matrix `%s` write requires `[row][col]` "
                    "not just `[row]`", base->u.ident);
                return 1;
            }
            if (!b || b->type != KFLT_VECTOR) {
                if (diag) kflc_diag_errorf(diag, e->line,
                    "lvalue: `%s` is not a vector in scope",
                    base->u.ident);
                return 1;
            }
            fprintf(out, "%s(%s, ", vec_macro, base->u.ident);
            if (kflc_emit_expr(out, e->u.index.idx, ctx, diag)) return 1;
            fputs(")", out);
            return 0;
        }
        /* m[i][j] matrix write: base is itself an INDEX whose base
         * is the matrix identifier. */
        if (base->kind == KFLE_INDEX &&
            base->u.index.base &&
            base->u.index.base->kind == KFLE_IDENT)
        {
            const char *mname = base->u.index.base->u.ident;
            const KflcExprBinding *b = lookup_binding(mname, ctx);
            if (!b || b->type != KFLT_MATRIX) {
                if (diag) kflc_diag_errorf(diag, e->line,
                    "lvalue: `%s[i][j]` requires `%s` to be a matrix in scope",
                    mname, mname);
                return 1;
            }
            /* Route matrix access through KFLC_MAT_R/W so bounds-
             * check builds catch OOB on both row and col. */
            fprintf(out, "%s(%s, ", mat_macro, mname);
            if (kflc_emit_expr(out, base->u.index.idx, ctx, diag)) return 1;
            fputs(", ", out);
            if (kflc_emit_expr(out, e->u.index.idx, ctx, diag)) return 1;
            fputs(")", out);
            return 0;
        }
        if (diag) kflc_diag_errorf(diag, e->line,
            "lvalue: index chains deeper than `m[i][j]` are not supported");
        return 1;
    }
    if (diag) kflc_diag_errorf(diag, e->line,
        "lvalue: expected identifier or index chain");
    return 1;
}

/* ---- Text serialization (inverse of parse_expr) ----------------- */

static const char *op_text(KflcOp op)
{
    switch (op) {
    case KFLOP_NEG: return "-";
    case KFLOP_POS: return "+";
    case KFLOP_NOT: return "!";
    case KFLOP_ADD: return "+";
    case KFLOP_SUB: return "-";
    case KFLOP_MUL: return "*";
    case KFLOP_DIV: return "/";
    case KFLOP_MOD: return "%";
    case KFLOP_LT:  return "<";
    case KFLOP_LE:  return "<=";
    case KFLOP_GT:  return ">";
    case KFLOP_GE:  return ">=";
    case KFLOP_EQ:  return "==";
    case KFLOP_NE:  return "!=";
    case KFLOP_AND: return "&&";
    case KFLOP_OR:  return "||";
    }
    return "?";
}

void kflc_expr_to_text(FILE *out, const KflcExpr *e)
{
    if (!out || !e) return;
    switch (e->kind) {
    case KFLE_INT_LIT:
        fprintf(out, "%ld", e->u.i);
        return;
    case KFLE_FLOAT_LIT: {
        /* %.17g preserves the double's bits round-trip-safe, but for
         * values that round to a whole number (1.0, -2.0) it drops the
         * trailing decimal point and the reparser picks it up as an
         * int literal — losing the FLOAT_LIT discriminant. Force a
         * `.0` suffix in that case so round-tripping preserves kind. */
        char buf[64];
        int  n = snprintf(buf, sizeof(buf), "%.17g", e->u.f);
        int  has_dot_or_exp = 0;
        for (int i = 0; i < n; i++) {
            if (buf[i] == '.' || buf[i] == 'e' || buf[i] == 'E' ||
                buf[i] == 'n' || buf[i] == 'N') {
                /* '.'/'e' for finite floats; 'n'/'N' for nan/inf */
                has_dot_or_exp = 1;
                break;
            }
        }
        fputs(buf, out);
        if (!has_dot_or_exp) fputs(".0", out);
        return;
    }
    case KFLE_IDENT:
        fputs(e->u.ident ? e->u.ident : "?", out);
        return;
    case KFLE_UNARY:
        fputc('(', out);
        fputs(op_text(e->u.un.op), out);
        kflc_expr_to_text(out, e->u.un.operand);
        fputc(')', out);
        return;
    case KFLE_BINARY:
        fputc('(', out);
        kflc_expr_to_text(out, e->u.bin.lhs);
        fprintf(out, " %s ", op_text(e->u.bin.op));
        kflc_expr_to_text(out, e->u.bin.rhs);
        fputc(')', out);
        return;
    case KFLE_CALL:
        fputs(e->u.call.name ? e->u.call.name : "?", out);
        fputc('(', out);
        for (int i = 0; i < e->u.call.n_args; i++) {
            if (i > 0) fputs(", ", out);
            kflc_expr_to_text(out, e->u.call.args[i]);
        }
        fputc(')', out);
        return;
    case KFLE_VEC_LIT:
        fputc('[', out);
        for (int i = 0; i < e->u.vec.n_elems; i++) {
            if (i > 0) fputs(", ", out);
            kflc_expr_to_text(out, e->u.vec.elems[i]);
        }
        fputc(']', out);
        return;
    case KFLE_INDEX:
        kflc_expr_to_text(out, e->u.index.base);
        fputc('[', out);
        kflc_expr_to_text(out, e->u.index.idx);
        fputc(']', out);
        return;
    }
}
