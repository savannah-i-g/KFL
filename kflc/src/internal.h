#ifndef KFLC_INTERNAL_H
#define KFLC_INTERNAL_H

#include "kflc.h"

#include <stddef.h>
#include <stdint.h>

typedef enum {
    T_EOF = 0,
    T_NEWLINE,
    T_IDENT,
    T_STRING,
    T_INT,
    T_FLOAT,
    T_SIZE,        /* WxH (int w, int h) */
    T_COLOR,       /* #rrggbb (uint32_t rgb) */
    T_ARROW,       /* "->" or "→" */
    T_SHORTCUT     /* opaque token consumed via lex_take_shortcut() */
} TokenKind;

typedef struct {
    TokenKind  kind;
    int        line;
    int        col;        /* 1-indexed column of the first byte */
    char      *str;        /* arena-owned: T_IDENT, T_STRING, T_SHORTCUT */
    long       i;          /* T_INT */
    double     f;          /* T_FLOAT */
    int        w, h;       /* T_SIZE */
    uint32_t   rgb;        /* T_COLOR */
} Token;

typedef struct {
    const char *src;
    size_t      pos;
    size_t      len;
    int         line;
    /* Byte offset of the first character on the current line. Derived
     * column = (pos - line_start_pos) + 1. Updated when the lexer
     * advances past a '\n'. Lets the diagnostic system emit
     * editor-friendly `<path>:<line>:<col>:` locations without paying
     * for per-byte column counting. */
    size_t      line_start_pos;
    KflcArena  *arena;
    KflcDiag   *diag;
} Lexer;

void  lex_init(Lexer *L, const char *src, size_t len,
               KflcArena *arena, KflcDiag *diag);
int   lex_next(Lexer *L, Token *tok);

/* lex_take_shortcut consumes the next non-whitespace run (stopping at
 * whitespace or newline) as a single shortcut literal. Use after a
 * statement keyword like `shortcut` that expects a free-form key combo. */
int   lex_take_shortcut(Lexer *L, Token *tok);

/* lex_take_word consumes the next non-whitespace run as an identifier-
 * shaped token. Same shape as lex_take_shortcut but emits T_IDENT. Used
 * for attributes whose value is an arbitrary identifier the parser
 * cannot pre-lex (rare). */
int   lex_take_word(Lexer *L, Token *tok);

/* Skip newlines (return next non-newline token in *tok). */
int   lex_next_non_nl(Lexer *L, Token *tok);

/* Fn-body statement parser + emitter (stmt.c). */

/* Parse statements off `L` (with `cur` as the lookahead Token), until
 * we hit `terminator` or any of the `also_break` keywords (NULL-
 * terminated array, or NULL itself). The terminator token is NOT
 * consumed — caller decides. Returns a parent placeholder node whose
 * `children` is the parsed statement list. */
KflcNode *kfl_parse_stmt_block(Lexer *L, Token *cur,
                                KflcArena *arena, KflcDiag *diag,
                                int *had_error,
                                const char *terminator,
                                const char **also_break);

/* Emit a single statement subtree to `out`. Returns 0 on success. */
int       kfl_emit_stmt(FILE *out, const KflcNode *s,
                         const KflcExprCtx *ctx, KflcDiag *diag,
                         int indent);

/* Returns 1 if `name` matches a keyword reserved for a future KFL
 * grammar version (thread / join / mutex / lock / unlock /
 * parallel_for / bind / import). Callers should pair the check with a
 * `kflc_diag_warnf` so authors get a deprecation-style heads-up before
 * the keyword starts being enforced. Used by the `arg` parser and the
 * `let`/`const` parser. */
int       kfl_is_reserved_future(const char *name);

/* Block-scope tracker for heap-typed `let`s. The fn-body root level
 * (depth 0) is folded in too, so this tracker is the single source of
 * truth for "what heap-typed locals are live and how do we free them
 * when a scope exits or a return fires."
 *
 * Lifecycle: emit.c calls `kfl_emit_stmt_reset_scopes(tmp_arena,
 * fn_return_type)` at the top of each fn-body emit loop. The
 * return-type is passed in so the RETURN handler can save a
 * non-void return value into a temporary BEFORE freeing locals
 * (avoiding use-after-free of returned heap content; see the
 * `_kfl_rv` pattern). Then `kfl_emit_stmt` manages push/pop and
 * per-block free emission internally as it walks the AST. emit.c
 * must call `kfl_emit_stmt_drain_root(out, indent)` at the natural
 * end of the fn body (when control falls through, not via return)
 * to clean up depth-0 lets. */
void      kfl_emit_stmt_reset_scopes(KflcArena *arena,
                                      KflcType   fn_return_type);
/* Opaque-subtype-aware variant. Pass NULL for fn_return_subtype when
 * fn_return_type isn't KFLT_OPAQUE. The non-_ex form is kept as a
 * thin wrapper so existing call sites continue to work; new code
 * (e.g. `fn world` emission) should call this. */
void      kfl_emit_stmt_reset_scopes_ex(KflcArena  *arena,
                                         KflcType    fn_return_type,
                                         const char *fn_return_subtype);
void      kfl_emit_stmt_drain_root(FILE *out, int indent);

#endif /* KFLC_INTERNAL_H */
