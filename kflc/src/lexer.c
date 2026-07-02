/* kflc — lexer. Hand-rolled; produces Tokens for the parser. */

#include "internal.h"

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void lex_init(Lexer *L, const char *src, size_t len,
              KflcArena *arena, KflcDiag *diag)
{
    L->src            = src;
    L->pos            = 0;
    L->len            = len;
    L->line           = 1;
    L->line_start_pos = 0;
    L->arena          = arena;
    L->diag           = diag;
}

static int peek(Lexer *L, size_t off)
{
    size_t p = L->pos + off;
    if (p >= L->len) return -1;
    return (unsigned char)L->src[p];
}

static int hex_value(int c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static int is_ident_start(int c)
{
    return (c == '_') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

static int is_ident_cont(int c)
{
    return is_ident_start(c) || (c >= '0' && c <= '9');
}

/* Skip horizontal whitespace + line-continuation comments (`#` to EOL,
 * but only if `#` isn't the start of a #rrggbb colour literal). */
static void skip_ws_and_comments(Lexer *L)
{
    for (;;) {
        int c = peek(L, 0);
        if (c == ' ' || c == '\t' || c == '\r') {
            L->pos++;
            continue;
        }
        if (c == '#') {
            /* Distinguish colour literal from comment. A colour is `#`
             * followed by exactly 6 hex digits and then a non-alnum
             * boundary. Otherwise it's a comment. */
            int is_colour = 1;
            for (int i = 1; i <= 6; i++) {
                int d = peek(L, i);
                if (hex_value(d) < 0) { is_colour = 0; break; }
            }
            if (is_colour) {
                int tail = peek(L, 7);
                if (tail >= 0 && (is_ident_cont(tail) || tail == '#')) {
                    is_colour = 0;
                }
            }
            if (is_colour) return;   /* colour literal; let the main lexer handle */
            while (L->pos < L->len && L->src[L->pos] != '\n') L->pos++;
            continue;
        }
        return;
    }
}

static int read_string(Lexer *L, Token *tok)
{
    int line = L->line;
    L->pos++;  /* consume opening " */
    /* Decode into a scratch buffer (over-allocated in the arena). */
    size_t start = L->pos;
    size_t cap = 64;
    char *buf = (char *)kflc_arena_alloc(L->arena, cap);
    size_t bn = 0;

    while (L->pos < L->len) {
        int c = (unsigned char)L->src[L->pos];
        if (c == '"') {
            L->pos++;
            char *out = (char *)kflc_arena_alloc(L->arena, bn + 1);
            memcpy(out, buf, bn);
            out[bn] = '\0';
            tok->kind = T_STRING;
            tok->line = line;
            tok->str  = out;
            return 1;
        }
        if (c == '\n') {
            kflc_diag_errorf(L->diag, line,
                "unterminated string literal");
            tok->kind = T_EOF;
            return 0;
        }
        if (c == '\\') {
            int esc = peek(L, 1);
            int decoded;
            switch (esc) {
                case 'n':  decoded = '\n'; break;
                case 't':  decoded = '\t'; break;
                case 'r':  decoded = '\r'; break;
                case '"':  decoded = '"';  break;
                case '\\': decoded = '\\'; break;
                default:
                    kflc_diag_errorf(L->diag, L->line,
                        "unknown escape \\%c in string", esc);
                    tok->kind = T_EOF;
                    return 0;
            }
            if (bn + 1 >= cap) {
                size_t newcap = cap * 2;
                char *nb = (char *)kflc_arena_alloc(L->arena, newcap);
                memcpy(nb, buf, bn);
                buf = nb;
                cap = newcap;
            }
            buf[bn++] = (char)decoded;
            L->pos += 2;
            continue;
        }
        if (bn + 1 >= cap) {
            size_t newcap = cap * 2;
            char *nb = (char *)kflc_arena_alloc(L->arena, newcap);
            memcpy(nb, buf, bn);
            buf = nb;
            cap = newcap;
        }
        buf[bn++] = (char)c;
        L->pos++;
    }
    kflc_diag_errorf(L->diag, line, "unterminated string literal");
    tok->kind = T_EOF;
    (void)start;
    return 0;
}

static int read_number_or_size(Lexer *L, Token *tok)
{
    int line = L->line;
    size_t start = L->pos;
    int negative = 0;
    if (L->src[L->pos] == '-') { negative = 1; L->pos++; }
    size_t int_start = L->pos;
    while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
    if (L->pos == int_start) {
        /* `-` without digits — not a number. Rewind. */
        L->pos = start;
        return 0;
    }

    /* WxH ? */
    if (L->pos < L->len && L->src[L->pos] == 'x' &&
        L->pos + 1 < L->len && isdigit((unsigned char)L->src[L->pos + 1]))
    {
        size_t w_end = L->pos;
        char wbuf[16];
        size_t wlen = w_end - int_start;
        if (negative || wlen >= sizeof(wbuf)) {
            kflc_diag_errorf(L->diag, line, "invalid size literal");
            tok->kind = T_EOF;
            return 0;
        }
        memcpy(wbuf, L->src + int_start, wlen);
        wbuf[wlen] = '\0';

        L->pos++;  /* consume 'x' */
        size_t h_start = L->pos;
        while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
        size_t hlen = L->pos - h_start;
        char hbuf[16];
        if (hlen == 0 || hlen >= sizeof(hbuf)) {
            kflc_diag_errorf(L->diag, line, "invalid size literal");
            tok->kind = T_EOF;
            return 0;
        }
        memcpy(hbuf, L->src + h_start, hlen);
        hbuf[hlen] = '\0';

        tok->kind = T_SIZE;
        tok->line = line;
        tok->w    = atoi(wbuf);
        tok->h    = atoi(hbuf);
        return 1;
    }

    /* Float? */
    if (L->pos < L->len && L->src[L->pos] == '.' &&
        L->pos + 1 < L->len && isdigit((unsigned char)L->src[L->pos + 1]))
    {
        L->pos++;
        while (L->pos < L->len && isdigit((unsigned char)L->src[L->pos])) L->pos++;
        char buf[64];
        size_t n = L->pos - start;
        if (n >= sizeof(buf)) {
            kflc_diag_errorf(L->diag, line, "float literal too long");
            tok->kind = T_EOF;
            return 0;
        }
        memcpy(buf, L->src + start, n);
        buf[n] = '\0';
        tok->kind = T_FLOAT;
        tok->line = line;
        tok->f    = strtod(buf, NULL);
        return 1;
    }

    /* Plain int. */
    char buf[32];
    size_t n = L->pos - start;
    if (n >= sizeof(buf)) {
        kflc_diag_errorf(L->diag, line, "integer literal too long");
        tok->kind = T_EOF;
        return 0;
    }
    memcpy(buf, L->src + start, n);
    buf[n] = '\0';
    tok->kind = T_INT;
    tok->line = line;
    tok->i    = strtol(buf, NULL, 10);
    return 1;
}

static int read_ident(Lexer *L, Token *tok)
{
    int line = L->line;
    size_t start = L->pos;
    while (L->pos < L->len && is_ident_cont((unsigned char)L->src[L->pos])) {
        L->pos++;
    }
    size_t n = L->pos - start;
    tok->kind = T_IDENT;
    tok->line = line;
    tok->str  = kflc_arena_strndup(L->arena, L->src + start, n);
    return 1;
}

static int read_color(Lexer *L, Token *tok)
{
    int line = L->line;
    /* skip_ws_and_comments has already confirmed 6 hex digits follow. */
    L->pos++;  /* '#' */
    uint32_t v = 0;
    for (int i = 0; i < 6; i++) {
        int d = hex_value((unsigned char)L->src[L->pos++]);
        v = (v << 4) | (uint32_t)d;
    }
    tok->kind = T_COLOR;
    tok->line = line;
    tok->rgb  = v;
    return 1;
}

static int read_arrow(Lexer *L, Token *tok)
{
    int line = L->line;
    /* '-' '>' */
    if (L->pos + 1 < L->len &&
        L->src[L->pos] == '-' && L->src[L->pos + 1] == '>') {
        L->pos += 2;
        tok->kind = T_ARROW;
        tok->line = line;
        return 1;
    }
    return 0;
}

/* UTF-8 right-arrow → is E2 86 92. */
static int read_unicode_arrow(Lexer *L, Token *tok)
{
    int line = L->line;
    if (L->pos + 2 < L->len &&
        (unsigned char)L->src[L->pos]     == 0xE2 &&
        (unsigned char)L->src[L->pos + 1] == 0x86 &&
        (unsigned char)L->src[L->pos + 2] == 0x92)
    {
        L->pos += 3;
        tok->kind = T_ARROW;
        tok->line = line;
        return 1;
    }
    return 0;
}

int lex_next(Lexer *L, Token *tok)
{
    skip_ws_and_comments(L);

    /* Column of the first byte of whatever token we're about to read.
     * Computed once per call after ws-skip so every successful return
     * path can splat it onto tok->col with no duplication. */
    int col = (int)(L->pos - L->line_start_pos) + 1;

    if (L->pos >= L->len) {
        tok->kind = T_EOF;
        tok->line = L->line;
        tok->col  = col;
        return 1;
    }

    int c = (unsigned char)L->src[L->pos];

    if (c == '\n') {
        L->pos++;
        /* Update line_start_pos here, not on the next call's ws-skip:
         * the column reported for a T_NEWLINE refers to the position
         * of the '\n' itself, but the NEXT token starts at the byte
         * after the newline, so the running line_start_pos must
         * advance now. */
        L->line_start_pos = L->pos;
        tok->kind = T_NEWLINE;
        tok->line = L->line;
        tok->col  = col;
        L->line++;
        return 1;
    }

    if (c == '"') {
        if (!read_string(L, tok)) return 0;
        tok->col = col;
        return 1;
    }
    if (c == '#') {
        if (!read_color(L, tok)) return 0;
        tok->col = col;
        return 1;
    }

    if (c == '-') {
        if (read_arrow(L, tok))            { tok->col = col; return 1; }
        if (!read_number_or_size(L, tok))  return 0;
        tok->col = col;
        return 1;
    }

    if (c == 0xE2) {
        if (read_unicode_arrow(L, tok))    { tok->col = col; return 1; }
        /* Otherwise an unrecognised non-ASCII byte. Treat as identifier
         * character (UTF-8 inside identifiers is not supported in v1;
         * but inside a "raw word" run e.g. for keyboard symbols, allow). */
    }

    if (isdigit(c)) {
        if (!read_number_or_size(L, tok)) return 0;
        tok->col = col;
        return 1;
    }

    if (is_ident_start(c)) {
        if (!read_ident(L, tok)) return 0;
        tok->col = col;
        return 1;
    }

    kflc_diag_errorf_col(L->diag, L->line, col,
        "unexpected character: '%c' (0x%02x)",
        (c >= 32 && c < 127) ? c : '?', c);
    L->pos++;
    tok->kind = T_EOF;
    tok->col  = col;
    return 0;
}

int lex_next_non_nl(Lexer *L, Token *tok)
{
    do {
        if (!lex_next(L, tok)) return 0;
    } while (tok->kind == T_NEWLINE);
    return 1;
}

/* Consume the next non-whitespace run (stopping at whitespace or
 * newline) as an opaque string. */
static int lex_take_word_like(Lexer *L, Token *tok, TokenKind out_kind)
{
    /* Skip horizontal whitespace but DO NOT cross a newline. */
    while (L->pos < L->len) {
        int c = (unsigned char)L->src[L->pos];
        if (c == ' ' || c == '\t' || c == '\r') { L->pos++; continue; }
        break;
    }
    int line = L->line;
    int col  = (int)(L->pos - L->line_start_pos) + 1;
    if (L->pos >= L->len || L->src[L->pos] == '\n') {
        kflc_diag_errorf_col(L->diag, line, col, "expected value");
        tok->kind = T_EOF;
        tok->col  = col;
        return 0;
    }
    size_t start = L->pos;
    while (L->pos < L->len) {
        int c = (unsigned char)L->src[L->pos];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') break;
        L->pos++;
    }
    size_t n = L->pos - start;
    tok->kind = out_kind;
    tok->line = line;
    tok->col  = col;
    tok->str  = kflc_arena_strndup(L->arena, L->src + start, n);
    return 1;
}

int lex_take_shortcut(Lexer *L, Token *tok)
{
    return lex_take_word_like(L, tok, T_SHORTCUT);
}

int lex_take_word(Lexer *L, Token *tok)
{
    return lex_take_word_like(L, tok, T_IDENT);
}
