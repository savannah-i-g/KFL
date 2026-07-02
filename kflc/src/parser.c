/* kflc — recursive-descent parser. Consumes Lexer tokens, builds AST. */

#include "internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    Lexer       L;
    Token       cur;
    KflcArena  *arena;
    KflcDiag   *diag;
    int         had_error;
} Parser;

/* ---- Helpers ----------------------------------------------------- */

static void advance(Parser *P)
{
    if (!lex_next(&P->L, &P->cur)) P->had_error = 1;
}

static int at_eof(Parser *P) { return P->cur.kind == T_EOF; }
static int at_nl (Parser *P) { return P->cur.kind == T_NEWLINE; }

static int is_ident_named(const Token *t, const char *name)
{
    return t->kind == T_IDENT && t->str && strcmp(t->str, name) == 0;
}

static void skip_newlines(Parser *P)
{
    while (at_nl(P)) advance(P);
}

static void expect_nl_or_eof(Parser *P, const char *ctx)
{
    if (at_nl(P) || at_eof(P)) {
        if (at_nl(P)) advance(P);
        return;
    }
    kflc_diag_errorf(P->diag, P->cur.line,
        "expected end of line after %s (got token kind %d)", ctx, P->cur.kind);
    P->had_error = 1;
    /* recovery: drain rest of line */
    while (!at_nl(P) && !at_eof(P)) advance(P);
    if (at_nl(P)) advance(P);
}

static KflcNode *new_node(Parser *P, KflcNodeKind kind, int line)
{
    KflcNode *n = (KflcNode *)kflc_arena_alloc(P->arena, sizeof(*n));
    n->kind = kind;
    n->line = line;
    return n;
}

static void append_child(KflcNode *parent, KflcNode *child)
{
    if (!parent->children) {
        parent->children = child;
        return;
    }
    KflcNode *p = parent->children;
    while (p->next) p = p->next;
    p->next = child;
}

static KflcAttr *new_attr(Parser *P, const char *name, int line)
{
    KflcAttr *a = (KflcAttr *)kflc_arena_alloc(P->arena, sizeof(*a));
    a->name = kflc_arena_strdup(P->arena, name);
    a->line = line;
    return a;
}

static void append_attr(KflcNode *target, KflcAttr *a)
{
    if (!target->attrs) {
        target->attrs = a;
        return;
    }
    KflcAttr *p = target->attrs;
    while (p->next) p = p->next;
    p->next = a;
}

/* Convert current Token into a KflcValue. Does NOT advance. */
static KflcValue tok_to_value(const Token *t)
{
    KflcValue v;
    memset(&v, 0, sizeof(v));
    switch (t->kind) {
    case T_INT:      v.kind = KFLV_INT;     v.u.i = t->i; break;
    case T_FLOAT:    v.kind = KFLV_FLOAT;   v.u.f = t->f; break;
    case T_STRING:   v.kind = KFLV_STR;     v.u.s = t->str; break;
    case T_IDENT:    v.kind = KFLV_IDENT;   v.u.s = t->str; break;
    case T_COLOR:    v.kind = KFLV_COLOR;   v.u.rgb = t->rgb; break;
    case T_SIZE:     v.kind = KFLV_SIZE;    v.u.size.w = t->w; v.u.size.h = t->h; break;
    case T_SHORTCUT: v.kind = KFLV_SHORTCUT; v.u.s = t->str; break;
    default:         v.kind = KFLV_NONE; break;
    }
    return v;
}

/* ---- Attribute parsing ------------------------------------------- */

/* Names of statement keywords that close a block. */
static int is_block_terminator(const Token *t)
{
    return is_ident_named(t, "end");
}

/* The only leaf widget in the artifact language is `plot`. */
static int is_widget_keyword(const Token *t)
{
    return is_ident_named(t, "plot");
}

/* Attributes are parsed by consuming the rest of the line. The
 * `target` node accumulates them as `<key> <value>` pairs. */
static void parse_attrs_until_nl(Parser *P, KflcNode *target)
{
    while (!at_nl(P) && !at_eof(P)) {
        if (P->cur.kind != T_IDENT) {
            kflc_diag_errorf(P->diag, P->cur.line,
                "expected attribute keyword");
            P->had_error = 1;
            while (!at_nl(P) && !at_eof(P)) advance(P);
            return;
        }
        char *key = P->cur.str;
        int   kln = P->cur.line;

        advance(P);

        /* Default: consume one value token. */
        KflcValue v = tok_to_value(&P->cur);
        if (v.kind == KFLV_NONE) {
            kflc_diag_errorf(P->diag, kln,
                "%s: expected value, got token kind %d", key, P->cur.kind);
            P->had_error = 1;
            while (!at_nl(P) && !at_eof(P)) advance(P);
            return;
        }
        advance(P);
        KflcAttr *a = new_attr(P, key, kln);
        a->value = v;
        append_attr(target, a);
    }
}

/* ---- Block bodies ------------------------------------------------- */

/* parse_widget_or_container — current token is the `plot` keyword.
 * Returns the corresponding KFLN_WIDGET node, or NULL on error. The
 * only leaf the artifact language recognises is the plot figure. */
static KflcNode *parse_widget_or_container(Parser *P)
{
    if (is_widget_keyword(&P->cur)) {
        int line = P->cur.line;
        advance(P);

        KflcNode *n = new_node(P, KFLN_WIDGET, line);
        n->widget = KFL_W_PLOT;

        /* Optional positional first value: string (label) or ident (name). */
        if (P->cur.kind == T_STRING || P->cur.kind == T_IDENT) {
            /* Don't consume an ident that is actually the next attribute
             * keyword. Distinguish by: positional value is allowed only
             * IMMEDIATELY after the widget keyword, before any attrs.
             * We treat the first STRING as a label; the first IDENT is
             * a name only if it isn't a known attribute keyword.
             *
             * Heuristic: if the ident matches a known attribute keyword
             * (on, paint, label, value, etc.), it's an attribute. Else
             * it's the positional name.
             */
            int treat_as_position = 1;
            if (P->cur.kind == T_IDENT) {
                static const char *attr_keywords[] = {
                    "on","data","paint","label","value","color",
                    "label_color","text_color","box","down_box","font",
                    "font_size","text_size","align","at","size",
                    "width","height","weight","tooltip","divider","toggle",
                    "shortcut","icon","min","max","step","thumb","orientation",
                    /* plot3d / scene3d attributes */
                    "scene","camera_pos","camera_target","background","fov",
                    "near","far","up","show_axes","wireframe","backface_cull",
                    "ssaa","nav_target","show_breadcrumb","animate_flag",
                    "select_index_arg","select_flag_arg",
                    "focus_x_arg","focus_y_arg","focus_z_arg",
                    /* runtime-toggle / command / dir-pane attrs */
                    "bind_visible","bind_hidden","target","cmd","sep",
                    /* instance / per-instance buffers */
                    "xforms","tints","navs",NULL
                };
                for (int i = 0; attr_keywords[i]; i++) {
                    if (strcmp(P->cur.str, attr_keywords[i]) == 0) {
                        treat_as_position = 0;
                        break;
                    }
                }
            }
            if (treat_as_position) {
                if (P->cur.kind == T_IDENT) n->name = P->cur.str;
                n->position = tok_to_value(&P->cur);
                advance(P);
            }
        } else if (P->cur.kind == T_INT) {
            /* Numeric positional — unusual but allowed (e.g. label 5). */
            n->position = tok_to_value(&P->cur);
            advance(P);
        }

        /* Attributes on the same line, OR continuation lines that begin
         * with an attribute keyword. The grammar doc shows multi-line
         * widget descriptors like:
         *
         *     output display
         *         height     56
         *         text_size  24
         *
         * To support that, after the first line's attrs we continue
         * eating attribute lines that DO NOT start with a known
         * widget/container/block-terminator keyword. */
        parse_attrs_until_nl(P, n);
        if (at_nl(P)) advance(P);

        for (;;) {
            /* Save position so we can rewind if it's not a continuation. */
            size_t saved_pos = P->L.pos;
            int    saved_ln  = P->L.line;
            Token  saved_cur = P->cur;
            skip_newlines(P);

            if (at_eof(P)) break;
            if (is_block_terminator(&P->cur)) break;
            if (P->cur.kind != T_IDENT) break;

            /* In continuation context, attribute keywords win over
             * widget/container keywords of the same name — EXCEPT for
             * the two cases (`box` and `label`) where both readings
             * are syntactically valid. For those, peek at the next
             * token: `box` followed by a known box style is an attr;
             * `label` followed by a string is an attr. Anything else
             * is a new widget statement (rewind). */
            static const char *attr_keywords[] = {
                "on","data","paint","label","value","color",
                "label_color","text_color","box","down_box","font",
                /* box_style / text are the new canonical names for
                 * the attribute spellings of `box` / `label`. The
                 * old names are accepted with a deprecation warning
                 * (see parse_attrs_until_nl); the new names are
                 * unambiguous (no widget kinds by those names).
                 * Widget-kind keywords `box` and `label` are
                 * unaffected by this rename. */
                "box_style","text",
                "font_size","text_size","align","at","size",
                "width","height","weight","tooltip","divider","toggle",
                "shortcut","icon","min","max","step","thumb","orientation",
                /* plot / compute */
                "kind","title","x_label","y_label","show_grid","show_legend",
                "x_min","x_max","y_min","y_max","expression","format",
                /* plot3d / scene3d */
                "scene","camera_pos","camera_target","background","fov",
                "near","far","up","show_axes","wireframe","backface_cull",
                "ssaa","nav_target","show_breadcrumb","animate_flag",
                "select_index_arg","select_flag_arg",
                "focus_x_arg","focus_y_arg","focus_z_arg",
                /* runtime-toggle / command / dir-pane attrs */
                "bind_visible","bind_hidden","target","cmd","sep",
                /* per-instance buffers */
                "xforms","tints","navs",NULL
            };
            int is_attr = 0;
            for (int i = 0; attr_keywords[i]; i++) {
                if (strcmp(P->cur.str, attr_keywords[i]) == 0) {
                    is_attr = 1; break;
                }
            }

            /* Peek-and-decide for the ambiguous ones. We need to look
             * at the token AFTER the keyword without permanently
             * consuming it. Snapshot full state, advance one token,
             * inspect, then rewind. */
            if (is_attr && (strcmp(P->cur.str, "box") == 0 ||
                            strcmp(P->cur.str, "label") == 0))
            {
                int is_box = (strcmp(P->cur.str, "box") == 0);
                size_t peek_saved_pos = P->L.pos;
                int    peek_saved_ln  = P->L.line;
                Token  peek_saved_cur = P->cur;
                advance(P);
                int reads_as_attr = 0;
                if (is_box) {
                    /* Attribute usage: `box <known_box_style_ident>`. */
                    static const char *box_styles[] = {
                        "flat","up","down","thin_up","thin_down","no_box",
                        "plastic_thin_up","plastic_thin_down","plastic_up","plastic_down",
                        NULL
                    };
                    if (P->cur.kind == T_IDENT) {
                        for (int i = 0; box_styles[i]; i++) {
                            if (strcmp(P->cur.str, box_styles[i]) == 0) {
                                reads_as_attr = 1; break;
                            }
                        }
                    }
                } else {
                    /* `label` as attribute requires a string value. */
                    if (P->cur.kind == T_STRING) reads_as_attr = 1;
                }
                /* Rewind the peek. */
                P->L.pos  = peek_saved_pos;
                P->L.line = peek_saved_ln;
                P->cur    = peek_saved_cur;
                if (!reads_as_attr) is_attr = 0;
            }

            if (!is_attr) {
                /* Rewind to saved state — next iteration of the parent
                 * parser will dispatch on the unconsumed token. */
                P->L.pos  = saved_pos;
                P->L.line = saved_ln;
                P->cur    = saved_cur;
                break;
            }
            parse_attrs_until_nl(P, n);
            if (at_nl(P)) advance(P);
        }

        return n;
    }

    kflc_diag_errorf(P->diag, P->cur.line,
        "expected widget or container keyword, got '%s'",
        P->cur.kind == T_IDENT ? P->cur.str : "(non-ident)");
    P->had_error = 1;
    while (!at_nl(P) && !at_eof(P)) advance(P);
    if (at_nl(P)) advance(P);
    return NULL;
}

/* ---- Form-level statements -------------------------------------- */

static void parse_form_attr(Parser *P, KflcNode *form, const char *name)
{
    int line = P->cur.line;
    advance(P);  /* consume keyword */
    if (P->cur.kind == T_EOF || at_nl(P)) {
        kflc_diag_errorf(P->diag, line, "%s: expected value", name);
        P->had_error = 1;
        return;
    }
    KflcValue v = tok_to_value(&P->cur);
    if (v.kind == KFLV_NONE) {
        kflc_diag_errorf(P->diag, line, "%s: invalid value", name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        if (at_nl(P)) advance(P);
        return;
    }
    advance(P);
    KflcAttr *a = new_attr(P, name, line);
    a->value = v;
    append_attr(form, a);
    expect_nl_or_eof(P, name);
}

/* `tick <handler_name> interval_ms <N>`: registers a periodic
 * timer callback at form scope. KFL grammar style: positional
 * identifier is the handler name; attribute keyword `interval_ms` (no
 * `=`) is followed by an integer literal. Emit wires the registration
 * into the existing kflc_plot3d_set_tick_cb seam. */
static void parse_tick(Parser *P, KflcNode *form)
{
    int line = P->cur.line;
    advance(P);   /* consume `tick` */
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line,
            "tick: expected handler identifier");
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcNode *n = new_node(P, KFLN_TICK, line);
    n->name = P->cur.str;
    advance(P);
    parse_attrs_until_nl(P, n);
    append_child(form, n);
    expect_nl_or_eof(P, "tick");
}

/* `frame <name> inertial` or `frame <name> body_fixed body <ident>`.
 * Declares a new frame newtype at form scope; the emitter registers
 * an opaque type whose KFL name is `<name>`. Inside fn bodies, a
 * `let p: <name>` declaration uses the frame as a type tag — values
 * tagged with one frame can't be passed to functions expecting a
 * different frame without an explicit `to_frame(...)` conversion. */
static void parse_frame(Parser *P, KflcNode *form)
{
    int line = P->cur.line;
    advance(P);   /* consume `frame` */
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line,
            "frame: expected frame name identifier");
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcNode *n = new_node(P, KFLN_FRAME, line);
    n->name = P->cur.str;
    advance(P);
    if (P->cur.kind != T_IDENT ||
        (strcmp(P->cur.str, "inertial")  != 0 &&
         strcmp(P->cur.str, "body_fixed") != 0)) {
        kflc_diag_errorf(P->diag, line,
            "frame %s: expected `inertial` or `body_fixed`", n->name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    /* Stash the kind in attrs[kind]. */
    KflcAttr *kind_a = new_attr(P, "kind", line);
    kind_a->value.kind = KFLV_IDENT;
    kind_a->value.u.s = P->cur.str;
    append_attr(n, kind_a);
    advance(P);
    parse_attrs_until_nl(P, n);   /* picks up `body <ident>` if present */
    append_child(form, n);
    expect_nl_or_eof(P, "frame");
}

/* `arena <name> capacity <int>[KB|MB|GB]` declares a form-level
 * bump allocator. Lives for the lifetime of the form; fn bodies
 * opt into the arena via `allocator <name>`. The `capacity` value
 * is parsed as a byte count, with optional unit suffix:
 *   1024, 64 KB, 16 MB, 1 GB; all decoded to bytes at parse time.
 * Optional trailing `reset_mode <fn|manual|form>` attr controls when
 * the arena resets; defaults to `fn` (reset on each fn entry that
 * binds it).
 *
 * Syntax note: KFL form-scope attributes are space-separated
 * positional pairs (no `=`), matching `tick fn=foo` style and
 * `frame X body_fixed body Y`. The form `capacity = 64MB` is
 * accepted at fn-body scope via the let/const = pathway when an
 * arena is locally allocated; the form-scope syntax is positional. */
static void parse_arena(Parser *P, KflcNode *form)
{
    int line = P->cur.line;
    advance(P);   /* consume `arena` */
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line,
            "arena: expected arena name identifier");
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcNode *n = new_node(P, KFLN_ARENA, line);
    n->name = P->cur.str;
    advance(P);
    if (P->cur.kind != T_IDENT || strcmp(P->cur.str, "capacity") != 0) {
        kflc_diag_errorf(P->diag, line,
            "arena %s: expected `capacity <bytes>`", n->name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    advance(P);   /* consume `capacity` */
    if (P->cur.kind != T_INT) {
        kflc_diag_errorf(P->diag, line,
            "arena %s: expected integer byte count after `capacity`", n->name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    long bytes = P->cur.i;
    /* Optional unit suffix on the next token: KB / MB / GB. */
    advance(P);
    if (P->cur.kind == T_IDENT && P->cur.str) {
        if      (strcmp(P->cur.str, "KB") == 0) { bytes *= 1024L;            advance(P); }
        else if (strcmp(P->cur.str, "MB") == 0) { bytes *= 1024L*1024L;      advance(P); }
        else if (strcmp(P->cur.str, "GB") == 0) { bytes *= 1024L*1024L*1024L;advance(P); }
        /* Else leave the ident to parse_attrs_until_nl below — it
         * may be `reset_mode` or similar. */
    }
    KflcAttr *cap_a = new_attr(P, "capacity", line);
    cap_a->value.kind = KFLV_INT;
    cap_a->value.u.i  = bytes;
    append_attr(n, cap_a);
    parse_attrs_until_nl(P, n);   /* picks up optional `reset_mode <kind>` */
    append_child(form, n);
    expect_nl_or_eof(P, "arena");
}

/* `epoch <name> "<ISO-8601>" <scale>`: declares a named compile-
 * time epoch constant. The ISO-8601 string lives in attrs["iso"]
 * (KFLV_STR); the scale ident in attrs["scale"]. The emitter parses
 * the ISO string at compile time and emits the
 * (int64 days_since_J2000, double seconds_of_day, uint8 scale)
 * triple as a K26AstroEpoch literal. */
static void parse_epoch(Parser *P, KflcNode *form)
{
    int line = P->cur.line;
    advance(P);   /* consume `epoch` */
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line,
            "epoch: expected binding name");
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcNode *n = new_node(P, KFLN_EPOCH_LITERAL, line);
    n->name = P->cur.str;
    advance(P);
    if (P->cur.kind != T_STRING) {
        kflc_diag_errorf(P->diag, line,
            "epoch %s: expected ISO-8601 string literal", n->name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcAttr *iso_a = new_attr(P, "iso", line);
    iso_a->value.kind = KFLV_STR;
    iso_a->value.u.s  = P->cur.str;
    append_attr(n, iso_a);
    advance(P);
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line,
            "epoch %s: expected time-scale identifier (TAI/UTC/UT1/TT/TDB)",
            n->name);
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        return;
    }
    KflcAttr *sc_a = new_attr(P, "scale", line);
    sc_a->value.kind = KFLV_IDENT;
    sc_a->value.u.s  = P->cur.str;
    append_attr(n, sc_a);
    advance(P);
    append_child(form, n);
    expect_nl_or_eof(P, "epoch");
}

/* v2: arg <name> [default <value>] */
static void parse_arg(Parser *P, KflcNode *form)
{
    int line = P->cur.line;
    advance(P);   /* consume 'arg' */
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, line, "arg: expected name identifier");
        P->had_error = 1;
        while (!at_nl(P) && !at_eof(P)) advance(P);
        if (at_nl(P)) advance(P);
        return;
    }
    KflcNode *a = new_node(P, KFLN_ARG, line);
    a->name = P->cur.str;
    if (kfl_is_reserved_future(a->name)) {
        kflc_diag_warnf(P->diag, line,
            "arg `%s` shadows a name reserved for a future KFL keyword; "
            "consider renaming to avoid breakage when that keyword lands",
            a->name);
    }
    advance(P);

    /* Optional `default <value>`. */
    if (is_ident_named(&P->cur, "default")) {
        advance(P);
        if (at_nl(P) || at_eof(P)) {
            kflc_diag_errorf(P->diag, line, "arg %s: `default` requires a value", a->name);
            P->had_error = 1;
        } else {
            KflcValue v = tok_to_value(&P->cur);
            /* Accept idents `true` / `false` as bool defaults. */
            if (v.kind == KFLV_IDENT && v.u.s) {
                if (strcmp(v.u.s, "true")  == 0) { v.kind = KFLV_INT; v.u.i = 1; }
                else if (strcmp(v.u.s, "false") == 0) { v.kind = KFLV_INT; v.u.i = 0; }
            }
            a->position = v;
            advance(P);
        }
    } else {
        /* Default: bool false. */
        a->position.kind = KFLV_INT;
        a->position.u.i  = 0;
    }

    /* Optional trailing mutability modifier.
     *   arg <name> [default <value>] [mutable | readonly]
     * Default is `mutable`. `readonly` marks the arg so the emit:
     *  - prefixes the global with `const`
     *  - rejects any `target <arg>` / `cmd toggle target <arg>` that
     *    would write to it. */
    if (is_ident_named(&P->cur, "readonly")) {
        a->flags |= KFL_NF_READONLY;
        advance(P);
    } else if (is_ident_named(&P->cur, "mutable")) {
        /* explicit mutable: default behaviour, just consume the token */
        advance(P);
    }

    append_child(form, a);
    expect_nl_or_eof(P, "arg");
}

/* ---- Top-level form ---------------------------------------------- */

static KflcNode *parse_form(Parser *P)
{
    skip_newlines(P);
    if (!is_ident_named(&P->cur, "form")) {
        kflc_diag_errorf(P->diag, P->cur.line, "expected `form` at top of file");
        P->had_error = 1;
        return NULL;
    }
    int form_line = P->cur.line;
    advance(P);
    if (P->cur.kind != T_IDENT) {
        kflc_diag_errorf(P->diag, form_line, "form: expected name identifier");
        P->had_error = 1;
        return NULL;
    }
    KflcNode *form = new_node(P, KFLN_FORM, form_line);
    form->name = P->cur.str;
    advance(P);
    expect_nl_or_eof(P, "form opener");

    for (;;) {
        skip_newlines(P);
        if (at_eof(P)) {
            kflc_diag_errorf(P->diag, form_line, "missing `end` for `form %s`", form->name);
            P->had_error = 1;
            return form;
        }
        if (is_block_terminator(&P->cur)) {
            advance(P);
            expect_nl_or_eof(P, "end");
            return form;
        }
        if (P->cur.kind != T_IDENT) {
            kflc_diag_errorf(P->diag, P->cur.line,
                "expected statement keyword inside form");
            P->had_error = 1;
            while (!at_nl(P) && !at_eof(P)) advance(P);
            if (at_nl(P)) advance(P);
            continue;
        }

        const char *kw = P->cur.str;
        if      (strcmp(kw, "title")    == 0) parse_form_attr(P, form, "title");
        else if (strcmp(kw, "size")     == 0) parse_form_attr(P, form, "size");
        else if (strcmp(kw, "cfg")      == 0) parse_form_attr(P, form, "cfg");
        /* Top-level plot declaration: `plot <name> data <fn> [title ...]
         * [x_label ...] ...`. Reuses the plot widget parse; artifact
         * emission renders each to <name>.png + <name>.svg. */
        else if (strcmp(kw, "plot")     == 0) {
            KflcNode *w = parse_widget_or_container(P);
            if (w) append_child(form, w);
        }
        else if (strcmp(kw, "tick")     == 0) parse_tick    (P, form);
        else if (strcmp(kw, "frame")    == 0) parse_frame   (P, form);
        else if (strcmp(kw, "epoch")    == 0) parse_epoch   (P, form);
        else if (strcmp(kw, "arena")    == 0) parse_arena   (P, form);
        else if (strcmp(kw, "arg")      == 0) parse_arg(P, form);
        else if (strcmp(kw, "fn")       == 0) {
            /* Inline function definition.
             *   fn <return_type> <name>(<typed_args>) NEWLINE
             *     <statements>
             *   end
             *
             * Data-callback shape:
             *   fn data <name> NEWLINE
             *     <statements>
             *   end
             * Emits a function with the libk26plot data-callback C-ABI
             * signature; body uses series_<kind> to populate the
             * out-parameters.
             *
             * The header line is consumed via the main token stream
             * (identifiers + the wordy ident form `name(arg1, arg2)`).
             * Since the lexer doesn't tokenise `(` `)` `,`, we capture
             * the parenthesised arg list by reading raw line bytes.
             */
            int  fn_line = P->cur.line;
            advance(P);   /* consume 'fn' */

            if (P->cur.kind != T_IDENT) {
                kflc_diag_errorf(P->diag, fn_line, "fn: expected return type");
                P->had_error = 1;
                while (!at_nl(P) && !at_eof(P)) advance(P);
                if (at_nl(P)) advance(P);
                continue;
            }

            /* `fn world <name>`: callback-shaped fn for an
             * astro/sim runtime. The body is plain statements (let /
             * const / assign / expression-stmt) plus calls to
             * registered astro builtins; codegen emits an
             * `extern "C" void <name>(<world_cxx> world, void *_kfl_user_data)`
             * where <world_cxx> comes from the "world" opaque type
             * registry. Round-trip works without "world" being
             * registered (parse + serialize are type-agnostic); full
             * emit requires the registration. */
            if (is_ident_named(&P->cur, "world")) {
                advance(P);   /* consume 'world' */
                if (P->cur.kind != T_IDENT) {
                    kflc_diag_errorf(P->diag, fn_line,
                        "fn world: expected function name");
                    P->had_error = 1;
                    while (!at_nl(P) && !at_eof(P)) advance(P);
                    if (at_nl(P)) advance(P);
                    continue;
                }
                char *fn_name = P->cur.str;
                advance(P);
                expect_nl_or_eof(P, "fn world opener");

                KflcNode *fnode = new_node(P, KFLN_FN_WORLD, fn_line);
                fnode->name = fn_name;
                fnode->type = KFLT_VOID;

                const char *brk[] = { "end", NULL };
                int had_inner = 0;
                KflcNode *body = kfl_parse_stmt_block(&P->L, &P->cur,
                                                       P->arena, P->diag, &had_inner,
                                                       "end", brk);
                if (had_inner) P->had_error = 1;
                if (body && body->children) {
                    fnode->children = body->children;
                }
                if (is_ident_named(&P->cur, "end")) {
                    advance(P);
                    expect_nl_or_eof(P, "end");
                } else {
                    kflc_diag_errorf(P->diag, fn_line,
                                     "fn world %s: expected `end`", fn_name);
                    P->had_error = 1;
                }
                append_child(form, fnode);
                continue;
            }

            /* `fn data <name>`: callback-shaped fn for plot data. */
            if (is_ident_named(&P->cur, "data")) {
                advance(P);   /* consume 'data' */
                if (P->cur.kind != T_IDENT) {
                    kflc_diag_errorf(P->diag, fn_line,
                        "fn data: expected function name");
                    P->had_error = 1;
                    while (!at_nl(P) && !at_eof(P)) advance(P);
                    if (at_nl(P)) advance(P);
                    continue;
                }
                char *fn_name = P->cur.str;
                advance(P);
                expect_nl_or_eof(P, "fn data opener");

                KflcNode *fnode = new_node(P, KFLN_FN_DATA, fn_line);
                fnode->name = fn_name;
                fnode->type = KFLT_VOID;

                const char *brk[] = { "end", NULL };
                int had_inner = 0;
                KflcNode *body = kfl_parse_stmt_block(&P->L, &P->cur,
                                                       P->arena, P->diag, &had_inner,
                                                       "end", brk);
                if (had_inner) P->had_error = 1;
                if (body && body->children) {
                    fnode->children = body->children;
                }
                if (is_ident_named(&P->cur, "end")) {
                    advance(P);
                    expect_nl_or_eof(P, "end");
                } else {
                    kflc_diag_errorf(P->diag, fn_line,
                                     "fn data %s: expected `end`", fn_name);
                    P->had_error = 1;
                }
                append_child(form, fnode);
                continue;
            }
            /* `fn [<lifetime>] <type> <name>(...)`. If the first token
             * after `fn` names a lifetime qualifier, the type follows;
             * otherwise the first token IS the type. */
            KflcLifetimeQualifier ret_lq = KFL_LQ_NONE;
            int ret_lq_check = kflc_lifetime_qualifier_from_str(P->cur.str);
            if (ret_lq_check >= 0) {
                ret_lq = (KflcLifetimeQualifier)ret_lq_check;
                advance(P);
                if (P->cur.kind != T_IDENT) {
                    kflc_diag_errorf(P->diag, fn_line,
                                     "fn: expected return type after `%s`",
                                     kflc_lifetime_qualifier_kfl_str(ret_lq));
                    P->had_error = 1;
                    while (!at_nl(P) && !at_eof(P)) advance(P);
                    if (at_nl(P)) advance(P);
                    continue;
                }
            }
            const char *rt_opaque = NULL;
            int rt = kflc_type_from_str(P->cur.str, &rt_opaque);
            if (rt < 0) {
                kflc_diag_errorf(P->diag, fn_line,
                                 "fn: unknown return type `%s`", P->cur.str);
                P->had_error = 1;
                while (!at_nl(P) && !at_eof(P)) advance(P);
                if (at_nl(P)) advance(P);
                continue;
            }
            KflcType ret_type = (KflcType)rt;
            /* Capture the opaque-subtype string into the arena so it
             * survives the parser's token churn. */
            const char *ret_opaque = rt_opaque
                ? kflc_arena_strdup(P->arena, rt_opaque) : NULL;
            advance(P);

            if (P->cur.kind != T_IDENT) {
                kflc_diag_errorf(P->diag, fn_line,
                                 "fn: expected function name");
                P->had_error = 1;
                while (!at_nl(P) && !at_eof(P)) advance(P);
                if (at_nl(P)) advance(P);
                continue;
            }
            char *fn_name = P->cur.str;
            /* Take the rest of the line as the parenthesised arg list. */
            size_t arg_line_start = P->L.pos;
            while (P->L.pos < P->L.len && P->L.src[P->L.pos] != '\n') P->L.pos++;
            size_t arg_line_len = P->L.pos - arg_line_start;
            char *arg_buf = (char *)kflc_arena_alloc(P->arena, arg_line_len + 1);
            if (arg_line_len > 0) memcpy(arg_buf, P->L.src + arg_line_start, arg_line_len);
            arg_buf[arg_line_len] = '\0';
            /* Refill cur to be the NEWLINE that follows. */
            advance(P);   /* past the fn-name token */
            if (at_nl(P)) advance(P);

            KflcNode *fnode = new_node(P, KFLN_FN, fn_line);
            fnode->name = fn_name;
            fnode->type = ret_type;
            fnode->type_subtype = ret_opaque;
            fnode->lifetime_qualifier = ret_lq;

            /* Trim leading whitespace, expect `(`. */
            char *p = arg_buf;
            while (*p == ' ' || *p == '\t') p++;
            if (*p != '(') {
                kflc_diag_errorf(P->diag, fn_line,
                                 "fn %s: expected `(` after name", fn_name);
                P->had_error = 1;
                /* still attempt to parse a body so we don't desync */
            } else {
                p++;
                /* Parse comma-separated `[<lifetime>] <type> <name>`
                 * pairs until `)`. The lifetime qualifier is recognised
                 * before the type token; absence leaves the arg with
                 * KFL_LQ_NONE. */
                while (*p && *p != ')') {
                    while (*p == ' ' || *p == '\t') p++;
                    if (*p == ')') break;
                    /* First word may be `own` / `borrow` / `ptr`
                     * (lifetime qualifier) or the type itself. */
                    char *first_start = p;
                    while (*p && !(*p == ' ' || *p == '\t')) p++;
                    char saved = *p; *p = '\0';
                    KflcLifetimeQualifier arg_lq = KFL_LQ_NONE;
                    int arg_lq_check = kflc_lifetime_qualifier_from_str(first_start);
                    char *type_start;
                    if (arg_lq_check >= 0) {
                        arg_lq = (KflcLifetimeQualifier)arg_lq_check;
                        *p = saved;
                        while (*p == ' ' || *p == '\t') p++;
                        type_start = p;
                        while (*p && !(*p == ' ' || *p == '\t')) p++;
                        saved = *p; *p = '\0';
                    } else {
                        type_start = first_start;
                    }
                    const char *arg_opaque_raw = NULL;
                    int t = kflc_type_from_str(type_start, &arg_opaque_raw);
                    /* Copy out the opaque name BEFORE restoring the
                     * delimiter — type_start points back into the
                     * caller's arg_buf. */
                    const char *arg_opaque = arg_opaque_raw
                        ? kflc_arena_strdup(P->arena, arg_opaque_raw) : NULL;
                    *p = saved;
                    if (t < 0) {
                        kflc_diag_errorf(P->diag, fn_line,
                                         "fn %s: unknown arg type `%s`",
                                         fn_name, type_start);
                        P->had_error = 1;
                        break;
                    }
                    /* name token */
                    while (*p == ' ' || *p == '\t') p++;
                    char *name_start = p;
                    while (*p && !(*p == ' ' || *p == '\t' || *p == ',' || *p == ')')) p++;
                    char saved2 = *p; *p = '\0';
                    char *aname = kflc_arena_strdup(P->arena, name_start);
                    *p = saved2;
                    KflcNode *an = new_node(P, KFLN_FN_ARG, fn_line);
                    an->name = aname;
                    an->type = (KflcType)t;
                    an->type_subtype = arg_opaque;
                    an->lifetime_qualifier = arg_lq;
                    append_child(fnode, an);
                    /* Eat optional comma. */
                    while (*p == ' ' || *p == '\t') p++;
                    if (*p == ',') p++;
                }
            }

            /* Body: a sequence of statements until `end`. The statement
             * parser uses the same Lexer/Token pair. We capture them as
             * the local L/cur, hand off, then sync back. */
            const char *brk[] = { "end", NULL };
            int had_inner = 0;
            KflcNode *body = kfl_parse_stmt_block(&P->L, &P->cur,
                                                   P->arena, P->diag, &had_inner,
                                                   "end", brk);
            if (had_inner) P->had_error = 1;
            /* Steal the body's children list as fnode's children, AFTER
             * the arg list. Walk fnode->children to the tail then
             * splice. */
            if (body && body->children) {
                if (!fnode->children) {
                    fnode->children = body->children;
                } else {
                    KflcNode *tail = fnode->children;
                    while (tail->next) tail = tail->next;
                    tail->next = body->children;
                }
            }
            /* Consume `end` + newline. */
            if (is_ident_named(&P->cur, "end")) {
                advance(P);
                expect_nl_or_eof(P, "end");
            } else {
                kflc_diag_errorf(P->diag, fn_line,
                                 "fn %s: expected `end`", fn_name);
                P->had_error = 1;
            }
            append_child(form, fnode);
        }
        else {
            kflc_diag_errorf(P->diag, P->cur.line,
                "unknown form-level keyword `%s`", kw);
            P->had_error = 1;
            while (!at_nl(P) && !at_eof(P)) advance(P);
            if (at_nl(P)) advance(P);
        }
    }
}

/* ---- Public entry points ----------------------------------------- */

KflcNode *kflc_parse(const char *src, size_t len,
                     KflcArena *arena, KflcDiag *diag)
{
    Parser P;
    P.arena     = arena;
    P.diag      = diag;
    P.had_error = 0;
    lex_init(&P.L, src, len, arena, diag);
    advance(&P);
    KflcNode *form = parse_form(&P);
    if (P.had_error || (diag && diag->errors)) return NULL;
    return form;
}

KflcNode *kflc_parse_file(const char *path,
                          KflcArena *arena, KflcDiag *diag)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        kflc_diag_errorf(diag, 0, "cannot open `%s`", path);
        return NULL;
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz < 0) {
        fclose(f);
        kflc_diag_errorf(diag, 0, "cannot size `%s`", path);
        return NULL;
    }
    char *buf = (char *)kflc_arena_alloc(arena, (size_t)sz + 1);
    if (!buf) {
        fclose(f);
        kflc_diag_errorf(diag, 0, "out of memory");
        return NULL;
    }
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[got] = '\0';
    return kflc_parse(buf, got, arena, diag);
}
