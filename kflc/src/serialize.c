/* kflc — AST serialiser. Inverse of the parser: walks a parsed
 * KflcNode form and writes textual .kfl source back to a FILE *.
 *
 * Pure C. Depends only on libkflc's public types. Used by embedding
 * tools (k26form, k26canvas) that load a form via kflc_parse_file,
 * edit the AST in place, and save back to disk. */

#include "kflc.h"

#include <stdio.h>
#include <string.h>

static const char *widget_kw(KflcWidgetKind w)
{
    switch (w) {
    case KFL_W_PLOT:            return "plot";
    }
    return "plot";
}

static void indent(FILE *out, int n)
{
    for (int i = 0; i < n; i++) fputs("    ", out);
}

static void emit_string(FILE *out, const char *s)
{
    fputc('"', out);
    for (const char *p = s; p && *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (c == '"' || c == '\\') fputc('\\', out);
        if (c == '\n') { fputs("\\n", out); continue; }
        if (c == '\t') { fputs("\\t", out); continue; }
        fputc((char)c, out);
    }
    fputc('"', out);
}

static void emit_value(FILE *out, const KflcValue *v)
{
    switch (v->kind) {
    case KFLV_INT:      fprintf(out, "%ld", v->u.i); break;
    case KFLV_FLOAT: {
        /* Same discriminant-preserving trick as kflc_expr_to_text:
         * %.17g drops the trailing decimal for whole-valued doubles
         * (1.0 → "1") which reparses as KFLV_INT. Force a `.0`
         * suffix to keep the float kind round-trip-safe. */
        char buf[64];
        int  n = snprintf(buf, sizeof(buf), "%.17g", v->u.f);
        int  has = 0;
        for (int i = 0; i < n; i++) {
            if (buf[i] == '.' || buf[i] == 'e' || buf[i] == 'E' ||
                buf[i] == 'n' || buf[i] == 'N') { has = 1; break; }
        }
        fputs(buf, out);
        if (!has) fputs(".0", out);
        break;
    }
    case KFLV_STR:      emit_string(out, v->u.s ? v->u.s : ""); break;
    case KFLV_IDENT:    fputs(v->u.s ? v->u.s : "", out); break;
    case KFLV_COLOR:    fprintf(out, "#%06x", v->u.rgb); break;
    case KFLV_SIZE:     fprintf(out, "%dx%d", v->u.size.w, v->u.size.h); break;
    case KFLV_SHORTCUT: fputs(v->u.s ? v->u.s : "", out); break;
    case KFLV_ALIGN:
    case KFLV_NONE:
    default: break;
    }
}

/* Triple-attr key bases. Parser expands `<key> X Y Z` into three
 * KflcAttrs named `<key>_x`, `<key>_y`, `<key>_z` (stmt.c:472-475).
 * Serializer needs the inverse: when iterating attrs, detect the
 * `_x` member and emit `<key> X Y Z` as a group; skip the `_y` /
 * `_z` members because they're already covered. */
static const char *const TRIPLE_BASES_[] = {
    "position", "direction",
    "transform_translate", "transform_scale", "transform_rotate",
    NULL
};

static const char *triple_base_for_suffix_(const char *name, char want)
{
    size_t n = strlen(name);
    if (n < 3) return NULL;
    if (name[n - 2] != '_' || name[n - 1] != want) return NULL;
    for (int i = 0; TRIPLE_BASES_[i]; i++) {
        size_t bl = strlen(TRIPLE_BASES_[i]);
        if (n == bl + 2 && memcmp(name, TRIPLE_BASES_[i], bl) == 0)
            return TRIPLE_BASES_[i];
    }
    return NULL;
}

static void emit_attr_value_raw_(FILE *out, const KflcAttr *a)
{
    /* For triple-attr components the parser stores either a FLOAT
     * (pure number) or an IDENT (verbatim source for in-scope refs).
     * Emit each in the form the parser will accept back. */
    if (a->value.kind == KFLV_FLOAT) {
        /* Match the kflc_expr_to_text FLOAT_LIT pattern — keep a
         * decimal point so reparse preserves the float discriminant. */
        char buf[64];
        int  n = snprintf(buf, sizeof(buf), "%.17g", a->value.u.f);
        int  has = 0;
        for (int i = 0; i < n; i++)
            if (buf[i] == '.' || buf[i] == 'e' || buf[i] == 'E') { has = 1; break; }
        fputs(buf, out);
        if (!has) fputs(".0", out);
    } else if (a->value.kind == KFLV_IDENT) {
        fputs(a->value.u.s ? a->value.u.s : "?", out);
    } else {
        emit_value(out, &a->value);
    }
}

/* Emit attribute (`key value` form, or just the flag word for boolean
 * flags like `divider`/`toggle`). */
static void emit_attr_inline(FILE *out, const KflcAttr *a)
{
    int is_flag = ((strcmp(a->name, "divider") == 0 || strcmp(a->name, "toggle") == 0) &&
                   a->value.kind == KFLV_INT && a->value.u.i == 1);
    if (is_flag) {
        fputs(a->name, out);
        return;
    }
    /* Widget-context `at` and `size` parse as `<int> <int>` (two
     * positional ints packed into a KFLV_SIZE), not as the `WxH`
     * T_SIZE token used at form scope. Emit them back in the
     * two-int form so the round-trip parses cleanly. The form-scope
     * `size <WxH>` path runs through emit_form_attrs (separate code
     * path) and stays unaffected. */
    if (a->value.kind == KFLV_SIZE &&
        (strcmp(a->name, "at") == 0 || strcmp(a->name, "size") == 0)) {
        fprintf(out, "%s %d %d", a->name, a->value.u.size.w, a->value.u.size.h);
        return;
    }
    fputs(a->name, out);
    fputc(' ', out);
    emit_value(out, &a->value);
}

/* Look up an attr by exact name in a linked list. NULL on miss. */
static const KflcAttr *find_attr_(const KflcAttr *list, const char *name)
{
    for (const KflcAttr *a = list; a; a = a->next) {
        if (strcmp(a->name, name) == 0) return a;
    }
    return NULL;
}

static int in_skip_list_(const char *name, const char *const *skip)
{
    if (!skip) return 0;
    for (int i = 0; skip[i]; i++) {
        if (strcmp(name, skip[i]) == 0) return 1;
    }
    return 0;
}

static void emit_attrs_inline_skip(FILE *out, const KflcAttr *list,
                                    const char *const *skip)
{
    for (const KflcAttr *a = list; a; a = a->next) {
        if (in_skip_list_(a->name, skip)) continue;
        /* Skip _y / _z members of triple-attrs — the _x emit (below)
         * already wrote all three as `<base> X Y Z`. */
        if (triple_base_for_suffix_(a->name, 'y') ||
            triple_base_for_suffix_(a->name, 'z')) {
            continue;
        }
        const char *base = triple_base_for_suffix_(a->name, 'x');
        if (base) {
            char ny[64], nz[64];
            snprintf(ny, sizeof(ny), "%s_y", base);
            snprintf(nz, sizeof(nz), "%s_z", base);
            const KflcAttr *ay = find_attr_(list, ny);
            const KflcAttr *az = find_attr_(list, nz);
            if (ay && az) {
                fputc(' ', out);
                fputs(base, out);
                fputc(' ', out); emit_attr_value_raw_(out, a);
                fputc(' ', out); emit_attr_value_raw_(out, ay);
                fputc(' ', out); emit_attr_value_raw_(out, az);
                continue;
            }
            /* Fall through to plain emit if the y/z partners are
             * missing (unusual; means an upstream parser pass
             * produced an inconsistent triple). */
        }
        fputc(' ', out);
        emit_attr_inline(out, a);
    }
}

static void emit_attrs_inline(FILE *out, const KflcAttr *list)
{
    emit_attrs_inline_skip(out, list, NULL);
}

static void emit_node(FILE *out, const KflcNode *n, int level);

static void emit_form_attrs(FILE *out, const KflcNode *form, int level)
{
    /* Order: title, size, cfg. Then args / children. */
    static const char *singletons[] = { "title", "size", "cfg" };
    for (size_t i = 0; i < sizeof(singletons) / sizeof(singletons[0]); i++) {
        for (const KflcAttr *a = form->attrs; a; a = a->next) {
            if (strcmp(a->name, singletons[i]) == 0) {
                indent(out, level);
                fprintf(out, "%s  ", a->name);
                emit_value(out, &a->value);
                fputc('\n', out);
            }
        }
    }
}

static void emit_node(FILE *out, const KflcNode *n, int level)
{
    switch (n->kind) {
    case KFLN_FORM:
        fprintf(out, "form %s\n", n->name ? n->name : "FORM");
        emit_form_attrs(out, n, level + 1);
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        fputs("end\n", out);
        break;

    case KFLN_ARG: {
        indent(out, level);
        int dflt = (n->position.kind == KFLV_INT && n->position.u.i != 0);
        int ro   = (n->flags & KFL_NF_READONLY) ? 1 : 0;
        fprintf(out, "arg %s default %s%s\n",
                n->name ? n->name : "?",
                dflt ? "true" : "false",
                ro   ? " readonly" : "");
        break;
    }

    case KFLN_WIDGET:
        indent(out, level);
        fputs(widget_kw(n->widget), out);
        if (n->position.kind != KFLV_NONE) {
            fputc(' ', out);
            emit_value(out, &n->position);
        }
        emit_attrs_inline(out, n->attrs);
        fputc('\n', out);
        break;

    /* Function definitions and their body statements. Optional
     * lifetime qualifier on the return type and on each fn-arg,
     * prefixed before the type ident. */
    case KFLN_FN: {
        const char *rt = kflc_type_kfl_str(n->type, n->type_subtype);
        if (!rt) rt = "double";
        const char *ret_lqs = kflc_lifetime_qualifier_kfl_str(n->lifetime_qualifier);
        indent(out, level);
        if (ret_lqs) {
            fprintf(out, "fn %s %s %s(", ret_lqs, rt, n->name ? n->name : "?");
        } else {
            fprintf(out, "fn %s %s(", rt, n->name ? n->name : "?");
        }
        int aidx = 0;
        for (const KflcNode *a = n->children; a; a = a->next) {
            if (a->kind != KFLN_FN_ARG) continue;
            if (aidx > 0) fputs(", ", out);
            const char *at = kflc_type_kfl_str(a->type, a->type_subtype);
            const char *arg_lqs = kflc_lifetime_qualifier_kfl_str(a->lifetime_qualifier);
            if (arg_lqs) {
                fprintf(out, "%s %s %s",
                        arg_lqs, at ? at : "double",
                        a->name ? a->name : "?");
            } else {
                fprintf(out, "%s %s", at ? at : "double", a->name ? a->name : "?");
            }
            aidx++;
        }
        fputs(")\n", out);
        for (const KflcNode *c = n->children; c; c = c->next) {
            if (c->kind == KFLN_FN_ARG) continue;
            emit_node(out, c, level + 1);
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    }
    case KFLN_FN_ARG:
        /* Emitted inside the parent fn header. */
        break;
    case KFLN_STMT_LET:
    case KFLN_STMT_CONST: {
        indent(out, level);
        const char *ty = kflc_type_kfl_str(n->type, n->type_subtype);
        /* Emit optional lifetime qualifier between the colon
         * and the type, mirroring the parser's accept order. NULL
         * from kflc_lifetime_qualifier_kfl_str (KFL_LQ_NONE) skips
         * the prefix entirely so v3.x round-trip is byte-exact. */
        const char *lqs = kflc_lifetime_qualifier_kfl_str(n->lifetime_qualifier);
        if (lqs) {
            fprintf(out, "%s %s: %s %s = ",
                    n->kind == KFLN_STMT_CONST ? "const" : "let",
                    n->name ? n->name : "?",
                    lqs,
                    ty ? ty : "double");
        } else {
            fprintf(out, "%s %s: %s = ",
                    n->kind == KFLN_STMT_CONST ? "const" : "let",
                    n->name ? n->name : "?",
                    ty ? ty : "double");
        }
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("0", out);
        fputc('\n', out);
        break;
    }
    case KFLN_STMT_ASSIGN:
        indent(out, level);
        fprintf(out, "%s = ", n->name ? n->name : "?");
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("0", out);
        fputc('\n', out);
        break;
    case KFLN_STMT_RETURN:
        indent(out, level);
        if (n->expr) { fputs("return ", out); kflc_expr_to_text(out, n->expr); fputc('\n', out); }
        else         fputs("return\n", out);
        break;
    case KFLN_STMT_EXPR:
        indent(out, level);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        fputc('\n', out);
        break;
    case KFLN_STMT_PRINT: {
        indent(out, level);
        fputs("print", out);
        int first = 1;
        for (const KflcNode *arg = n->children; arg; arg = arg->next) {
            fputs(first ? " " : ", ", out);
            first = 0;
            if (arg->position.kind == KFLV_STR) {
                emit_string(out, arg->position.u.s);
            } else if (arg->expr) {
                kflc_expr_to_text(out, arg->expr);
            }
        }
        fputc('\n', out);
        break;
    }
    case KFLN_STMT_IF:
        indent(out, level);
        fputs("if ", out);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        fputc('\n', out);
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        if (n->else_children) {
            indent(out, level);
            fputs("else\n", out);
            for (const KflcNode *c = n->else_children; c; c = c->next) {
                emit_node(out, c, level + 1);
            }
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    case KFLN_STMT_WHILE:
        indent(out, level);
        fputs("while ", out);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        fputc('\n', out);
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    case KFLN_FN_DATA: {
        indent(out, level);
        fprintf(out, "fn data %s\n", n->name ? n->name : "?");
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    }
    case KFLN_STMT_INDEX_ASSIGN:
        indent(out, level);
        fprintf(out, "%s[", n->name ? n->name : "?");
        if (n->expr2) kflc_expr_to_text(out, n->expr2);
        fputs("] = ", out);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("0", out);
        fputc('\n', out);
        break;
    case KFLN_STMT_LVALUE_ASSIGN:
        /* B2: LHS is a full KflcExpr — re-emit via kflc_expr_to_text
         * so any future lvalue shape (struct fields, deeper chains)
         * round-trips without further serialize edits. */
        indent(out, level);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("?", out);
        fputs(" = ", out);
        if (n->expr2) kflc_expr_to_text(out, n->expr2);
        else          fputs("0", out);
        fputc('\n', out);
        break;
    case KFLN_STMT_SERIES: {
        const char *kn =
            (n->position.u.i == 0) ? "series_line" :
            (n->position.u.i == 1) ? "series_scatter" :
            (n->position.u.i == 2) ? "series_errorbar" :
            (n->position.u.i == 3) ? "series_histogram" :
            (n->position.u.i == 4) ? "series_box" :
                                     "series_heatmap";
        if (n->position.u.i == 5) {
            /* Heatmap: single matrix identifier. */
            const char *mat_name = "?";
            for (const KflcAttr *a = n->attrs; a; a = a->next) {
                if (strcmp(a->name, "mat") == 0 && a->value.kind == KFLV_IDENT)
                    mat_name = a->value.u.s;
            }
            indent(out, level);
            fprintf(out, "%s ", kn);
            emit_string(out, n->name ? n->name : "");
            fprintf(out, " %s\n", mat_name);
            break;
        }
        const char *xs_name = "?", *ys_name = "?";
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "xs") == 0 && a->value.kind == KFLV_IDENT)
                xs_name = a->value.u.s;
            if (strcmp(a->name, "ys") == 0 && a->value.kind == KFLV_IDENT)
                ys_name = a->value.u.s;
        }
        indent(out, level);
        fprintf(out, "%s ", kn);
        emit_string(out, n->name ? n->name : "");
        fprintf(out, " %s %s\n", xs_name, ys_name);
        break;
    }

    case KFLN_FN_WORLD: {
        indent(out, level);
        fprintf(out, "fn world %s\n", n->name ? n->name : "?");
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    }

    case KFLN_TICK: {
        indent(out, level);
        fprintf(out, "tick %s", n->name ? n->name : "?");
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            fputc(' ', out);
            fputs(a->name, out);
            fputc(' ', out);
            switch (a->value.kind) {
            case KFLV_INT:   fprintf(out, "%ld", a->value.u.i); break;
            case KFLV_FLOAT: fprintf(out, "%.17g", a->value.u.f); break;
            case KFLV_IDENT: fputs(a->value.u.s, out); break;
            case KFLV_STR:   fprintf(out, "\"%s\"", a->value.u.s); break;
            default: fputs("?", out); break;
            }
        }
        fputc('\n', out);
        break;
    }

    case KFLN_FRAME: {
        indent(out, level);
        fprintf(out, "frame %s", n->name ? n->name : "?");
        /* `kind` is the second positional in the source surface
         * (inertial / body_fixed); other attrs (e.g. `body`) trail. */
        const KflcAttr *kind_a = NULL;
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "kind") == 0) { kind_a = a; break; }
        }
        if (kind_a && kind_a->value.kind == KFLV_IDENT) {
            fprintf(out, " %s", kind_a->value.u.s);
        } else {
            fputs(" inertial", out);
        }
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (a == kind_a) continue;
            fputc(' ', out);
            fputs(a->name, out);
            fputc(' ', out);
            switch (a->value.kind) {
            case KFLV_IDENT: fputs(a->value.u.s, out); break;
            case KFLV_STR:   fprintf(out, "\"%s\"", a->value.u.s); break;
            case KFLV_INT:   fprintf(out, "%ld", a->value.u.i); break;
            default: fputs("?", out); break;
            }
        }
        fputc('\n', out);
        break;
    }

    case KFLN_EPOCH_LITERAL: {
        indent(out, level);
        const KflcAttr *iso_a = NULL, *sc_a = NULL;
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if      (strcmp(a->name, "iso")   == 0) iso_a = a;
            else if (strcmp(a->name, "scale") == 0) sc_a  = a;
        }
        fprintf(out, "epoch %s \"%s\" %s\n",
                n->name ? n->name : "?",
                (iso_a && iso_a->value.kind == KFLV_STR) ? iso_a->value.u.s : "?",
                (sc_a  && sc_a->value.kind  == KFLV_IDENT) ? sc_a->value.u.s : "UTC");
        break;
    }

    /* ---- Astro statements ------------------------------------- */

    case KFLN_STMT_ASTRO_BODY: {
        indent(out, level);
        fprintf(out, "astro_body %s", n->name ? n->name : "?");
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            const char *v = (a->value.kind == KFLV_IDENT && a->value.u.s)
                            ? a->value.u.s : "?";
            fprintf(out, " %s=%s", a->name, v);
        }
        fputc('\n', out);
        break;
    }

    case KFLN_STMT_STEP: {
        indent(out, level);
        fputs("step ", out);
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("?", out);
        fputc('\n', out);
        break;
    }

    case KFLN_STMT_PROPAGATE: {
        indent(out, level);
        fprintf(out, "propagate %s for ", n->name ? n->name : "?");
        if (n->expr) kflc_expr_to_text(out, n->expr);
        else         fputs("?", out);
        fputc('\n', out);
        break;
    }

    case KFLN_STMT_FOR_EACH: {
        const char *world = NULL;
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "world") == 0) world = a->value.u.s;
        }
        indent(out, level);
        fprintf(out, "for_each %s in %s\n",
                n->name ? n->name : "?", world ? world : "?");
        for (const KflcNode *c = n->children; c; c = c->next) {
            emit_node(out, c, level + 1);
        }
        indent(out, level);
        fputs("end\n", out);
        break;
    }

    case KFLN_STMT_OBSERVE: {
        const char *observer = NULL;
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "observer") == 0) observer = a->value.u.s;
        }
        indent(out, level);
        fprintf(out, "observe %s from %s",
                n->name ? n->name : "?",
                observer ? observer : "?");
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "observer") == 0) continue;
            const char *v = (a->value.kind == KFLV_IDENT && a->value.u.s)
                            ? a->value.u.s : "?";
            fprintf(out, " %s=%s", a->name, v);
        }
        fputc('\n', out);
        break;
    }

    /* KFL memory model round-trip surface:
     *   arena <name> capacity = <bytes>
     *   allocator = <arena_name>
     * Capacity is emitted as raw bytes; the parser accepts both raw
     * int and unit-suffixed forms (KB / MB / GB) but the serializer's
     * canonical form is unsuffixed so round-trip is exact for any
     * byte count. */
    case KFLN_ARENA: {
        long capacity = 0;
        const char *reset_mode = NULL;
        for (const KflcAttr *a = n->attrs; a; a = a->next) {
            if (strcmp(a->name, "capacity") == 0 && a->value.kind == KFLV_INT)
                capacity = a->value.u.i;
            else if (strcmp(a->name, "reset_mode") == 0 &&
                     a->value.kind == KFLV_IDENT && a->value.u.s)
                reset_mode = a->value.u.s;
        }
        indent(out, level);
        /* Positional KFL form-scope syntax: `arena <name> capacity
         * <bytes>`. Round-trip emits raw bytes so the parser doesn't
         * have to round-trip unit suffixes — any whole-byte capacity
         * survives the parse → serialize → parse fixed-point. */
        fprintf(out, "arena %s capacity %ld",
                n->name ? n->name : "?", capacity);
        if (reset_mode) fprintf(out, " reset_mode %s", reset_mode);
        fputc('\n', out);
        break;
    }

    case KFLN_ALLOCATOR_BIND: {
        indent(out, level);
        fprintf(out, "allocator = %s\n", n->name ? n->name : "?");
        break;
    }
    }
}

int kflc_serialize(FILE *out, const KflcNode *form)
{
    if (!out || !form || form->kind != KFLN_FORM) return 1;
    fputs("# kflc-generated KFL source.\n", out);
    emit_node(out, form, 0);
    return ferror(out) ? 1 : 0;
}
