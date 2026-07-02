/* kflc — C++ emitter.
 *
 * Walks the parsed KFL AST and emits a complete standalone .cc source
 * that builds into an artifact-producing (non-interactive) program:
 * compute + simulation results print to stdout, and plot widgets render
 * to <name>.png / <name>.svg.
 *
 * The emitted program's main():
 *   - Parses `arg <name>` declarations as `--<name>` / `--no-<name>`
 *     command-line overrides (readonly args keep their declared default).
 *   - Runs each `fn world` body against a freshly created astronomical
 *     world (PORTABLE mode for bit-identical reproducibility); the body's
 *     `observe` statements print to stdout.
 *   - Renders every plot widget via libk26plot.
 *   - Calls the conventional parameterless entry fn `run` if defined.
 *
 * Emitted TU structure: runtime includes (only those the program uses),
 * the B5 bounds-check macros, form-arg globals, form-level bump arenas,
 * user `fn` / `fn world` / `fn data` bodies, then main().
 */

#include "kflc.h"
#include "internal.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Attr lookup helpers ----------------------------------------- */

static const KflcAttr *attr_find(const KflcNode *n, const char *name)
{
    for (const KflcAttr *a = n->attrs; a; a = a->next) {
        if (strcmp(a->name, name) == 0) return a;
    }
    return NULL;
}

/* Count and enumerate args in declaration order. */
static int count_args(const KflcNode *form)
{
    int n = 0;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind == KFLN_ARG) n++;
    }
    return n;
}

/* Count KFLN_FN children at form scope (user-declared functions). */
static int count_user_fns(const KflcNode *form)
{
    int n = 0;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind == KFLN_FN) n++;
    }
    return n;
}

/* Look up a KFLN_FN_DATA by name at form scope. Used to decide
 * whether a `plot data <name>` resolves to an inline KFL data
 * provider (in which case no extern decl is needed) or an external
 * C function in the handlers .cc. */
static const KflcNode *find_fn_data(const KflcNode *form, const char *name)
{
    if (!name) return NULL;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind == KFLN_FN_DATA && c->name && strcmp(c->name, name) == 0) {
            return c;
        }
    }
    return NULL;
}

/* Count fn-arg children of a KFLN_FN node. */
static int count_fn_args(const KflcNode *fn)
{
    int n = 0;
    for (const KflcNode *c = fn->children; c; c = c->next) {
        if (c->kind == KFLN_FN_ARG) n++;
    }
    return n;
}

/* Walk a fn-body subtree and append every let / const binding
 * (any nesting depth) into *live. The expression emitter checks
 * `is_bound` to decide whether KFLE_IDENT resolves, so a binding
 * declared inside an `if` / `while` must appear in this list or
 * any reference to it in the loop body will be rejected. We
 * accept the slight semantic looseness — KFL won't catch
 * "use before declare" itself; C++ does, at the downstream compile. */
/* Pre-populate the fn-body binding table with the form's
 * KFLN_ARG declarations so expressions inside `fn` / `fn data` /
 * `fn scene` can reference `kfl_arg_<name>` by its KFL name. String
 * args get type KFLT_STRING; the legacy bool args are exposed as
 * KFLT_BOOL. Vector / matrix arg types are not currently supported. */
static void collect_form_args(const KflcNode *form,
                              KflcArena *arena,
                              KflcExprBinding **live,
                              int *live_n, int *live_cap)
{
    if (!form) return;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind != KFLN_ARG || !c->name) continue;
        if (*live_n == *live_cap) {
            int nc = (*live_cap ? *live_cap * 2 : 4);
            KflcExprBinding *nl = (KflcExprBinding *)kflc_arena_alloc(
                arena, sizeof(KflcExprBinding) * (size_t)nc);
            if (*live_n > 0) {
                memcpy(nl, *live, sizeof(KflcExprBinding) * (size_t)(*live_n));
            }
            *live = nl;
            *live_cap = nc;
        }
        /* Zero-init the slot so KflcExprBinding fields
         * (type_subtype, lifetime_qualifier) default cleanly even
         * when not set per-arg below. The arena is malloc-backed,
         * not calloc, so this matters. */
        memset(&(*live)[*live_n], 0, sizeof(KflcExprBinding));
        (*live)[*live_n].name = c->name;
        /* Type heuristic for KFLN_ARG. String defaults give
         * KFLT_STRING. Float defaults → KFLT_DOUBLE. Int defaults
         * that aren't 0/1 → KFLT_INT (negative or > 1). Else KFLT_BOOL
         * for the legacy true/false pattern. */
        if (c->position.kind == KFLV_STR) {
            (*live)[*live_n].type = KFLT_STRING;
        } else if (c->position.kind == KFLV_FLOAT) {
            (*live)[*live_n].type = KFLT_DOUBLE;
        } else if (c->position.kind == KFLV_INT &&
                   (c->position.u.i < 0 || c->position.u.i > 1)) {
            (*live)[*live_n].type = KFLT_INT;
        } else {
            (*live)[*live_n].type = KFLT_BOOL;
        }
        (*live)[*live_n].is_form_arg      = 1;
        (*live)[*live_n].is_observed_cell =
            (c->flags & KFL_NF_OBSERVED_CELL) ? 1 : 0;
        /* Form args carry no lifetime qualifier (scalar/string args
         * have no ownership semantics). Explicit assignment for
         * clarity; the zero-init above already leaves this at
         * KFL_LQ_NONE. */
        (*live)[*live_n].lifetime_qualifier = KFL_LQ_NONE;
        /* Ownership tracking defaults; form args are neither moved
         * nor borrow-sourced from anywhere. */
        (*live)[*live_n].moved_from        = 0;
        (*live)[*live_n].borrow_source_idx = -1;
        (*live)[*live_n].scope_depth       = 0;
        (*live_n)++;
    }
}

/* Count statements of `kind` in the tree rooted at `n`, recursing
 * through children AND else_children so statements nested inside
 * if/while bodies are tallied. Used by the fn_scene / fn_data
 * preamble emitters to size their static caches to the actual node
 * count rather than the historical 256/16 ceilings. Cheap — one walk
 * over a small AST at codegen time, no allocation. */
static int count_stmt_kind(const KflcNode *n, KflcNodeKind kind)
{
    if (!n) return 0;
    int total = (n->kind == kind) ? 1 : 0;
    for (const KflcNode *c = n->children; c; c = c->next) {
        total += count_stmt_kind(c, kind);
    }
    for (const KflcNode *c = n->else_children; c; c = c->next) {
        total += count_stmt_kind(c, kind);
    }
    return total;
}

/* Sum count_stmt_kind across a body's siblings. Use when the body
 * head is the first stmt of an fn (siblings linked via .next) rather
 * than a single subtree root. */
static int count_stmt_kind_body(const KflcNode *body_first, KflcNodeKind kind)
{
    int total = 0;
    for (const KflcNode *s = body_first; s; s = s->next) {
        total += count_stmt_kind(s, kind);
    }
    return total;
}

static void collect_let_bindings(const KflcNode *n,
                                  KflcArena *arena,
                                  KflcExprBinding **live,
                                  int *live_n, int *live_cap)
{
    if (!n) return;
    if (n->kind == KFLN_STMT_LET || n->kind == KFLN_STMT_CONST) {
        if (*live_n == *live_cap) {
            int nc = (*live_cap ? *live_cap * 2 : 4);
            KflcExprBinding *nl = (KflcExprBinding *)kflc_arena_alloc(
                arena, sizeof(KflcExprBinding) * (size_t)nc);
            if (*live_n > 0) {
                memcpy(nl, *live, sizeof(KflcExprBinding) * (size_t)(*live_n));
            }
            *live = nl;
            *live_cap = nc;
        }
        /* Zero-init for defence-in-depth, then propagate the full
         * binding surface (type, subtype, lifetime) from the
         * let/const node. */
        memset(&(*live)[*live_n], 0, sizeof(KflcExprBinding));
        (*live)[*live_n].name         = n->name;
        (*live)[*live_n].type         = n->type;
        (*live)[*live_n].type_subtype = n->type_subtype;
        (*live)[*live_n].is_form_arg  = 0;
        (*live)[*live_n].is_observed_cell = 0;
        /* Propagate lifetime qualifier from the let/const declaration
         * so downstream emit can see own/borrow/ptr through
         * expression contexts. */
        (*live)[*live_n].lifetime_qualifier = n->lifetime_qualifier;
        /* Ownership tracking defaults: moved_from defaults to 0 (not
         * moved); borrow_source_idx defaults to -1 (no source
         * resolved yet; emit-time identifier resolver fills it for
         * borrow bindings). scope_depth left at 0; stmt.c records
         * the actual block depth on the binding at scope-add time
         * once that hook lands. */
        (*live)[*live_n].moved_from        = 0;
        (*live)[*live_n].borrow_source_idx = -1;
        (*live)[*live_n].scope_depth       = 0;
        (*live_n)++;
    }
    for (const KflcNode *c = n->children; c; c = c->next) {
        collect_let_bindings(c, arena, live, live_n, live_cap);
    }
    if (n->kind == KFLN_STMT_IF) {
        for (const KflcNode *c = n->else_children; c; c = c->next) {
            collect_let_bindings(c, arena, live, live_n, live_cap);
        }
    }
}

/* ---- Token translation ------------------------------------------- */

static void emit_string_literal(FILE *out, const char *s)
{
    fputc('"', out);
    for (const unsigned char *p = (const unsigned char *)s; p && *p; p++) {
        unsigned char c = *p;
        switch (c) {
        case '"':  fputs("\\\"", out); break;
        case '\\': fputs("\\\\", out); break;
        case '\n': fputs("\\n",  out); break;
        case '\t': fputs("\\t",  out); break;
        case '\r': fputs("\\r",  out); break;
        default:
            if (c < 0x20) fprintf(out, "\\x%02x", c);
            else          fputc(c, out);
            break;
        }
    }
    fputc('"', out);
}

/* ---- Handler / painter collection -------------------------------- */

typedef struct {
    char       **names;
    int          count;
    int          cap;
    KflcArena   *arena;
} NameSet;

static void nset_init(NameSet *s, KflcArena *arena)
{
    s->names = NULL;
    s->count = 0;
    s->cap   = 0;
    s->arena = arena;
}

static void nset_add(NameSet *s, const char *name)
{
    if (!name) return;
    for (int i = 0; i < s->count; i++) {
        if (strcmp(s->names[i], name) == 0) return;
    }
    if (s->count == s->cap) {
        int newcap = s->cap ? s->cap * 2 : 8;
        char **fresh = (char **)kflc_arena_alloc(s->arena,
                                                  sizeof(char *) * (size_t)newcap);
        if (s->count > 0) memcpy(fresh, s->names, sizeof(char *) * (size_t)s->count);
        s->names = fresh;
        s->cap   = newcap;
    }
    s->names[s->count++] = (char *)name;
}

/* Collect the `data <ident>` provider names referenced by every plot
 * widget in the form; each names a series-providing C function with
 * the libk26plot data-callback signature (K26PlotConfig + K26PSeries
 * out-pointers), which may also be satisfied by an inline `fn data`. */
static void collect_handlers(const KflcNode *n, NameSet *plot_data)
{
    if (!n) return;

    if (n->kind == KFLN_WIDGET && n->widget == KFL_W_PLOT) {
        const KflcAttr *a = attr_find(n, "data");
        if (a && a->value.kind == KFLV_IDENT) nset_add(plot_data, a->value.u.s);
    }
    for (const KflcNode *c = n->children; c; c = c->next) {
        collect_handlers(c, plot_data);
    }
}

/* True if the form contains at least one plot widget. */
static int form_has_plot(const KflcNode *n)
{
    if (!n) return 0;
    if (n->kind == KFLN_WIDGET && n->widget == KFL_W_PLOT) return 1;
    for (const KflcNode *c = n->children; c; c = c->next) {
        if (form_has_plot(c)) return 1;
    }
    return 0;
}

/* Recursive scan for any memory-model surface in the AST. Triggers
 * the conditional `#include "kflc_runtime.h"` in the emit prologue
 * and gates the form-level arena declaration/init/destroy sweep. */
static int node_uses_w8_1_(const KflcNode *n)
{
    if (!n) return 0;
    if (n->kind == KFLN_ARENA || n->kind == KFLN_ALLOCATOR_BIND) return 1;
    if (n->lifetime_qualifier != KFL_LQ_NONE) return 1;
    for (const KflcNode *c = n->children; c; c = c->next) {
        if (node_uses_w8_1_(c)) return 1;
    }
    for (const KflcNode *c = n->else_children; c; c = c->next) {
        if (node_uses_w8_1_(c)) return 1;
    }
    return 0;
}

static int form_uses_w8_1(const KflcNode *form)
{
    if (!form) return 0;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (node_uses_w8_1_(c)) return 1;
    }
    return 0;
}

static int form_has_arena(const KflcNode *form)
{
    if (!form) return 0;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind == KFLN_ARENA) return 1;
    }
    return 0;
}

/* ---- Top-level emit ---------------------------------------------- */

/* Headless plot rendering. Walks the form tree for plot widgets and,
 * for each `plot <name> data <fn>`, emits a block that runs the data
 * function and writes <name>.png + <name>.svg via k26plot_render, using
 * the widget's title / axis labels for the figure config. Recursion
 * mirrors form_has_plot (children + KFLN_IF else_children). */
static void emit_headless_plots_(FILE *out, const KflcNode *n, int *counter)
{
    if (!n) return;
    if (n->kind == KFLN_WIDGET && n->widget == KFL_W_PLOT) {
        const KflcAttr *data_a  = attr_find(n, "data");
        const KflcAttr *title_a = attr_find(n, "title");
        const KflcAttr *xlab_a  = attr_find(n, "x_label");
        const KflcAttr *ylab_a  = attr_find(n, "y_label");
        const char *data_fn = (data_a && data_a->value.kind == KFLV_IDENT)
                              ? data_a->value.u.s : NULL;
        if (data_fn) {
            char base[128];
            if (n->name) snprintf(base, sizeof base, "%s", n->name);
            else         snprintf(base, sizeof base, "plot_%d", *counter);
            (*counter)++;
            fputs("    {\n", out);
            fputs("        K26PlotConfig _pc; k26plot_config_defaults(&_pc);\n", out);
            if (title_a && title_a->value.kind == KFLV_STR) {
                fputs("        _pc.title = ", out);
                emit_string_literal(out, title_a->value.u.s);
                fputs(";\n", out);
            }
            if (xlab_a && xlab_a->value.kind == KFLV_STR) {
                fputs("        _pc.xlabel = ", out);
                emit_string_literal(out, xlab_a->value.u.s);
                fputs(";\n", out);
            }
            if (ylab_a && ylab_a->value.kind == KFLV_STR) {
                fputs("        _pc.ylabel = ", out);
                emit_string_literal(out, ylab_a->value.u.s);
                fputs(";\n", out);
            }
            fputs("        const K26PSeries *_ps = NULL; size_t _pn = 0;\n", out);
            fprintf(out, "        %s(&_pc, &_ps, &_pn, nullptr);\n", data_fn);
            fprintf(out,
                "        K26PStatus _pst = k26plot_render(&_pc, _ps, _pn, "
                "\"%s.png\", \"%s.svg\");\n", base, base);
            fprintf(out,
                "        if (_pst == K26P_OK) fprintf(stdout, "
                "\"wrote %s.png and %s.svg (%%zu series)\\n\", _pn);\n",
                base, base);
            fprintf(out,
                "        else fprintf(stderr, \"kflc: plot render failed "
                "for %s (status %%d)\\n\", (int)_pst);\n", base);
            fputs("    }\n", out);
        }
    }
    for (const KflcNode *c = n->children; c; c = c->next)
        emit_headless_plots_(out, c, counter);
}

int kflc_emit_cxx(FILE *out, const KflcNode *form, KflcDiag *diag)
{
    if (!form || form->kind != KFLN_FORM) {
        kflc_diag_errorf(diag, 0, "emit: not a form node");
        return 1;
    }
    /* A form needs no window header (title / size / cfg) — artifact
     * programs have no window. Any such attrs are ignored if present. */
    const char *form_name = form->name ? form->name : "FORM";

    /* Collect args in declaration order. */
    int nargs = count_args(form);

    /* Collect referenced plot-data provider names. */
    KflcArena *tmp_arena = kflc_arena_create();
    NameSet plot_data;
    nset_init(&plot_data, tmp_arena);
    collect_handlers(form, &plot_data);

    /* Banner + includes. */
    fprintf(out,
        "/* AUTO-GENERATED by kflc. DO NOT EDIT BY HAND.\n"
        " * Source form: %s\n"
        " */\n\n", form_name);

    /* NB: fputs (not fprintf) because the B5 KFLC_BOUNDS_CHECK macro
     * body below contains literal `%zu` / `%s` / `%d` printf-format
     * tokens (they're inside the emitted std::fprintf macro body,
     * meant for the user's runtime printf, NOT for this emit-time
     * call). Passing them through fprintf without matching args is
     * UB and segfaults on musl (silently consumed by glibc). */
    /* Pull in only the runtime headers the program actually uses — the
     * astronomical runtime for `fn world`, the plotting library for
     * `fn data` / plot widgets, and the compute library for vector
     * maths. */
    int hl_world = 0;
    for (const KflcNode *c = form->children; c; c = c->next) {
        if (c->kind == KFLN_FN_WORLD) { hl_world = 1; break; }
    }
    int hl_plot = form_has_plot(form);
    if (hl_world) {
        fputs(
            "#include <k26astro_rt/world.h>\n"
            "#include <k26astro_rt/observer.h>\n"
            "#include <k26astro_body/body.h>\n"
            "#include <k26astro_core/pos.h>\n",
            out);
    }
    if (hl_plot) {
        fputs("#include \"k26plot.h\"\n", out);
    }
    if (hl_plot) {
        fputs("#include \"k26compute.h\"\n", out);
    }
    if (hl_world || hl_plot) {
        fputs("#include <k26m3d.h>\n", out);
    }
    fputc('\n', out);
    /* A pure-compute program (no world / plot / compute vectors)
     * needs only the C standard headers emitted below. */
    fputs(
        "#include <stdint.h>\n"
        "#include <stdio.h>\n"
        "#include <string.h>\n"
        "#include <cmath>\n"
        "\n",
        out);

    /* When the form uses any memory-model keyword (own / borrow /
     * ptr / arena / allocator), pull in the runtime header so
     * generated calls to k26kfl_arena_alloc / k26kfl_move compile.
     * Conditional so unannotated forms keep the same trim include
     * surface. */
    if (form_uses_w8_1(form)) {
        fputs("#include \"kflc_runtime.h\"\n\n", out);
    }

    fputs(
        "/* B5: vector-element access. Release-build path is the bare\n"
        " * `.data[i]` deref — zero overhead, same instruction sequence\n"
        " * kflc has always emitted. Compile with -DKFLC_BOUNDS_CHECK\n"
        " * to get an index < length guard that aborts on OOB with a\n"
        " * `<file>:<line>` diagnostic. Use during form-author dev or\n"
        " * when chasing a crash in an unfamiliar .kfl. The lambda form\n"
        " * is the lvalue-friendly pattern that keeps i_expr from being\n"
        " * evaluated twice (matters for any future side-effecting\n"
        " * index, e.g. `xs[next_i()] = ...`). */\n"
        "#ifdef KFLC_BOUNDS_CHECK\n"
        "#  include <cstdio>\n"
        "#  include <cstdlib>\n"
        "#  define KFLC_VEC_R(v, i_expr) ([&]() -> double {\\\n"
        "      size_t _kfl_i = (size_t)((int)(i_expr));\\\n"
        "      if (_kfl_i >= (v).n) {\\\n"
        "          std::fprintf(stderr, \"kfl: bounds R %zu >= n %zu \"\\\n"
        "                       \"at %s:%d\\n\",\\\n"
        "                       _kfl_i, (v).n, __FILE__, __LINE__);\\\n"
        "          std::abort();\\\n"
        "      }\\\n"
        "      return (v).data[_kfl_i];\\\n"
        "  }())\n"
        "#  define KFLC_VEC_W(v, i_expr) ([&]() -> double& {\\\n"
        "      size_t _kfl_i = (size_t)((int)(i_expr));\\\n"
        "      if (_kfl_i >= (v).n) {\\\n"
        "          std::fprintf(stderr, \"kfl: bounds W %zu >= n %zu \"\\\n"
        "                       \"at %s:%d\\n\",\\\n"
        "                       _kfl_i, (v).n, __FILE__, __LINE__);\\\n"
        "          std::abort();\\\n"
        "      }\\\n"
        "      return (v).data[_kfl_i];\\\n"
        "  }())\n"
        "#else\n"
        "#  define KFLC_VEC_R(v, i_expr) ((v).data[(size_t)((int)(i_expr))])\n"
        "#  define KFLC_VEC_W(v, i_expr) ((v).data[(size_t)((int)(i_expr))])\n"
        "#endif\n"
        "\n"
        "/* Matrix-element access. Mirrors the KFLC_VEC_R/W shape but\n"
        " * with two compares (row vs rows, col vs cols). Release path\n"
        " * is the bare row-major deref; bounds path aborts on OOB with\n"
        " * `<file>:<line>` + the offending indices and the matrix\n"
        " * shape. */\n"
        "#ifdef KFLC_BOUNDS_CHECK\n"
        "#  define KFLC_MAT_R(m, r_expr, c_expr) ([&]() -> double {\\\n"
        "      size_t _kfl_r = (size_t)((int)(r_expr));\\\n"
        "      size_t _kfl_c = (size_t)((int)(c_expr));\\\n"
        "      if (_kfl_r >= (m).rows || _kfl_c >= (m).cols) {\\\n"
        "          std::fprintf(stderr, \"kfl: bounds R [%zu][%zu] vs shape \"\\\n"
        "                       \"[%zu][%zu] at %s:%d\\n\",\\\n"
        "                       _kfl_r, _kfl_c, (m).rows, (m).cols,\\\n"
        "                       __FILE__, __LINE__);\\\n"
        "          std::abort();\\\n"
        "      }\\\n"
        "      return (m).data[_kfl_r * (m).cols + _kfl_c];\\\n"
        "  }())\n"
        "#  define KFLC_MAT_W(m, r_expr, c_expr) ([&]() -> double& {\\\n"
        "      size_t _kfl_r = (size_t)((int)(r_expr));\\\n"
        "      size_t _kfl_c = (size_t)((int)(c_expr));\\\n"
        "      if (_kfl_r >= (m).rows || _kfl_c >= (m).cols) {\\\n"
        "          std::fprintf(stderr, \"kfl: bounds W [%zu][%zu] vs shape \"\\\n"
        "                       \"[%zu][%zu] at %s:%d\\n\",\\\n"
        "                       _kfl_r, _kfl_c, (m).rows, (m).cols,\\\n"
        "                       __FILE__, __LINE__);\\\n"
        "          std::abort();\\\n"
        "      }\\\n"
        "      return (m).data[_kfl_r * (m).cols + _kfl_c];\\\n"
        "  }())\n"
        "#else\n"
        "#  define KFLC_MAT_R(m, r_expr, c_expr) \\\n"
        "      ((m).data[(size_t)((int)(r_expr)) * (m).cols + (size_t)((int)(c_expr))])\n"
        "#  define KFLC_MAT_W(m, r_expr, c_expr) \\\n"
        "      ((m).data[(size_t)((int)(r_expr)) * (m).cols + (size_t)((int)(c_expr))])\n"
        "#endif\n"
        "\n", out);

    /* Form arguments, emitted as extern "C" globals so downstream
     * translation units can read them via their own
     * `extern "C" bool kfl_arg_<name>;` declarations.
     *
     * String args back onto a mutable static buffer. The
     * const-char* alias keeps existing readers unchanged. */
    if (nargs > 0) {
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_ARG) continue;
            if (c->position.kind == KFLV_STR) {
                /* C1: readonly string args back onto a const buffer
                 * so the compiler enforces immutability of the
                 * contents too, not just the alias pointer. */
                int ro = (c->flags & KFL_NF_READONLY) ? 1 : 0;
                fprintf(out,
                    "static %schar kfl_arg_%s_buf[2048] = ",
                    ro ? "const " : "", c->name);
                emit_string_literal(out, c->position.u.s ? c->position.u.s : "");
                fputs(";\n", out);
            }
        }
        fputs("extern \"C\" {\n", out);
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_ARG) continue;
            /* C1: const-qualify readonly args so any C++ assignment
             * fails at compile time. For strings we also const the
             * pointer alias itself (`const char * const`) so the
             * pointer can't be re-pointed even though the contents
             * are already const. */
            int ro = (c->flags & KFL_NF_READONLY) ? 1 : 0;
            if (c->position.kind == KFLV_STR) {
                fprintf(out,
                    "    const char *%skfl_arg_%s = kfl_arg_%s_buf;\n",
                    ro ? "const " : "", c->name, c->name);
            } else if (c->position.kind == KFLV_FLOAT) {
                /* Double form-arg from a float default. */
                fprintf(out, "    %sdouble kfl_arg_%s = %.17g;\n",
                        ro ? "const " : "", c->name, c->position.u.f);
            } else if (c->position.kind == KFLV_INT &&
                       (c->position.u.i < 0 || c->position.u.i > 1)) {
                /* Int form-arg: KFLV_INT with value outside {0, 1}
                 * disambiguates from the legacy bool default (which
                 * uses true/false / 1/0). Negatives are a common
                 * sentinel ("no selection" = -1). */
                fprintf(out, "    %sint kfl_arg_%s = %ld;\n",
                        ro ? "const " : "", c->name, c->position.u.i);
            } else {
                int dflt = (c->position.kind == KFLV_INT) ? (int)c->position.u.i : 0;
                fprintf(out, "    %sbool kfl_arg_%s = %s;\n",
                        ro ? "const " : "", c->name, dflt ? "true" : "false");
            }
        }
        fputs("}\n\n", out);

    }

    /* Form-level bump arenas: static handle per `arena <name>` decl
     * plus a one-shot init. Form ctor calls `_kfl_init_arenas()`
     * after textures. Process exit relies on the OS to reap the
     * malloc'd buffers (no atexit destruction for simplicity; the
     * embedding process model is single-form-lives-for-the-process).
     * Future hot-reload use-cases may add explicit destroy. */
    if (form_has_arena(form)) {
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_ARENA || !c->name) continue;
            fprintf(out, "static K26KflArena *_kfl_arena_%s = NULL;\n",
                    c->name);
        }
        /* Cleanup helper for `reset_mode fn_exit` (the default). Used
         * by KFLN_ALLOCATOR_BIND emit (stmt.c) via
         * __attribute__((cleanup(...))) so the arena resets on every
         * fn return path automatically. GCC/Clang feature; the target
         * toolchain is gcc-everywhere so portability is fine. Manual
         * opt-out (`reset_mode manual`) skips the cleanup attribute
         * and leaves the arena lifetime to the user. */
        fputs("static void _kfl_arena_reset_cleanup_(K26KflArena **a) {\n",
              out);
        fputs("    if (a && *a) k26kfl_arena_reset(*a);\n", out);
        fputs("}\n", out);
        fputs("static void _kfl_init_arenas(void) {\n", out);
        fputs("    static int _ready = 0;\n", out);
        fputs("    if (_ready) return;\n", out);
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_ARENA || !c->name) continue;
            long capacity = 0;
            for (const KflcAttr *a = c->attrs; a; a = a->next) {
                if (strcmp(a->name, "capacity") == 0 &&
                    a->value.kind == KFLV_INT) {
                    capacity = a->value.u.i;
                    break;
                }
            }
            fprintf(out,
                "    _kfl_arena_%s = k26kfl_arena_create((size_t)%ldL);\n",
                c->name, capacity);
        }
        fputs("    _ready = 1;\n", out);
        fputs("}\n\n", out);
    }

    /* Extern decls. */
    fputs("extern \"C\" {\n", out);
    for (int i = 0; i < plot_data.count; i++) {
        /* Skip externs for plot-data callbacks resolved by an inline
         * `fn data <name>`; they're defined in this same TU just
         * below. */
        if (find_fn_data(form, plot_data.names[i])) continue;
        fprintf(out,
            "    void %s(K26PlotConfig *cfg, const K26PSeries **series,\n"
            "            size_t *n_series, void *user_data);\n",
            plot_data.names[i]);
    }
    fputs("}\n\n", out);

    /* Build user-fn table for expression resolution. */
    int n_user_fns = count_user_fns(form);
    KflcExprFn *user_fn_arr = NULL;
    if (n_user_fns > 0) {
        user_fn_arr = (KflcExprFn *)kflc_arena_alloc(tmp_arena,
            sizeof(KflcExprFn) * (size_t)n_user_fns);
        int idx = 0;
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_FN) continue;
            user_fn_arr[idx].name  = c->name;
            user_fn_arr[idx].arity = count_fn_args(c);
            idx++;
        }
    }
    KflcExprCtx form_expr_ctx;
    memset(&form_expr_ctx, 0, sizeof(form_expr_ctx));
    form_expr_ctx.bindings   = NULL;
    form_expr_ctx.n_bindings = 0;
    form_expr_ctx.fns        = user_fn_arr;
    form_expr_ctx.n_fns      = n_user_fns;
    form_expr_ctx.form       = form;   /* arena/borrow lookup */
    form_expr_ctx.headless   = 1;      /* artifact-only emitter */

    /* Emit each user fn as a TU-scope C++ function. */
    for (const KflcNode *fn = form->children; fn; fn = fn->next) {
        if (fn->kind != KFLN_FN) continue;
        const char *ret_cxx = kflc_type_cxx(fn->type, fn->type_subtype);
        if (!ret_cxx) ret_cxx = "double";
        fputs("extern \"C\" ", out);
        fprintf(out, "%s %s(", ret_cxx, fn->name ? fn->name : "_anon");

        /* Args. */
        int aidx = 0;
        for (const KflcNode *a = fn->children; a; a = a->next) {
            if (a->kind != KFLN_FN_ARG) continue;
            if (aidx > 0) fputs(", ", out);
            const char *aty = kflc_type_cxx(a->type, a->type_subtype);
            fprintf(out, "%s %s", aty ? aty : "double", a->name ? a->name : "_a");
            aidx++;
        }
        if (aidx == 0) fputs("void", out);
        fputs(") {\n", out);

        /* Build expression context: fn args as typed bindings, plus
         * form-level fn names, plus every let / const name declared
         * anywhere in the body (recursive — see collect_let_bindings).
         * Vector/matrix locals also accumulate into `pending_frees`
         * so we can emit k26c_vec_free / k26c_mat_free before each
         * return + at the natural end of the body. */
        int n_args = count_fn_args(fn);
        KflcExprBinding *live = NULL;
        int live_cap = n_args + 8;
        int live_n   = 0;
        live = (KflcExprBinding *)kflc_arena_alloc(tmp_arena,
            sizeof(KflcExprBinding) * (size_t)live_cap);
        for (const KflcNode *a = fn->children; a; a = a->next) {
            if (a->kind != KFLN_FN_ARG) continue;
            /* Zero-init then populate from fn_arg. Fn args are
             * caller-provided; borrow_source_idx stays -1 (no in-fn
             * source binding to track). is_form_arg stays 0; these
             * are FN args, not FORM args. */
            memset(&live[live_n], 0, sizeof(KflcExprBinding));
            live[live_n].name               = a->name;
            live[live_n].type               = a->type;
            live[live_n].type_subtype       = a->type_subtype;
            live[live_n].lifetime_qualifier = a->lifetime_qualifier;
            live[live_n].moved_from         = 0;
            live[live_n].borrow_source_idx  = -1;
            live[live_n].scope_depth        = 0;
            live_n++;
        }
        for (const KflcNode *s = fn->children; s; s = s->next) {
            if (s->kind == KFLN_FN_ARG) continue;
            collect_let_bindings(s, tmp_arena, &live, &live_n, &live_cap);
        }

        KflcExprCtx body_ctx;
        memset(&body_ctx, 0, sizeof(body_ctx));
        body_ctx.bindings   = live;
        body_ctx.n_bindings = live_n;
        body_ctx.fns        = user_fn_arr;
        body_ctx.n_fns      = n_user_fns;
        body_ctx.form       = form;   /* arena/borrow lookup */

        /* stmt.c's block-scope tracker owns every heap-typed let,
         * depth 0 (fn-body root) included. Its KFLN_STMT_RETURN case
         * drains all scopes and applies the `_kfl_rv` temp-save when
         * needed (knowing the return type via this reset call). At
         * natural fn-body end we still need to release any depth-0
         * lets that didn't escape via a return: kfl_emit_stmt_drain_root
         * does that. */
        kfl_emit_stmt_reset_scopes(tmp_arena, fn->type);

        int last_was_return = 0;
        for (const KflcNode *s = fn->children; s; s = s->next) {
            if (s->kind == KFLN_FN_ARG) continue;
            last_was_return = (s->kind == KFLN_STMT_RETURN);
            if (kfl_emit_stmt(out, s, &body_ctx, diag, 4)) {
                /* emit error already raised */
            }
        }
        if (!last_was_return) {
            kfl_emit_stmt_drain_root(out, 4);
        }
        fputs("}\n\n", out);
    }

    /* Emit `fn data <name>` callbacks. Signature matches the
     * libk26plot data callback (KflcPlot::DataFn). Body uses the same
     * statement palette as a regular fn plus the `series_<kind>`
     * statement, which pokes its out-parameters via the local names
     * `_kfl_out_series` and `_kfl_n_out` referenced by the emitter. */
    for (const KflcNode *fn = form->children; fn; fn = fn->next) {
        if (fn->kind != KFLN_FN_DATA) continue;
        /* Per-fn-data static array of series + their caches. Each
         * series_<kind> statement in the body slots its K26PSeries
         * into _kfl_all[_kfl_count++] using a runtime counter; at
         * function exit we hand _kfl_all + _kfl_count to libk26plot
         * via the out-parameters.
         *
         * Size the array from the actual statement count in the fn
         * body. Headroom of +4 covers conditional-branch series that
         * the AST count already captures, plus a small safety margin.
         * Floor of 4 keeps the array allocation legal when a fn
         * declares no series (compiles to a no-op callback). */
        int n_series_stmts = count_stmt_kind_body(fn->children,
                                                   KFLN_STMT_SERIES);
        int max_series = n_series_stmts + 4;
        if (max_series < 4) max_series = 4;
        fprintf(out,
            "extern \"C\" void %s(K26PlotConfig *cfg,\n"
            "                     const K26PSeries **_kfl_out_series,\n"
            "                     size_t *_kfl_n_out,\n"
            "                     void *_kfl_user_data) {\n"
            "    (void)cfg; (void)_kfl_user_data;\n"
            "    enum { _KFL_MAX_SERIES = %d };\n"
            "    static K26PSeries _kfl_all[_KFL_MAX_SERIES] = {};\n"
            "    static K26CVector _kfl_xh[_KFL_MAX_SERIES]  = {};\n"
            "    static K26CVector _kfl_yh[_KFL_MAX_SERIES]  = {};\n"
            "    static K26CVector _kfl_hh[_KFL_MAX_SERIES]  = {};\n"
            "    int _kfl_count = 0;\n",
            fn->name, max_series);

        /* No fn args for data callbacks; the C-ABI is fixed. Bindings
         * come entirely from let/const statements anywhere in the body
         * (recursive collect, as for regular fns). */
        KflcExprBinding *live = NULL;
        int live_cap = 8;
        int live_n   = 0;
        live = (KflcExprBinding *)kflc_arena_alloc(tmp_arena,
            sizeof(KflcExprBinding) * (size_t)live_cap);
        for (const KflcNode *s = fn->children; s; s = s->next) {
            collect_let_bindings(s, tmp_arena, &live, &live_n, &live_cap);
        }
        /* Form-level KFLN_ARGs are accessible inside any fn body as
         * KFLT_STRING / KFLT_BOOL identifiers. Appended after the let
         * bindings so a `let foo: int = ...` inside the body shadows
         * a like-named `arg foo` declaration (first-match lookup
         * direction). */
        collect_form_args(form, tmp_arena, &live, &live_n, &live_cap);

        KflcExprCtx body_ctx;
        memset(&body_ctx, 0, sizeof(body_ctx));
        body_ctx.bindings   = live;
        body_ctx.n_bindings = live_n;
        body_ctx.fns        = user_fn_arr;
        body_ctx.n_fns      = n_user_fns;
        body_ctx.form       = form;   /* arena/borrow lookup */

        /* Fn-data callbacks are void-returning, so the RETURN
         * temp-save dance is moot; reset with KFLT_VOID. Fn-data
         * bodies don't typically `return` either, so cleanup happens
         * once at fn-body fall-through via kfl_emit_stmt_drain_root
         * below. */
        kfl_emit_stmt_reset_scopes(tmp_arena, KFLT_VOID);

        for (const KflcNode *s = fn->children; s; s = s->next) {
            if (kfl_emit_stmt(out, s, &body_ctx, diag, 4)) {
                /* error already raised */
            }
        }
        /* fn-data has no return path; emit end-of-fn cleanup
         * unconditionally. The series_<kind> statement's static
         * caches hold *copies* of the vectors so this cleanup is
         * safe — see stmt.c emit. */
        kfl_emit_stmt_drain_root(out, 4);
        /* Hand the gathered series array to libk26plot. */
        fputs("    *_kfl_out_series = _kfl_all;\n", out);
        fputs("    *_kfl_n_out = (size_t)_kfl_count;\n", out);
        fputs("}\n\n", out);
    }

    /* Emit `fn world <name>` callbacks. Signature is
     *   extern "C" void <name>(<world_cxx> world, void *_kfl_user_data);
     * where <world_cxx> comes from the "world" opaque registration.
     * If "world" is not registered we emit `void *` and warn; the
     * generated code still compiles standalone (e.g. for the
     * round-trip test); a host linking it against an actual runtime
     * (libk26astro_rt) will register the proper opaque before
     * compilation. */
    for (const KflcNode *fn = form->children; fn; fn = fn->next) {
        if (fn->kind != KFLN_FN_WORLD) continue;
        const char *world_cxx = kflc_opaque_cxx("world");
        if (!world_cxx) {
            kflc_diag_warnf(diag, fn->line,
                "fn world %s: opaque type `world` is not registered; "
                "emitting `void *` as the handle type; a host runtime "
                "must register the proper opaque (kflc_opaque_register) "
                "before this code is linked against it",
                fn->name ? fn->name : "?");
            world_cxx = "void *";
        }
        fprintf(out,
            "extern \"C\" void %s(%s world, void *_kfl_user_data) {\n"
            "    (void)world;\n"
            "    (void)_kfl_user_data;\n",
            fn->name ? fn->name : "_anon",
            world_cxx);

        KflcExprBinding *live = NULL;
        int live_n = 0, live_cap = 0;
        collect_form_args(form, tmp_arena, &live, &live_n, &live_cap);

        /* Pre-pass: collect all `astro_body NAME` declarations in
         * this fn body, emit `int _kfl_body_<NAME>_idx` at fn
         * prologue so subsequent statements can deref via the cached
         * idx instead of runtime `k26astro_world_find_body`. */
        const char **known_body_names = NULL;
        int n_known_bodies = 0, cap_known = 0;
        for (const KflcNode *s = fn->children; s; s = s->next) {
            if (s->kind != KFLN_STMT_ASTRO_BODY || !s->name) continue;
            /* Dedup: same name twice in source → emit one idx slot. */
            int dup = 0;
            for (int i = 0; i < n_known_bodies; i++) {
                if (strcmp(known_body_names[i], s->name) == 0) { dup = 1; break; }
            }
            if (dup) continue;
            if (n_known_bodies >= cap_known) {
                cap_known = cap_known ? cap_known * 2 : 8;
                known_body_names = kflc_arena_alloc(tmp_arena,
                    (size_t)cap_known * sizeof(const char *));
            }
            known_body_names[n_known_bodies++] = s->name;
        }
        for (int i = 0; i < n_known_bodies; i++) {
            fprintf(out, "    int _kfl_body_%s_idx = -1;\n",
                    known_body_names[i]);
            fprintf(out, "    (void)_kfl_body_%s_idx;\n",
                    known_body_names[i]);
        }

        /* Pre-pass: collect all `let` / `const` bindings in this fn
         * world body so expressions can resolve forward references to
         * them. Mirrors the existing KFLN_FN body emit. Without this,
         * `let x: double = ...` followed by `let y = f(x)` reports
         * "expr: unknown identifier x" because the body_ctx binding
         * table doesn't include x yet when the y expression is
         * evaluated. */
        for (const KflcNode *s = fn->children; s; s = s->next) {
            collect_let_bindings(s, tmp_arena, &live, &live_n, &live_cap);
        }

        /* Register the implicit `world` parameter as an opaque
         * binding so KFL expressions can reference it (e.g.
         * `enable_outer_planets(world, eph)`). The C signature is
         * `void fn(K26AstroWorld *world, void *_kfl_user_data)`;
         * KFL sees `world` as type opaque<world>. */
        if (live_n >= live_cap) {
            live_cap = live_cap ? live_cap * 2 : 8;
            KflcExprBinding *grown = kflc_arena_alloc(tmp_arena,
                (size_t)live_cap * sizeof(KflcExprBinding));
            if (live_n > 0) {
                memcpy(grown, live, (size_t)live_n * sizeof(KflcExprBinding));
            }
            live = grown;
        }
        /* Zero-init then populate so lifetime_qualifier defaults to
         * KFL_LQ_NONE. The implicit world binding is unqualified;
         * fn-world body code that wants ownership semantics declares
         * a separate own/borrow binding via `let w: borrow world =
         * world` or similar. */
        memset(&live[live_n], 0, sizeof(KflcExprBinding));
        live[live_n].name             = "world";
        live[live_n].type             = KFLT_OPAQUE;
        live[live_n].type_subtype     = "world";
        live[live_n].is_form_arg      = 0;
        live[live_n].is_observed_cell = 0;
        live[live_n].lifetime_qualifier = KFL_LQ_NONE;
        /* Implicit `world` binding is not moved and not a borrow. */
        live[live_n].moved_from        = 0;
        live[live_n].borrow_source_idx = -1;
        live[live_n].scope_depth       = 0;
        live_n++;

        KflcExprCtx body_ctx;
        memset(&body_ctx, 0, sizeof(body_ctx));
        body_ctx.bindings           = live;
        body_ctx.n_bindings         = live_n;
        body_ctx.fns                = user_fn_arr;
        body_ctx.n_fns              = n_user_fns;
        body_ctx.known_body_names   = known_body_names;
        body_ctx.n_known_body_names = n_known_bodies;
        body_ctx.form               = form;   /* arena/borrow lookup */
        body_ctx.headless           = 1;  /* observe → stdout */

        kfl_emit_stmt_reset_scopes(tmp_arena, KFLT_VOID);

        for (const KflcNode *s = fn->children; s; s = s->next) {
            kfl_emit_stmt(out, s, &body_ctx, diag, 4);
        }
        kfl_emit_stmt_drain_root(out, 4);
        fputs("}\n\n", out);
    }


    /* main(). */
    fputs("int main(int argc, char **argv) {\n", out);
    if (nargs > 0) {
        fputs("    for (int i = 1; i < argc; i++) {\n", out);
        for (const KflcNode *c = form->children; c; c = c->next) {
            if (c->kind != KFLN_ARG) continue;
            /* C1: skip CLI override for readonly args. The global is
             * `const`-qualified and the assignment would fail to
             * compile; even if it didn't, exposing a `--<arg>` flag
             * for a readonly arg would mis-advertise the contract.
             * The arg keeps its declared default for the lifetime
             * of the process. */
            if (c->flags & KFL_NF_READONLY) continue;
            if (c->position.kind == KFLV_STR) {
                fprintf(out,
                    "        if (strcmp(argv[i], \"--%s\") == 0 && i+1 < argc)"
                    " { kfl_arg_%s = argv[++i]; continue; }\n",
                    c->name, c->name);
            } else {
                fprintf(out, "        if (strcmp(argv[i], \"--%s\") == 0)    kfl_arg_%s = true;\n",
                        c->name, c->name);
                fprintf(out, "        if (strcmp(argv[i], \"--no-%s\") == 0) kfl_arg_%s = false;\n",
                        c->name, c->name);
            }
        }
        fputs("    }\n", out);
    }
    /* Initialise any form-level arenas before running program code
     * that allocates from them. */
    if (form_has_arena(form)) {
        fputs("    _kfl_init_arenas();\n", out);
    }
    /* No window. Run each `fn world` against a freshly created
     * astronomical world (PORTABLE mode for bit-identical
     * reproducibility across CPUs), then tear it down. The world
     * body's `observe` statements print to stdout (see stmt.c).
     * Programs with no `fn world` fall through to a clean exit. */
    int any_world = 0;
    for (const KflcNode *fn = form->children; fn; fn = fn->next) {
        if (fn->kind == KFLN_FN_WORLD) { any_world = 1; break; }
    }
    if (any_world) {
        fputs("    K26AstroWorld *_kfl_world = k26astro_world_create("
              "K26ASTRO_MODE_PORTABLE, K26ASTRO_COORDS_SECTOR_GRID);\n",
              out);
        fputs("    if (!_kfl_world) { "
              "fprintf(stderr, \"kflc: world create failed\\n\"); "
              "return 1; }\n", out);
        for (const KflcNode *fn = form->children; fn; fn = fn->next) {
            if (fn->kind != KFLN_FN_WORLD || !fn->name) continue;
            fprintf(out, "    %s(_kfl_world, nullptr);\n", fn->name);
        }
        fputs("    k26astro_world_destroy(_kfl_world);\n", out);
    }
    /* Render every plot widget to <name>.png + <name>.svg. */
    if (form_has_plot(form)) {
        int plot_counter = 0;
        emit_headless_plots_(out, form, &plot_counter);
    }
    /* Call the conventional entry fn `run` (parameterless) if the
     * program defines one — the place for a pure-compute program to
     * do its work and print results. */
    for (const KflcNode *fn = form->children; fn; fn = fn->next) {
        if (fn->kind == KFLN_FN && fn->name &&
            strcmp(fn->name, "run") == 0) {
            fputs("    (void)run();\n", out);
            break;
        }
    }
    fputs("    return 0;\n", out);
    fputs("}\n", out);

    kflc_arena_release(tmp_arena);
    return diag->errors ? 1 : 0;
}
