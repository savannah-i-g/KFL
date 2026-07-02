#ifndef KFLC_H
#define KFLC_H

/* kflc — KFL compiler.
 *
 * Reads a .kfl source file, produces C++ source, and invokes the
 * system C++ toolchain to produce a native binary. See GRAMMAR.md
 * for the language reference.
 *
 * Public API: parse a file into an AST, emit C++ from the AST, drive
 * the system compiler to produce a binary.
 *
 * All allocations come from a per-compile arena. Free the whole arena
 * with kflc_arena_release() at end of compile; do not free individual
 * nodes / strings.
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Arena -------------------------------------------------------- */

typedef struct KflcArena KflcArena;

KflcArena  *kflc_arena_create(void);
void        kflc_arena_release(KflcArena *a);
void       *kflc_arena_alloc(KflcArena *a, size_t n);
char       *kflc_arena_strdup(KflcArena *a, const char *s);
char       *kflc_arena_strndup(KflcArena *a, const char *s, size_t n);

/* ---- Value -------------------------------------------------------- */

typedef enum {
    KFLV_NONE = 0,
    KFLV_INT,
    KFLV_FLOAT,
    KFLV_STR,        /* quoted string, escapes already decoded */
    KFLV_IDENT,      /* bare identifier (e.g. font name) */
    KFLV_COLOR,      /* 0x00RRGGBB */
    KFLV_SIZE,       /* WxH */
    KFLV_SHORTCUT,   /* opaque shortcut string; emitter converts to FL_* int */
    KFLV_ALIGN,      /* alignment bitmask (legacy attribute value) */
} KflcValueKind;

typedef struct {
    KflcValueKind kind;
    union {
        long      i;
        double    f;
        char     *s;        /* heap-allocated in arena (STR/IDENT/SHORTCUT) */
        uint32_t  rgb;
        struct { int w, h; } size;
        int       align;    /* see kflc_align_bits */
    } u;
} KflcValue;

/* ---- AST ---------------------------------------------------------- */

typedef enum {
    KFLN_FORM = 1,        /* top-level form */
    KFLN_WIDGET,          /* leaf widget (plot figure) */
    KFLN_ARG,             /* arg <name> [default <value>] */
    /* Inline functions + statement-level constructs. Statement nodes
     * live inside KFLN_FN bodies.
     */
    KFLN_FN,              /* fn <ret_type> <name>(<args>) ... end */
    KFLN_FN_ARG,          /* typed function parameter (in KFLN_FN.children) */
    KFLN_STMT_LET,        /* let <name>: <type> = <expr> */
    KFLN_STMT_CONST,      /* const <name>: <type> = <expr> */
    KFLN_STMT_ASSIGN,     /* <name> = <expr> */
    KFLN_STMT_RETURN,     /* return [<expr>] */
    KFLN_STMT_EXPR,       /* expression-as-statement */
    KFLN_STMT_IF,         /* if <expr> ... [else] end (fn-body) */
    KFLN_STMT_WHILE,      /* while <expr> ... end (fn-body) */
    /* print <arg>[, <arg>...] — write each argument to stdout followed
     * by a newline. Each arg is a string literal (printed verbatim) or a
     * numeric expression (printed as a double). Args are held as child
     * nodes: a string arg carries KFLV_STR in `position`, a numeric arg
     * carries its `expr`. */
    KFLN_STMT_PRINT,
    /* Self-contained .kfl. `fn data <name>` declares an inline plot
     * data provider with the libk26plot callback signature.
     * `xs[i] = expr` writes a vector element. `series_<kind> <label>
     * <xs> <ys>` packages two vectors into a K26PSeries and assigns it
     * to the plot's out-parameters. */
    KFLN_FN_DATA,         /* fn data <name> ... end */
    KFLN_STMT_INDEX_ASSIGN, /* <name>[<expr>] = <expr> (single-level vector) */
    KFLN_STMT_SERIES,     /* series_<kind> <label> <xs_ident> <ys_ident> */
    /* Generalised lvalue assignment. Carries the full LHS KflcExpr in
     * `expr` and the RHS in `expr2`. Used for any lvalue shape that
     * doesn't fit STMT_ASSIGN (bare ident) or STMT_INDEX_ASSIGN
     * (single-level vector index, kept distinct so the KFLC_VEC_W
     * bounds-check macro applies there). Today this means matrix
     * `m[row][col] = expr`; future struct-field writes (`p.x = expr`)
     * follow the same node shape. */
    KFLN_STMT_LVALUE_ASSIGN,
    /* Astro / simulation surface. `fn world <name>` declares a
     * runtime callback that an external runtime (libk26astro_rt or any
     * other host) calls to populate / advance an opaque "world"
     * handle. Codegen mirrors `fn scene`: extern "C" callback at file
     * scope, body is plain statements (let/const/assign/expr +
     * registered astro builtins). The world handle type is whatever
     * the "world" opaque registration maps to (typically
     * `K26AstroWorld *`); if "world" is not registered at emit time,
     * the emitter reports an error. */
    KFLN_FN_WORLD,
    /* Form-level `tick fn=<handler> interval_ms=<n>` declares a
     * timer callback. Emits a registration call to the existing
     * kflc_plot3d_set_tick_cb seam at compiler-emitted form-ctor
     * scope. `name` carries the handler identifier; the interval
     * lives in attrs as an int literal under "interval_ms". */
    KFLN_TICK,
    /* Form-level `frame <name> inertial` or `frame <name>
     * body_fixed body=<ident>` declares a frame newtype. The compiler
     * registers an opaque type whose KFL name is the frame
     * identifier; attempts to mix frames at call sites produce a
     * type mismatch. The `attrs` list carries kind=("inertial" |
     * "body_fixed") and optional `body=<ident>`. */
    KFLN_FRAME,
    /* Form-level `epoch <name> <ISO-8601> <scale>` binds an
     * ISO-8601 timestamp + time-scale ID (TAI/UTC/UT1/TT/TDB) to a
     * named constant. The parser converts to (days_since_J2000,
     * seconds_of_day, scale) at compile time so the runtime never
     * re-parses ISO-8601 strings. `name` carries the binding name;
     * attrs hold the parsed components. */
    KFLN_EPOCH_LITERAL,
    /* Astro statements live inside `fn world` bodies. Parsed by
     * stmt.c, serialised by serialize.c, emitted by stmt.c's
     * kfl_emit_stmt dispatch into calls against libk26astro_rt
     * builtins. Named-argument form (`name=expr`) is accepted in
     * statement-arg context only; form-attribute syntax (`tick fn=foo`)
     * is unaffected. */
    KFLN_STMT_ASTRO_BODY,  /* astro_body <name> { gm=..., pos=..., vel=..., ... } */
    KFLN_STMT_STEP,        /* step <dt_expr> */
    KFLN_STMT_PROPAGATE,   /* propagate <body_ident> for <dt_expr> */
    KFLN_STMT_FOR_EACH,    /* for_each <ident> in <world_ident> { <stmts> } */
    KFLN_STMT_OBSERVE,     /* observe <target_ident> from <observer_ident> mode=<kw> */
    /* KFL memory model nodes. Form-level arena declaration and
     * fn-prologue allocator binding. The arena is a bump-allocator
     * region whose lifetime extends across all fn invocations that
     * declare `allocator = <arena_name>`; capacity is the upper bound
     * on total live allocations within the region. Resets implicitly
     * on each fn entry when fn-scoped, or never if form-scoped
     * (controlled by reset_mode attr; default fn-scoped).
     *
     * KFLN_ARENA `name` carries the arena identifier; attrs hold
     * `capacity` (size_t, in bytes — KFLV_INT post-suffix-decoded),
     * optional `reset_mode` (KFLV_IDENT "fn" | "manual" | "form").
     *
     * KFLN_ALLOCATOR_BIND is parsed inside fn bodies as a prologue
     * attribute (must appear before any non-binding statement in the
     * fn body). `name` is the arena identifier being bound; the parser
     * links this to the form-level KFLN_ARENA so the emit knows the
     * arena's C symbol name. */
    KFLN_ARENA,
    KFLN_ALLOCATOR_BIND
} KflcNodeKind;

/* Type system. The base scalar kinds are joined by KFLT_VECTOR /
 * KFLT_MATRIX (K26CVector / K26CMatrix value structs with heap-owned
 * data buffers, freed at fn-scope exit + before each return) and by
 * KFLT_OPAQUE, a named handle type registered by external libraries
 * through `.kflbi` builtin manifests. Opaque values are pointer-sized;
 * the KFL surface treats them as nominal (no field access, only
 * function calls), and the emitter maps the KFL name (e.g. "world",
 * "body", "ephem") to a registered C++ name (e.g. "K26AstroWorld *",
 * "K26AstroBody *"). */
typedef enum {
    KFLT_VOID   = 0,
    KFLT_DOUBLE = 1,
    KFLT_INT    = 2,
    KFLT_BOOL   = 3,
    KFLT_STRING = 4,        /* `const char *` in emitted C */
    KFLT_VECTOR = 5,
    KFLT_MATRIX = 6,
    KFLT_OPAQUE = 7         /* named external handle */
} KflcType;

/* Lifetime qualifier orthogonal to KflcType. Set on KflcNode (for
 * let/const/fn-arg/fn-return slots) and on KflcExprBinding (for
 * in-scope variable lookups). NONE preserves value-typed scalars and
 * opaque handles emitted as `T *` with no ownership tracking.
 *
 *   OWN    , value, owned. On opaques, assignment is a pointer move
 *            with the RHS source nulled out at C-emit time; only one
 *            live binding to the underlying object at a time. On
 *            scalars, no effect (compile-time error at emit time).
 *   BORROW , value, borrowed. On opaques, emitted as `T *` but with
 *            a compile-time scope check: a borrow binding cannot
 *            escape the source binding's scope (return, store in
 *            arg, capture in arena-allocated struct). On scalars,
 *            no effect.
 *   PTR    , pointer constructor. Adds one level of C-emit
 *            indirection: `ptr f64` emits `double *`, `ptr body`
 *            emits `K26AstroBody *` (same as opaque-default for
 *            opaques, but explicitly raw, no ownership). Possibly
 *            NULL; receiver must null-check.
 *
 * NONE is the default; programs that do not annotate bindings parse
 * with NONE on every binding. */
typedef enum {
    KFL_LQ_NONE   = 0,
    KFL_LQ_OWN    = 1,
    KFL_LQ_BORROW = 2,
    KFL_LQ_PTR    = 3
} KflcLifetimeQualifier;

/* Recognise a lifetime-qualifier keyword. Returns the enum value or
 * -1 on no match. Called by parser.c at the type-prefix slot of
 * let/const/fn-arg/fn-return parses (NOT by kflc_type_from_str —
 * lifetime qualifiers are a separate parse axis). */
int kflc_lifetime_qualifier_from_str(const char *name);

/* Inverse — serializer-side. Returns "own"/"borrow"/"ptr" or NULL
 * for KFL_LQ_NONE (so serialize.c omits the prefix entirely when
 * none). */
const char *kflc_lifetime_qualifier_kfl_str(KflcLifetimeQualifier lq);

/* Opaque registry. Lookups are process-global; populate at compiler
 * init time (manifest loader scans K26_KFL_BUILTIN_PATH).
 * `kfl_name` is the surface identifier in .kfl source; `cxx_name` is
 * the C++ type emitted in generated code (typically a pointer to a
 * lib-provided struct, e.g. "K26AstroWorld *"). Both strings must
 * outlive the compiler invocation — pass arena-allocated or static
 * literals.
 *
 * Returns 0 on success, nonzero on collision or NULL inputs. */
int         kflc_opaque_register(const char *kfl_name, const char *cxx_name);
/* Extended registration with an optional copy-function name.
 * `copy_fn_cxx` is the C symbol of a deep-copy helper with signature
 *   `<cxx_name> (<copy_fn_cxx>)(<cxx_name>)`
 * (one argument, source; returns a freshly-allocated duplicate that
 * the receiver owns). Used by emit when the KFL surface contains
 * `<binding>.copy()`. NULL `copy_fn_cxx` means the opaque has no
 * copy support — `.copy()` calls on such bindings are compile-time
 * errors. Otherwise same as kflc_opaque_register. */
int         kflc_opaque_register_with_copy(const char *kfl_name,
                                           const char *cxx_name,
                                           const char *copy_fn_cxx);
/* Look up an already-registered opaque by KFL name. Returns the
 * cxx_name string passed to kflc_opaque_register, or NULL if the
 * name isn't registered. */
const char *kflc_opaque_cxx(const char *kfl_name);
/* Returns the copy-function C symbol registered for `kfl_name`, or
 * NULL if the opaque has no copy support (or isn't registered).
 * Consulted by emit when lowering `<binding>.copy()`. */
const char *kflc_opaque_copy_fn(const char *kfl_name);
/* Clear the registry. Call between independent compile invocations
 * if reusing the same process; called automatically at the start of
 * kflc_register_builtins_reload(). */
void        kflc_opaque_clear(void);

/* Dynamic builtin registry. Companion to kflc_opaque_register: a
 * library publishes its function bindings here so the expression
 * emitter can resolve KFL-level call names to their C++ counterparts.
 * Strings must outlive the compiler invocation (static literals or
 * arena-allocated). Arity is checked at emit time. Returns 0 on
 * success, 1 on NULL inputs, 2 on collision, 3 on table full. */
int         kflc_register_builtin(const char *kfl_name,
                                  const char *cxx_name,
                                  int         arity);
/* Empty the dynamic registry. Static built-ins (libm, libk26compute)
 * are unaffected — they live in a separate const table. */
void        kflc_clear_builtins(void);

/* Astro bootstrap. Registers the `world` / `starfield` opaque types
 * and all libk26astro_rt + libk26astro_render builtins (mirrors the
 * .kflbi manifests as an in-process fallback for embedders that do
 * not ship the manifest path). Idempotent; safe to call from any
 * embedding tool. */
void        kflc_register_astro_builtins(void);

/* Builtin manifest loader. Scans K26_KFL_BUILTIN_PATH (default
 * /usr/share/kflc/builtins/) for *.kflbi files (sorted) and registers
 * their opaque types + builtins via the C API. Silent no-op if the
 * directory is missing; the hardcoded astro registrar remains the
 * fallback for host / no-install builds. Collisions are reported on
 * stderr and resolved as "first registration wins". */
void        kflc_load_manifests(void);

/* Map a (KflcType, opaque-subtype-name) pair to its C++ surface
 * type. `opaque_name` may be NULL when `t` is not KFLT_OPAQUE;
 * required (non-NULL) when `t == KFLT_OPAQUE`. NULL return means
 * unknown type or unregistered opaque name. */
const char *kflc_type_cxx(KflcType t, const char *opaque_name);
/* Parse a type identifier into a KflcType. Returns the KflcType
 * (non-negative) on a recognised scalar/vector/matrix name, or
 * KFLT_OPAQUE when the name matches a registered opaque type — in
 * which case `*out_opaque_name` is set to the registry's stored
 * `kfl_name` pointer. Returns -1 on no match. `out_opaque_name` may
 * be NULL, in which case opaque types are not recognised. */
int         kflc_type_from_str(const char *name,
                               const char **out_opaque_name);
/* Inverse of kflc_type_from_str — map a (KflcType, opaque-subtype)
 * back to its KFL surface name. For KFLT_OPAQUE returns
 * `opaque_name`. Used by the serializer so round-tripped .kfl source
 * uses KFL names, not the C++ surface that kflc_type_cxx returns.
 * NULL on unknown. */
const char *kflc_type_kfl_str(KflcType t, const char *opaque_name);

/* Legacy container discriminant. The layout grammar that produced
 * container nodes has been removed; this enum and the KflcNode
 * `container` field it types are retained only so structural AST
 * comparison keeps a stable field to diff. No parse produces these. */
typedef enum {
    KFL_CONT_VBOX = 1,
    KFL_CONT_HBOX,
    KFL_CONT_GRID,
    KFL_CONT_GROUP,
    KFL_CONT_PANEL,
    KFL_CONT_TILE_H,
    KFL_CONT_TABS,
    KFL_CONT_TAB,
    KFL_CONT_SCROLL
} KflcContainerKind;

typedef enum {
    /* The only leaf widget: a plot figure drawn via libk26plot. Its
     * `data <fn>` provider populates one or more series; artifact
     * emission renders each plot to <name>.png / <name>.svg. */
    KFL_W_PLOT = 1
} KflcWidgetKind;

/* Forward decl; full definition appears later. KflcNode holds an
 * optional KflcExpr * for statement values in fn bodies. */
typedef struct KflcExpr KflcExpr;

typedef struct KflcAttr KflcAttr;
struct KflcAttr {
    char       *name;     /* arena-owned */
    KflcValue   value;
    int         line;
    /* Optional pre-parsed expression carried alongside `value`.
     * Used by attributes whose right-hand side is a paren-balanced
     * expression (e.g. mesh `label <call(...)>`). NULL for simple
     * literal / identifier attributes. */
    KflcExpr   *expr;
    KflcAttr   *next;
};

typedef struct KflcNode KflcNode;
struct KflcNode {
    KflcNodeKind        kind;
    int                 line;

    /* Form: form name (e.g. "ORBITS").
     * Widget: optional plot name (ident).
     * Fn: function name. Fn_arg / let / const / assign: declared name.
     */
    char               *name;

    /* Positional first value:
     *   plot p_inner ...  -> position is IDENT "p_inner" (the name).
     */
    KflcValue           position;

    KflcContainerKind   container;   /* legacy; unused (see KflcContainerKind) */
    KflcWidgetKind      widget;      /* for KFLN_WIDGET */

    /* Scalar type of a fn return / fn arg / let / const binding.
     * Statement nodes set this; non-statement nodes leave it KFLT_VOID. */
    KflcType            type;
    /* Opaque-subtype name for KFLT_OPAQUE nodes; points into the
     * opaque registry (kflc_opaque_register) and is the KFL surface
     * name (e.g. "world", "body", "ephem"). NULL when `type` isn't
     * KFLT_OPAQUE. */
    const char         *type_subtype;
    /* Lifetime qualifier on the type slot. KFL_LQ_NONE for unannotated
     * bindings (default). KFL_LQ_OWN / BORROW for lifetime-tracked
     * opaques; KFL_LQ_PTR for raw-pointer indirection. See
     * KflcLifetimeQualifier docs above for the full semantic axes. */
    KflcLifetimeQualifier lifetime_qualifier;

    KflcAttr           *attrs;       /* linked list */
    KflcNode           *children;    /* first child (then-branch for
                                       *  KFLN_STMT_IF / KFLN_STMT_WHILE; body
                                       *  for KFLN_FN) */
    KflcNode           *else_children; /* else-branch (KFLN_STMT_IF only) */
    KflcNode           *next;        /* next sibling */

    /* Expression carried by statement nodes (LET / CONST / ASSIGN /
     * RETURN / EXPR / IF condition / WHILE condition). NULL on
     * non-statement nodes and on bare `return`. */
    KflcExpr           *expr;

    /* Secondary expression slot used by STMT_INDEX_ASSIGN to carry
     * the index. KflcNode.name = base vector ident; expr = RHS
     * (value); expr2 = index. NULL on other nodes. */
    KflcExpr           *expr2;

    /* Per-node bit flags. Currently used by KFLN_ARG to mark `readonly`
     * form-args (emitted as `static const` and rejected as
     * `target`/`cmd toggle` write destinations). Future uses
     * anticipated: deprecation markers, `persistent` arg flag. Cleared
     * to 0 by new_node so unused-by-this-kind nodes are unaffected. */
    uint16_t            flags;
};

/* KflcNode::flags bit definitions. Only meaningful for the node kind
 * specified in the comment; other kinds ignore the bits. */
#define KFL_NF_READONLY        (1u << 0)   /* KFLN_ARG: arg is read-only */
/* Set on KFLN_ARG nodes that are referenced by any of the three
 * reactive-binding mechanisms: `target <arg>`, `bind_visible <arg>` /
 * `bind_hidden <arg>`, or `cmd toggle target <arg>`. Marks the arg as
 * needing a `KflCell_<name>` substrate at form scope so writes can
 * fan out to subscribers (see 12-binding-model.md §3.4). Set post-parse
 * by a single AST walk in the emitter before form-arg emission. The
 * walker (kflc_emit_lvalue) reads this to decide whether to append an
 * `_kfl_cell_<name>.epoch++` after each write to the arg's storage
 * (12-binding-model.md §4.3). */
#define KFL_NF_OBSERVED_CELL   (1u << 1)   /* KFLN_ARG: cell substrate active */

/* ---- Expressions ------------------------------------------------- */

/* Expression sub-language. Parsed from a string source (the contents
 * of an `expression "..."` attribute, or the body of an inline `fn`
 * statement). The expression is itself an AST that the emitter
 * translates to a C++ expression. The whole language is statically
 * typed and pure-functional; identifier resolution happens at codegen
 * time, against form args + fn parameters + the libm / libk26compute
 * whitelist.
 *
 * Supported grammar includes literals, unary +/-/!, the arithmetic /
 * comparison / logical operators, parenthesised sub-expressions,
 * function calls, vector literals, and indexed reads (scalar from a
 * vector binding). See expr.c for the precedence climbing detail.
 */

typedef enum {
    KFLE_INT_LIT = 1,
    KFLE_FLOAT_LIT,
    KFLE_IDENT,           /* bare identifier; resolved at emit time */
    KFLE_UNARY,
    KFLE_BINARY,
    KFLE_CALL,
    /* Vector literal `[1.0, 2.0, 3.0]`. Holds a list of
     * sub-expressions, each evaluated to a double at construction
     * time. Only valid as a let/const initialiser (must produce a
     * vector lvalue, not embed in a scalar context). */
    KFLE_VEC_LIT,
    /* Indexed read into a vector binding. `xs[i]` evaluates to a
     * `double`. The base sub-expression MUST resolve to a vector
     * binding by name (no chained / anonymous bases yet). */
    KFLE_INDEX
} KflcExprKind;

typedef enum {
    KFLOP_NEG = 1,    /* unary - */
    KFLOP_POS,        /* unary + (no-op, kept for AST symmetry) */
    KFLOP_NOT,        /* unary ! */
    KFLOP_ADD,
    KFLOP_SUB,
    KFLOP_MUL,
    KFLOP_DIV,
    KFLOP_MOD,
    KFLOP_LT,
    KFLOP_LE,
    KFLOP_GT,
    KFLOP_GE,
    KFLOP_EQ,
    KFLOP_NE,
    KFLOP_AND,
    KFLOP_OR
} KflcOp;

struct KflcExpr {
    KflcExprKind kind;
    int          line;
    union {
        long       i;
        double     f;
        char      *ident;
        struct { KflcOp op; KflcExpr *operand; } un;
        struct { KflcOp op; KflcExpr *lhs, *rhs; } bin;
        struct { char *name; KflcExpr **args; int n_args; } call;
        /* KFLE_VEC_LIT: vector literal `[<expr>, <expr>, ...]`. */
        struct { KflcExpr **elems; int n_elems; } vec;
        /* KFLE_INDEX: <base>[<index>]. base is typically a KFLE_IDENT
         * referring to an in-scope vector binding. */
        struct { KflcExpr *base; KflcExpr *idx; } index;
    } u;
};

/* ---- Diagnostics -------------------------------------------------- */

typedef struct {
    const char *path;     /* source path for error prefix */
    int         errors;   /* count of errors raised on this context */
    int         warnings; /* count of warnings raised on this context */
    FILE       *errstream;
} KflcDiag;

void kflc_diag_init(KflcDiag *d, const char *path, FILE *errstream);

/* Plain line-only error. Kept for back-compat with existing call-sites
 * that haven't been updated to pass column information. Equivalent to
 * `kflc_diag_errorf_col(d, line, 0, fmt, ...)`. */
void kflc_diag_errorf(KflcDiag *d, int line, const char *fmt, ...)
    __attribute__((format(printf, 3, 4)));

/* Error with column. Prints `<path>:<line>:<col>: error: ...` when
 * `col > 0`, falling back to `<path>:<line>: error: ...` when `col`
 * is zero or negative. Lexer + parser sites that track column should
 * prefer this; emit-time diagnostics that only know line stay on the
 * line-only variant above. */
void kflc_diag_errorf_col(KflcDiag *d, int line, int col,
                          const char *fmt, ...)
    __attribute__((format(printf, 4, 5)));

/* Non-error diagnostics. Increment `d->warnings`, not `d->errors`,
 * so a clean parse can still complete with warnings present. Used by
 * deprecation notices (e.g. legacy `box`/`label` attribute spelling),
 * reserved-word shadowing, and other future-fragile constructs. */
void kflc_diag_warnf(KflcDiag *d, int line, const char *fmt, ...)
    __attribute__((format(printf, 3, 4)));
void kflc_diag_warnf_col(KflcDiag *d, int line, int col,
                         const char *fmt, ...)
    __attribute__((format(printf, 4, 5)));

/* ---- Parse -------------------------------------------------------- */

/* Parse a .kfl source string. Returns the root KFLN_FORM node, or NULL
 * on error (errors reported via KflcDiag). The arena owns all returned
 * memory. */
KflcNode *kflc_parse(const char *src, size_t len,
                     KflcArena *arena, KflcDiag *diag);

/* Convenience: read a file into memory and parse. */
KflcNode *kflc_parse_file(const char *path,
                          KflcArena *arena, KflcDiag *diag);

/* ---- Dump (debug) ------------------------------------------------- */

void kflc_dump_node(FILE *out, const KflcNode *n, int indent);

/* ---- Emit (codegen) ----------------------------------------------- */

/* Emit C++ source for the parsed form into `out`. Returns 0 on success
 * or nonzero on emit error (errors reported via KflcDiag).
 *
 * Emits a standalone console/batch program: no widget or window
 * scaffolding, and a plain `int main()` that runs the program's compute
 * and `fn world` logic, routing results to stdout / files. The emitted
 * source compiles against KFL_Stack alone. */
int kflc_emit_cxx(FILE *out, const KflcNode *form, KflcDiag *diag);

/* ---- Expression parse + emit ------------------------------------ */

/* Parse a self-contained expression from a NUL-terminated string.
 * Used by the compute widget's `expression` attribute and by inline
 * `fn` body parsing. Returns the root KflcExpr or NULL on parse
 * error; diagnostics go through KflcDiag.
 *
 * `line_base` is the source line of the surrounding statement; error
 * messages report this line. */
KflcExpr *kflc_parse_expr(const char *src, KflcArena *arena,
                          KflcDiag *diag, int line_base);

/* Per-emit context. `bindings` lists in-scope variables with their
 * KFL types — scalars emit with a `(double)` cast; vector/matrix
 * bindings emit as a bare struct lvalue so the call emitter can
 * splat the `.data` / `.n` fields into a libk26compute signature.
 * `fns` lists user-declared function names + arities; calls to these
 * are emitted as direct C++ calls (the user fn is emitted in the
 * same translation unit). Either array may be NULL with the
 * corresponding count set to 0. */
typedef struct {
    const char *name;
    int         arity;
} KflcExprFn;

typedef struct {
    const char *name;
    KflcType    type;
    /* Opaque-subtype name for KFLT_OPAQUE bindings. NULL when
     * `type` isn't KFLT_OPAQUE. */
    const char *type_subtype;
    /* Form-arg bindings (synthesised by collect_form_args) carry
     * is_form_arg = 1 so write-back assignments emit
     * `kfl_arg_<name>` rather than the bare identifier. Local let-
     * bindings stay at 0. */
    int         is_form_arg;
    /* Set on form-arg bindings whose KFLN_ARG node has
     * KFL_NF_OBSERVED_CELL. The lvalue walker (and the statement-case
     * assignment paths in stmt.c) emit a
     * `kfl_cell_notify(&_kfl_cell_<name>)` after writes to observed
     * args so subscribers (target subscribers, bind_visible
     * subscribers, show_hidden cascade, future reactive consumers)
     * fan out coherently. Non-form-arg bindings stay at 0. */
    int         is_observed_cell;
    /* Lifetime qualifier propagated from the binding's source
     * declaration. The emit-side expression walker consults this when
     * generating reads/writes so own/borrow/ptr semantics flow through
     * expression contexts. KFL_LQ_NONE on unannotated bindings. */
    KflcLifetimeQualifier lifetime_qualifier;
    /* Compile-time ownership tracking.
     *
     * `moved_from`: set when this `own` binding has been moved out of
     * via `move(x)` or own-to-own assignment. Subsequent reads in the
     * same scope are compile errors (intra-block enforcement only;
     * cross-block flow analysis is a future extension).
     *
     * `borrow_source_idx`: for `borrow` bindings, the index into the
     * same bindings[] array of the source binding. -1 for non-borrow
     * or when the source isn't in the local binding table (form-arg
     * borrows etc.).
     *
     * `scope_depth`: the BlockScope depth at which the binding was
     * declared. Compared against the source binding's depth on every
     * borrow read; error if the source has popped out of scope. The
     * no-cross-fn-borrow-escape rule additionally forbids a borrow
     * from outliving its source's enclosing fn. */
    int         moved_from;
    int         borrow_source_idx;
    int         scope_depth;
} KflcExprBinding;

typedef struct {
    const KflcExprBinding *bindings;
    int                    n_bindings;
    const KflcExprFn      *fns;
    int                    n_fns;
    /* Parse-time-known `astro_body NAME` declarations collected by
     * the fn-world prologue. Populated only inside `fn world` bodies;
     * NULL elsewhere. Each entry is a body name as declared in the KFL
     * source. `astro_body` emit assigns the corresponding `int
     * _kfl_body_<NAME>_idx` (declared at fn prologue);
     * `propagate`/`observe` emit looks up via idx instead of runtime
     * `k26astro_world_find_body`. Unknown names fall through to the
     * runtime path (handles `for_each` iters and runtime input
     * strings). */
    const char *const     *known_body_names;
    int                    n_known_body_names;
    /* Form pointer for arena reset_mode lookup at KFLN_ALLOCATOR_BIND
     * emit. Walk form->children to find the matching KFLN_ARENA decl
     * and read its reset_mode attr. NULL when ctx is built outside
     * any form (top-level expressions). */
    const KflcNode        *form;
    /* Borrow scope tracking: current block depth pushed by stmt.c's
     * BlockScope machinery. Updated as the emitter enters and leaves
     * block scopes (if/while/for bodies). The borrow resolver in
     * expr.c compares this against each binding's scope_depth to
     * decide if the binding is still in scope.
     *
     * Mutable in practice — emit code casts away const when updating
     * (the const-on-pointer signals "don't realloc bindings array"
     * which still holds; the scope int is metadata about traversal
     * state).
     *
     * Defaults to 0 (root fn scope). */
    int                    current_scope_depth;
    /* When non-zero, statement emitters route a result to stdout / a
     * file (e.g. `observe` prints its range + direction). Set to 1 by
     * kflc_emit_cxx for every emitted program. Defaults 0. */
    int                    headless;
} KflcExprCtx;

/* Emit a KflcExpr as a C++ expression to `out`. Returns 0 on
 * success, nonzero on emit error. `ctx` may be NULL for an empty
 * binding/fn-table context (e.g. compute widget expressions at
 * form scope outside any fn). */
int kflc_emit_expr(FILE *out, const KflcExpr *e,
                   const KflcExprCtx *ctx, KflcDiag *diag);

/* Returns non-zero if `name` was registered as a parse-time-known
 * body in the current fn-world body, via the `known_body_names` list
 * of `KflcExprCtx`. Returns 0 for unknown names (so callers fall
 * through to the runtime `find_body` path for `for_each` iters and
 * runtime input strings). NULL ctx or NULL name returns 0. */
int kflc_body_idx_known(const KflcExprCtx *ctx, const char *name);

/* kflc_emit_lvalue mode. The unified lvalue walker uses this to pick
 * between bounds-checked READ (KFLC_VEC_R / KFLC_MAT_R) and WRITE
 * (KFLC_VEC_W / KFLC_MAT_W) macro wrappers, and to know whether to
 * append an `_kfl_cell_<name>.epoch++` after observed-cell writes
 * (gated on KFL_NF_OBSERVED_CELL). Callers in KFLN_STMT_*_ASSIGN paths
 * pass KFLC_LV_WRITE; READ consumers (e.g. expression rvalues that
 * route through the unified walker rather than the older KFLE_INDEX
 * path) pass KFLC_LV_READ. */
typedef enum {
    KFLC_LV_READ  = 0,
    KFLC_LV_WRITE = 1
} KflcLvMode;

/* Emit a KflcExpr as a C++ lvalue (assignment target) or rvalue (when
 * mode is KFLC_LV_READ). Distinct from kflc_emit_expr because index
 * chains emit the `.data[...]` deref instead of the scalar-context
 * rvalue, and scalar idents skip the (double) cast that kflc_emit_expr
 * applies for bool/int rvalues. Supported shapes:
 *   - bare IDENT       → `kfl_arg_<name>` for form-args, bare for locals
 *   - vector[idx]      → `KFLC_VEC_W(<v>, idx)` / `KFLC_VEC_R(<v>, idx)`
 *   - matrix[row][col] → `KFLC_MAT_W(<m>, row, col)` / `KFLC_MAT_R(<m>, row, col)`
 * Anything else is a parse-time-or-emit-time error.
 *
 * `ctx` is required (binding lookup tells us whether the base is
 * a vector or matrix and whether an ident is a form-arg). */
int kflc_emit_lvalue(FILE *out, const KflcExpr *e,
                     const KflcExprCtx *ctx, KflcDiag *diag,
                     KflcLvMode mode);

/* Write a KflcExpr back out as KFL source text (the inverse of
 * kflc_parse_expr — used by serialize.c to round-trip fn bodies).
 * The output is human-readable and re-parsable by kflc_parse_expr. */
void kflc_expr_to_text(FILE *out, const KflcExpr *e);

/* ---- Serialize ---------------------------------------------------- */

/* Write a parsed form's AST back out as textual .kfl source. The
 * inverse of kflc_parse. Used by embedding tools (k26form, k26canvas)
 * that load a form, edit the AST in place, and save back to disk.
 *
 * Returns 0 on success, nonzero on I/O failure or an invalid form
 * root (form->kind must be KFLN_FORM). */
int kflc_serialize(FILE *out, const KflcNode *form);

/* Structural AST equality. Walks two trees comparing node kind, name,
 * position, container/widget/type discriminants, attrs (list-equal),
 * children + else_children + expr/expr2 subtrees, and the `flags`
 * bitfield. Source-line numbers are NOT compared (round-trip can
 * renumber). Float values use bitwise equality (==); if a round-trip
 * test surfaces a divergence here, the fix belongs in the serializer's
 * precision (use %.17g for floats), not in this helper. Returns 1 if
 * equal, 0 if not. NULLs compare equal to NULLs only.
 *
 * Used by the round-trip test harness to verify the parse / serialize
 * / parse semantic-equivalence promise. Reusable from k26form for
 * "is the in-memory form different from disk" checks, and from any
 * future LSP for change detection. */
int kflc_ast_equal(const KflcNode *a, const KflcNode *b);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* KFLC_H */
