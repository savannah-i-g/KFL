/* opaque.c: opaque-type registry.
 *
 * Holds a small process-global table mapping KFL surface names
 * ("world", "body", "ephem", ...) to their emitted C++ type names
 * ("K26AstroWorld *", "K26AstroBody *", ...). Populated by the
 * builtin-manifest loader (`builtin_loader.c`) when scanning
 * K26_KFL_BUILTIN_PATH at compiler startup; can also be populated
 * directly via kflc_opaque_register from in-process callers (tests,
 * embedding tools that bypass the manifest loader).
 *
 * Lookups are linear; the table is expected to stay small
 * (~20-50 entries across all loaded manifests). If it grows past
 * that, replace with a hash and revisit. */
#include "kflc.h"

#include <stdlib.h>
#include <string.h>

#define KFLC_OPAQUE_MAX 256

typedef struct {
    const char *kfl_name;
    const char *cxx_name;
    /* Optional deep-copy helper C symbol. NULL when the opaque has
     * no `.copy()` support; emit reports a compile-time error if
     * KFL source calls `.copy()` on such a binding. */
    const char *copy_fn_cxx;
} OpaqueEntry;

static OpaqueEntry g_table[KFLC_OPAQUE_MAX];
static int         g_count = 0;

int kflc_opaque_register(const char *kfl_name, const char *cxx_name)
{
    return kflc_opaque_register_with_copy(kfl_name, cxx_name, NULL);
}

int kflc_opaque_register_with_copy(const char *kfl_name,
                                   const char *cxx_name,
                                   const char *copy_fn_cxx)
{
    if (!kfl_name || !cxx_name) return 1;
    /* Reject collisions. The manifest loader is responsible for
     * coalescing duplicate names from sibling .kflbi files; an
     * unexpected duplicate here is a registration bug worth surfacing. */
    for (int i = 0; i < g_count; i++) {
        if (strcmp(g_table[i].kfl_name, kfl_name) == 0) return 2;
    }
    if (g_count >= KFLC_OPAQUE_MAX) return 3;
    g_table[g_count].kfl_name     = kfl_name;
    g_table[g_count].cxx_name     = cxx_name;
    g_table[g_count].copy_fn_cxx  = copy_fn_cxx;  /* may be NULL */
    g_count++;
    return 0;
}

const char *kflc_opaque_cxx(const char *kfl_name)
{
    if (!kfl_name) return NULL;
    for (int i = 0; i < g_count; i++) {
        if (strcmp(g_table[i].kfl_name, kfl_name) == 0)
            return g_table[i].cxx_name;
    }
    return NULL;
}

const char *kflc_opaque_copy_fn(const char *kfl_name)
{
    if (!kfl_name) return NULL;
    for (int i = 0; i < g_count; i++) {
        if (strcmp(g_table[i].kfl_name, kfl_name) == 0)
            return g_table[i].copy_fn_cxx;  /* may be NULL */
    }
    return NULL;
}

void kflc_opaque_clear(void)
{
    g_count = 0;
}

/* ------------------------------------------------------------------ */
/* Lifetime-qualifier surface.                                        */
/* Lives here rather than stmt.c because lifetime qualifiers compose  */
/* with the opaque registry's type semantics; both concerns share the */
/* "type-surface" naming layer that this file already owns.           */
/* ------------------------------------------------------------------ */

int kflc_lifetime_qualifier_from_str(const char *name)
{
    if (!name) return -1;
    if (strcmp(name, "own")    == 0) return KFL_LQ_OWN;
    if (strcmp(name, "borrow") == 0) return KFL_LQ_BORROW;
    if (strcmp(name, "ptr")    == 0) return KFL_LQ_PTR;
    return -1;
}

const char *kflc_lifetime_qualifier_kfl_str(KflcLifetimeQualifier lq)
{
    switch (lq) {
    case KFL_LQ_OWN:    return "own";
    case KFL_LQ_BORROW: return "borrow";
    case KFL_LQ_PTR:    return "ptr";
    case KFL_LQ_NONE:   /* fallthrough; serializer omits prefix */
    default:            return NULL;
    }
}
