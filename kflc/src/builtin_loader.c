/* builtin_loader.c: .kflbi manifest loader.
 *
 * At kflc startup, scan K26_KFL_BUILTIN_PATH (default
 * /usr/share/kflc/builtins/) for *.kflbi manifest files and register
 * their opaque types + builtins via the kflc C API. The loader is
 * additive: it runs alongside the hardcoded `builtin_register_astro`
 * registrar (which stays in place as a fallback for embedders that
 * don't ship the manifest path).
 *
 * Manifest format (v0.1, schema 1):
 *   schema  <int>                               -- required first line
 *   library <name>                              -- informational
 *   version <semver>                            -- informational
 *   opaque  <kfl_name> <cxx_type...>            -- opaque registration
 *   builtin <kfl_name> <cxx_name> <arity>       -- builtin registration
 *   # comments + blank lines ignored
 *
 * Error handling:
 *   - Each malformed line emits a stderr warning ("file:line:" prefix);
 *     loader continues to the next line / next manifest.
 *   - Schema != 1 → manifest rejected with a clear stderr message.
 *   - Collisions (same kfl_name registered twice) emit a stderr note;
 *     first registration wins (kflc_opaque/builtin_register semantics).
 *
 * String lifetime: kflc_opaque_register and kflc_register_builtin
 * require the string arguments to outlive the compiler invocation.
 * We strdup() everything once and never free — the loader's
 * registrations are effectively process-lifetime.
 *
 * Determinism: directory scan results are sorted alphabetically so
 * registration order is stable across runs / file systems. */
#include "kflc.h"

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define KFLC_BUILTIN_PATH_ENV     "K26_KFL_BUILTIN_PATH"
#define KFLC_BUILTIN_PATH_DEFAULT "/usr/share/kflc/builtins"
#define KFLC_KFLBI_SUFFIX         ".kflbi"
#define KFLC_MANIFEST_SCHEMA_MAX  1

static char *xstrdup_(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *p = (char *)malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n + 1);
    return p;
}

/* Trim leading whitespace; return pointer into the same buffer. */
static char *ltrim_(char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

/* Trim trailing whitespace in place. */
static void rtrim_(char *s)
{
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) s[--n] = '\0';
}

/* Consume the next whitespace-separated token from *cursor, writing a
 * NUL after it. Returns the token start or NULL if the cursor is at
 * end-of-line. */
static char *next_tok_(char **cursor)
{
    char *p = *cursor;
    while (*p && isspace((unsigned char)*p)) p++;
    if (!*p) { *cursor = p; return NULL; }
    char *tok = p;
    while (*p && !isspace((unsigned char)*p)) p++;
    if (*p) { *p = '\0'; p++; }
    *cursor = p;
    return tok;
}

/* Rest-of-line after stripping leading whitespace. Returns NULL if the
 * line has no remaining content. The returned pointer aliases into
 * the same line buffer; rtrim trailing whitespace before strdup. */
static char *rest_of_line_(char **cursor)
{
    char *p = *cursor;
    while (*p && isspace((unsigned char)*p)) p++;
    if (!*p) return NULL;
    return p;
}

static int parse_manifest_(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "kflc: warning: cannot open manifest %s: %s\n",
                path, strerror(errno));
        return -1;
    }

    int schema_seen = 0;
    int errors      = 0;
    char buf[1024];
    int line_no = 0;

    while (fgets(buf, sizeof buf, f)) {
        line_no++;
        /* Strip comment + trailing whitespace. */
        char *hash = strchr(buf, '#');
        if (hash) *hash = '\0';
        rtrim_(buf);
        char *line = ltrim_(buf);
        if (!*line) continue;

        char *cursor = line;
        char *kind   = next_tok_(&cursor);
        if (!kind) continue;

        if (strcmp(kind, "schema") == 0) {
            char *tok = next_tok_(&cursor);
            if (!tok) {
                fprintf(stderr, "%s:%d: schema needs a version\n",
                        path, line_no); errors++; continue;
            }
            int v = atoi(tok);
            if (v <= 0 || v > KFLC_MANIFEST_SCHEMA_MAX) {
                fprintf(stderr,
                        "%s:%d: unsupported schema %d (this kflc accepts up to %d)\n",
                        path, line_no, v, KFLC_MANIFEST_SCHEMA_MAX);
                fclose(f);
                return -1;  /* whole file rejected */
            }
            schema_seen = 1;
        } else if (strcmp(kind, "library") == 0
                || strcmp(kind, "version") == 0) {
            /* Informational; consume + ignore. */
        } else if (strcmp(kind, "opaque") == 0) {
            if (!schema_seen) {
                fprintf(stderr, "%s:%d: opaque before schema\n",
                        path, line_no); errors++; continue;
            }
            char *kfl_name = next_tok_(&cursor);
            char *cxx_rest = rest_of_line_(&cursor);
            if (!kfl_name || !cxx_rest) {
                fprintf(stderr, "%s:%d: opaque needs kfl_name + cxx_type\n",
                        path, line_no); errors++; continue;
            }
            rtrim_(cxx_rest);
            char *kfl_dup = xstrdup_(kfl_name);
            char *cxx_dup = xstrdup_(cxx_rest);
            if (!kfl_dup || !cxx_dup) { errors++; continue; }
            int rc = kflc_opaque_register(kfl_dup, cxx_dup);
            if (rc != 0 && rc != 2 /* collision */) {
                fprintf(stderr, "%s:%d: opaque register %s -> %s rc=%d\n",
                        path, line_no, kfl_dup, cxx_dup, rc);
                errors++;
            }
        } else if (strcmp(kind, "builtin") == 0) {
            if (!schema_seen) {
                fprintf(stderr, "%s:%d: builtin before schema\n",
                        path, line_no); errors++; continue;
            }
            char *kfl_name = next_tok_(&cursor);
            char *cxx_name = next_tok_(&cursor);
            char *arity_s  = next_tok_(&cursor);
            if (!kfl_name || !cxx_name || !arity_s) {
                fprintf(stderr,
                        "%s:%d: builtin needs kfl_name + cxx_name + arity\n",
                        path, line_no); errors++; continue;
            }
            int arity = atoi(arity_s);
            char *kfl_dup = xstrdup_(kfl_name);
            char *cxx_dup = xstrdup_(cxx_name);
            if (!kfl_dup || !cxx_dup) { errors++; continue; }
            int rc = kflc_register_builtin(kfl_dup, cxx_dup, arity);
            if (rc != 0 && rc != 2 /* collision */) {
                fprintf(stderr, "%s:%d: builtin register %s -> %s/%d rc=%d\n",
                        path, line_no, kfl_dup, cxx_dup, arity, rc);
                errors++;
            }
        } else {
            fprintf(stderr, "%s:%d: unknown directive '%s'\n",
                    path, line_no, kind);
            errors++;
        }
    }
    fclose(f);
    return errors;
}

/* qsort comparator for char *strings. */
static int strptr_cmp_(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

void kflc_load_manifests(void)
{
    const char *path = getenv(KFLC_BUILTIN_PATH_ENV);
    if (!path || !*path) path = KFLC_BUILTIN_PATH_DEFAULT;

    DIR *d = opendir(path);
    if (!d) {
        /* Silent miss is fine: hosts often have no install, embedders
         * may rely on the hardcoded registrar. */
        return;
    }

    /* Collect *.kflbi entries, sort alphabetically. */
    enum { MAX_ENTRIES = 256 };
    char *entries[MAX_ENTRIES];
    int   n_entries = 0;
    struct dirent *de;
    while ((de = readdir(d)) != NULL && n_entries < MAX_ENTRIES) {
        const char *n = de->d_name;
        size_t nl = strlen(n);
        size_t sl = strlen(KFLC_KFLBI_SUFFIX);
        if (nl <= sl) continue;
        if (strcmp(n + nl - sl, KFLC_KFLBI_SUFFIX) != 0) continue;
        size_t full_len = strlen(path) + 1 + nl + 1;
        char *full = (char *)malloc(full_len);
        if (!full) continue;
        snprintf(full, full_len, "%s/%s", path, n);
        entries[n_entries++] = full;
    }
    closedir(d);
    qsort(entries, (size_t)n_entries, sizeof entries[0], strptr_cmp_);

    for (int i = 0; i < n_entries; i++) {
        (void)parse_manifest_(entries[i]);
        free(entries[i]);
    }
}
