#ifndef K26_STR_H
#define K26_STR_H

/* libk26util: tiny string helpers exposed to KFL's expression
 * sub-language.
 *
 * The kflc compiler's BUILTINS table (expr.c) maps a small set of
 * KFL-callable names (strlen / streq / starts_with / ends_with /
 * concat) onto these wrappers. Kept separate from <string.h> so the
 * KFL surface is stable across libc evolutions and the names are
 * unambiguously K26-flavoured in emitted C++ output.
 *
 * Each function tolerates NULL inputs by returning the "empty" answer
 * for its return type (len 0, equality false, etc.); KFL form-args
 * can be unset and we don't want a crash there. */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

size_t k26_str_len(const char *s);
int    k26_str_eq(const char *a, const char *b);
int    k26_str_starts_with(const char *s, const char *prefix);
int    k26_str_ends_with  (const char *s, const char *suffix);

/* Concatenate `a` and `b` into `out` (NUL-terminated, capped at
 * `cap`). Caller-owned buffer; the KFL `concat` builtin uses a
 * per-callsite `static char[256]` via a GCC statement-expression
 * wrapper at emit time. Returns the length written (excluding NUL),
 * or 0 if cap == 0. */
size_t k26_str_concat(char *out, size_t cap, const char *a, const char *b);

#ifdef __cplusplus
}
#endif

#endif /* K26_STR_H */
