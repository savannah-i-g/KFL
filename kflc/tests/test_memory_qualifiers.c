/* test_memory_qualifiers.c: memory-model enforcement gate.
 *
 * Exercises the compile-time enforcement of own / borrow / ptr / arena
 * qualifiers in kflc. Positive cases confirm valid surface emits the
 * expected C++ shape; negative cases confirm the enforcement fires
 * with the right error message.
 *
 * Pattern: write a small .kfl fixture to a tmpfile, run
 *   ./bin/kflc --emit <tmpfile>
 * capture stderr + stdout + exit code, assert on the captured state.
 *
 * Wire: see kflc/Makefile MEMORY_QUAL_TEST + test target.
 */
#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

static void write_fixture_(const char *path, const char *content)
{
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(2); }
    fputs(content, f);
    fclose(f);
}

/* Runs `kflc --check <fixture>` (silent if clean; errors to stderr).
 * Returns exit code; writes stderr text to *err_out (caller frees). */
static int run_check_(const char *fixture, char **err_out)
{
    char cmd[512];
    snprintf(cmd, sizeof cmd,
             "./bin/kflc --check %s 2>/tmp/kflc_check_err.log", fixture);
    int rc = system(cmd);
    FILE *f = fopen("/tmp/kflc_check_err.log", "rb");
    if (!f) { *err_out = strdup(""); return WEXITSTATUS(rc); }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    if (sz < 0) sz = 0;
    fseek(f, 0, SEEK_SET);
    *err_out = (char *)malloc((size_t)sz + 1);
    if (sz > 0) (void)!fread(*err_out, 1, (size_t)sz, f);
    (*err_out)[sz] = '\0';
    fclose(f);
    return WEXITSTATUS(rc);
}

/* Runs `kflc --emit <fixture>` and returns the emitted C++ text
 * (stdout). Stderr is suppressed. Caller frees. */
static char *run_emit_(const char *fixture)
{
    char cmd[512];
    snprintf(cmd, sizeof cmd,
             "./bin/kflc --emit %s 2>/dev/null > /tmp/kflc_emit.cc",
             fixture);
    (void)system(cmd);
    FILE *f = fopen("/tmp/kflc_emit.cc", "rb");
    if (!f) return strdup("");
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    if (sz < 0) sz = 0;
    fseek(f, 0, SEEK_SET);
    char *buf = (char *)malloc((size_t)sz + 1);
    if (sz > 0) (void)!fread(buf, 1, (size_t)sz, f);
    buf[sz] = '\0';
    fclose(f);
    return buf;
}

static const char *FIX_POS_MOVE =
    "form POSMOVE\n"
    "title \"positive move\"\n"
    "size 320x200\n"
    "cfg \"pos-move\"\n"
    "fn double demo(own body x)\n"
    "    let y: own body = move(x)\n"
    "    return 1.0\n"
    "end\n"
    "end\n";

static const char *FIX_POS_BORROW =
    "form POSBORROW\n"
    "title \"positive borrow\"\n"
    "size 320x200\n"
    "cfg \"pos-borrow\"\n"
    "fn double demo(borrow body x)\n"
    "    let v: borrow body = x\n"
    "    return 2.0\n"
    "end\n"
    "end\n";

static const char *FIX_POS_PTR =
    "form POSPTR\n"
    "title \"positive ptr\"\n"
    "size 320x200\n"
    "cfg \"pos-ptr\"\n"
    "fn double demo(ptr body x)\n"
    "    let r: ptr body = x\n"
    "    return 3.0\n"
    "end\n"
    "end\n";

static const char *FIX_POS_ARENA =
    "form POSARENA\n"
    "title \"positive arena\"\n"
    "size 320x200\n"
    "cfg \"pos-arena\"\n"
    "arena scratch capacity 4096\n"
    "arena manual_buf capacity 4096 reset_mode manual\n"
    "fn double demo()\n"
    "    allocator = scratch\n"
    "    return 4.0\n"
    "end\n"
    "fn double demo2()\n"
    "    allocator = manual_buf\n"
    "    return 5.0\n"
    "end\n"
    "end\n";

static const char *FIX_NEG_BARE_OWN =
    "form NEGBARE\n"
    "title \"neg bare\"\n"
    "size 320x200\n"
    "cfg \"neg-bare\"\n"
    "fn double demo(own body x)\n"
    "    let y: own body = x\n"
    "    return 1.0\n"
    "end\n"
    "end\n";

static const char *FIX_NEG_DOUBLE_MOVE =
    "form NEGDM\n"
    "title \"neg dm\"\n"
    "size 320x200\n"
    "cfg \"neg-dm\"\n"
    "fn double demo(own body x)\n"
    "    let a: own body = move(x)\n"
    "    let b: own body = move(x)\n"
    "    return 1.0\n"
    "end\n"
    "end\n";

static const char *FIX_NEG_MOVE_NONOWN =
    "form NEGMNO\n"
    "title \"neg move-non-own\"\n"
    "size 320x200\n"
    "cfg \"neg-mno\"\n"
    "fn double demo(borrow body x)\n"
    "    let y: own body = move(x)\n"
    "    return 1.0\n"
    "end\n"
    "end\n";

int main(void)
{
    /* ---- POSITIVE cases ------------------------------------------- */

    write_fixture_("/tmp/kflc_pos_move.kfl", FIX_POS_MOVE);
    char *err = NULL;
    int rc = run_check_("/tmp/kflc_pos_move.kfl", &err);
    assert(rc == 0);
    free(err);
    char *cxx = run_emit_("/tmp/kflc_pos_move.kfl");
    /* move() emits the GCC stmt-expr with null-out of source. */
    assert(strstr(cxx,
        "({ K26AstroBody * _kfl_mv = x; x = NULL; _kfl_mv; })") != NULL);
    free(cxx);

    write_fixture_("/tmp/kflc_pos_borrow.kfl", FIX_POS_BORROW);
    rc = run_check_("/tmp/kflc_pos_borrow.kfl", &err);
    assert(rc == 0);
    free(err);
    cxx = run_emit_("/tmp/kflc_pos_borrow.kfl");
    /* borrow emits as const T *. */
    assert(strstr(cxx, "const K26AstroBody * v = x;") != NULL);
    free(cxx);

    write_fixture_("/tmp/kflc_pos_ptr.kfl", FIX_POS_PTR);
    rc = run_check_("/tmp/kflc_pos_ptr.kfl", &err);
    assert(rc == 0);
    free(err);
    cxx = run_emit_("/tmp/kflc_pos_ptr.kfl");
    /* ptr emits as raw T *. */
    assert(strstr(cxx, "K26AstroBody * r = x;") != NULL);
    free(cxx);

    write_fixture_("/tmp/kflc_pos_arena.kfl", FIX_POS_ARENA);
    rc = run_check_("/tmp/kflc_pos_arena.kfl", &err);
    assert(rc == 0);
    free(err);
    cxx = run_emit_("/tmp/kflc_pos_arena.kfl");
    /* Default (fn_exit) emits with cleanup attribute. */
    assert(strstr(cxx,
        "__attribute__((cleanup(_kfl_arena_reset_cleanup_))) "
        "= _kfl_arena_scratch") != NULL);
    /* `reset_mode manual` skips the cleanup attribute. */
    assert(strstr(cxx,
        "K26KflArena *_kfl_active_arena = _kfl_arena_manual_buf") != NULL);
    free(cxx);

    /* ---- NEGATIVE cases ------------------------------------------- */

    write_fixture_("/tmp/kflc_neg_bare.kfl", FIX_NEG_BARE_OWN);
    rc = run_check_("/tmp/kflc_neg_bare.kfl", &err);
    assert(rc != 0);
    assert(strstr(err, "RHS is a bare `own` binding") != NULL);
    assert(strstr(err, "Wrap with `move(") != NULL);
    free(err);

    write_fixture_("/tmp/kflc_neg_dm.kfl", FIX_NEG_DOUBLE_MOVE);
    rc = run_check_("/tmp/kflc_neg_dm.kfl", &err);
    assert(rc != 0);
    assert(strstr(err, "already moved-from") != NULL);
    free(err);

    write_fixture_("/tmp/kflc_neg_mno.kfl", FIX_NEG_MOVE_NONOWN);
    rc = run_check_("/tmp/kflc_neg_mno.kfl", &err);
    assert(rc != 0);
    assert(strstr(err, "must be `own` or `ptr` qualified") != NULL);
    free(err);

    printf("test_memory_qualifiers: OK "
           "(4 positive + 3 negative cases, all qualifier checks fire)\n");
    return 0;
}
