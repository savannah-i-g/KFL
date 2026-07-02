#ifndef KFLC_RUNTIME_H
#define KFLC_RUNTIME_H

/* kflc_runtime.h , runtime support for the KFL memory model.
 *
 * Linked into kflc-emitted C++ programs when the source uses any of
 * the memory-model keywords (`own` / `borrow` / `ptr` / `arena` /
 * `allocator`). The runtime is small (~150 LOC) and self-contained;
 * built into libkflc.a and pulled in by the generated binary's link
 * line via `-lkflc`.
 *
 * Three concerns live here:
 *
 *   1. Arena bump allocator. KFL form-level `arena <name> capacity
 *      <bytes>` declarations emit a `K26KflArena *_kfl_arena_<name>`
 *      lifetime-of-form handle; `allocator = <name>` inside a fn
 *      body opts subsequent allocations into the arena via
 *      k26kfl_arena_alloc instead of the default malloc. Arena
 *      lifetime tracks the form (or longer; manual reset_mode);
 *      individual allocations cannot be freed in isolation.
 *
 *   2. Move-semantics primitive. `own` bindings emit assignments
 *      via k26kfl_move() which null-out the source after the move.
 *      Compile-time use-after-move detection lives in the kflc
 *      emitter; this primitive is the runtime substrate.
 *
 *   3. Borrow / ptr have no runtime; they are compile-time
 *      contracts enforced by the kflc emitter (borrow scope check,
 *      ptr null-check obligation on the receiver). They appear
 *      in this header only as documentation.
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Arena bump allocator -------------------------------------- */

typedef struct K26KflArena K26KflArena;

/* Create an arena with the given capacity in bytes. Returns NULL on
 * allocation failure or capacity == 0. The arena's internal buffer
 * is allocated once at create time; subsequent k26kfl_arena_alloc
 * calls bump a pointer within that buffer. */
K26KflArena *k26kfl_arena_create(size_t capacity_bytes);

/* Release the arena and all allocations within it. After this call
 * any pointer returned from k26kfl_arena_alloc(a, ...) is invalid. */
void k26kfl_arena_destroy(K26KflArena *a);

/* Allocate `bytes` bytes from the arena. Returns NULL if the
 * remaining capacity is insufficient. Returned pointer is aligned
 * to alignof(max_align_t); valid until k26kfl_arena_reset() or
 * k26kfl_arena_destroy() is called. */
void *k26kfl_arena_alloc(K26KflArena *a, size_t bytes);

/* Reset the arena's bump pointer to 0, invalidating every allocation
 * within it (without releasing the underlying buffer). Used by the
 * `reset_mode fn` arena binding which resets on fn entry. */
void k26kfl_arena_reset(K26KflArena *a);

/* Return the number of bytes currently used (bump pointer position). */
size_t k26kfl_arena_used(const K26KflArena *a);

/* Return the arena's total capacity (set at create time). */
size_t k26kfl_arena_capacity(const K26KflArena *a);

/* ---- Move semantics primitive ---------------------------------- */

/* Move the void * value pointed to by `src` into `*dst`, nulling out
 * `*src` after the move. The kflc emitter generates this call for
 * `own`-qualified opaque assignments. The primitive itself does not
 * encode ownership — that contract is enforced at compile time by
 * the emitter's use-after-move check. */
void k26kfl_move(void **dst, void **src);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* KFLC_RUNTIME_H */
