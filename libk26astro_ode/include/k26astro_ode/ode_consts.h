/* ode_consts.h — DLSODA tolerance defaults + workspace sizing as
 * hex-literal IEEE-754 doubles for cross-platform determinism.
 *
 * References:
 *   - Petzold 1983 (LSODA original paper, SIAM J. Sci. Stat. Comput.)
 *   - Hindmarsh, ODEPACK user manual (LLNL 1980)
 *   - jacobwilliams/odepack modernized fork, commit 94bbe0d5 (2023);
 *     see src/upstream/odepack/UPSTREAM.md */
#ifndef K26ASTRO_ODE_CONSTS_H
#define K26ASTRO_ODE_CONSTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline double k26astro_ode_hex_d_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* LSODA conventional defaults: rtol = 1e-6, atol = 1e-9. Slightly
 * tighter than LSODA's internal "low precision" tier; appropriate for
 * physics integration where we typically want < 1 ppm error.
 *
 * Hex: 1e-6 = 0x3EB0C6F7A0B5ED8D
 * Hex: 1e-9 = 0x3E112E0BE826D695
 *
 * Tighter problems (Van der Pol stiff μ=1000, Robertson) usually need
 * rtol = 1e-8 / atol = 1e-12; loose problems (low-stiff harmonic)
 * tolerate rtol = 1e-3. K26 doesn't define those — pass them as
 * literals at the call site. */
#define K26ASTRO_ODE_RTOL_DEFAULT_BITS  0x3EB0C6F7A0B5ED8DULL
#define K26ASTRO_ODE_ATOL_DEFAULT_BITS  0x3E112E0BE826D695ULL

#define K26ASTRO_ODE_RTOL_DEFAULT  k26astro_ode_hex_d_(K26ASTRO_ODE_RTOL_DEFAULT_BITS)
#define K26ASTRO_ODE_ATOL_DEFAULT  k26astro_ode_hex_d_(K26ASTRO_ODE_ATOL_DEFAULT_BITS)

/* LSODA workspace sizing for full-Jacobian (JT=2, internally
 * generated). Per the DLSODA manual:
 *   LRW = 22 + NEQ * MAX(16, NEQ + 9)
 *   LIW = 20 + NEQ
 * K26 sizes workspace dynamically at call time based on `n`; see
 * src/k26astro_ode_iface.f90.
 *
 * For the canned opaque-handle surface (NEQ <= K26ASTRO_ODE_MAX_NEQ_HANDLE
 * = 16), the worst case is NEQ=16:
 *   LRW = 22 + 16 * MAX(16, 25) = 22 + 16*25 = 422
 *   LIW = 20 + 16 = 36
 * These bounds are baked into the Fortran iface's stack arrays so the
 * handle surface allocator-free in the integration path. */
#define K26ASTRO_ODE_LRW_HANDLE_MAX  422
#define K26ASTRO_ODE_LIW_HANDLE_MAX   36

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_ODE_CONSTS_H */
