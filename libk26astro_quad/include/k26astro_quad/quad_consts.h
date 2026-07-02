/* quad_consts.h — adaptive-quadrature tolerance defaults as hex-literal
 * IEEE-754 doubles for cross-platform determinism.
 *
 * References:
 *   - Piessens, de Doncker, Kahaner: QUADPACK manual, 1983
 *   - jacobwilliams/quadpack modernized fork, commit 702abfd5
 *     (2024-01-27); see src/upstream/quadpack/UPSTREAM.md */
#ifndef K26ASTRO_QUAD_CONSTS_H
#define K26ASTRO_QUAD_CONSTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Hex-literal double constructor (avoids strtod rounding noise across
 * libc implementations). */
static inline double k26astro_quad_hex_d_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* QUADPACK's classic default absolute tolerance is approximately
 * 1.49e-8 ≈ sqrt(epsilon_double). Same value used by SciPy's
 * scipy.integrate.quad. Hex: 0x3E45798EE2308C3A decodes to
 * 1.4901161193847656e-08 (exactly 2**-26, the sqrt of
 * machine epsilon 2**-52). */
#define K26ASTRO_QUAD_EPSABS_DEFAULT_BITS 0x3E45798EE2308C3AULL
#define K26ASTRO_QUAD_EPSREL_DEFAULT_BITS 0x3E45798EE2308C3AULL

#define K26ASTRO_QUAD_EPSABS_DEFAULT \
    k26astro_quad_hex_d_(K26ASTRO_QUAD_EPSABS_DEFAULT_BITS)
#define K26ASTRO_QUAD_EPSREL_DEFAULT \
    k26astro_quad_hex_d_(K26ASTRO_QUAD_EPSREL_DEFAULT_BITS)

/* Workspace sizing — these must match the Limit / Lenw passed by the
 * Fortran wrapper. K26's wrapper passes Limit=500, Lenw=2000 which
 * is sufficient for most well-behaved integrands. Increase here +
 * in src/k26astro_quad_iface.f90 in lock-step if you hit the
 * K26ASTRO_QUAD_E_LIMIT_REACHED error frequently. */
#define K26ASTRO_QUAD_WORKSPACE_LIMIT     500
#define K26ASTRO_QUAD_WORKSPACE_LENW      (4 * K26ASTRO_QUAD_WORKSPACE_LIMIT)

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_QUAD_CONSTS_H */
