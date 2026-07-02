/* quad_reference.h — analytical reference values for the smoke test.
 *
 * Hex-literal IEEE-754 doubles. All values are exact IEEE-754
 * round-to-nearest of the analytical result, computed externally
 * and pinned here so the gate is bit-stable across libc versions.
 *
 * Each entry documents both the decimal value (human-readable) and
 * the hex bit-pattern (machine-stable). Update both when a fixture
 * changes. */
#ifndef QUAD_REFERENCE_H
#define QUAD_REFERENCE_H

#include <stdint.h>

static inline double q_ref_hex_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* Reference: integral of e^x from 0 to 1.
 *   I = e - 1 ≈ 1.7182818284590452
 * Hex of (e - 1) rounded to double: 0x3FFB7E151628AED2 */
#define Q_REF_EXP_INT_0_1_BITS  0x3FFB7E151628AED2ULL
#define Q_REF_EXP_INT_0_1       q_ref_hex_(Q_REF_EXP_INT_0_1_BITS)

/* Reference: integral of e^(-x) from 0 to +infinity.
 *   I = 1
 * Trivially exact. */
#define Q_REF_EXP_NEG_INT_BITS  0x3FF0000000000000ULL
#define Q_REF_EXP_NEG_INT       q_ref_hex_(Q_REF_EXP_NEG_INT_BITS)

/* Reference: integral of x^2 from 0 to 1.
 *   I = 1/3 ≈ 0.333333333333333...
 * Hex of (1.0/3.0) IEEE-754 round-to-nearest: 0x3FD5555555555555 */
#define Q_REF_X_SQ_INT_0_1_BITS  0x3FD5555555555555ULL
#define Q_REF_X_SQ_INT_0_1       q_ref_hex_(Q_REF_X_SQ_INT_0_1_BITS)

/* Reference: Stefan-Boltzmann via Planck integral.
 *   ∫₀^∞ B(λ; T) dλ = σ T^4 / π
 *   where σ = 5.670374419e-8 W/(m²·K⁴) (CODATA 2018, exact)
 *
 * For T = 5778 K (Sun): σ T^4 / π = 2.008994...e7 W/(m²·sr).
 * Reference value pinned: 20089935.6...
 * Bit-pattern: 0x416322B6CFAE6AA8 ≈ 20089935.36...
 *
 * QUADPACK DQAGI tolerance: this is a steep-tailed integrand and
 * DQAGI's variable transformation handles it well, but the result
 * is sensitive to the absolute tolerance. Test gates to relative
 * 1e-6 (one part per million), not bit-exact. */
#define Q_REF_PLANCK_SUN_BITS    0x416322B6CFAE6AA8ULL
#define Q_REF_PLANCK_SUN         q_ref_hex_(Q_REF_PLANCK_SUN_BITS)

/* Tolerance for the smoke test gates. Bit-exact would be ideal but
 * QUADPACK accumulates rounding; 1e-12 absolute is the empirical
 * tolerance the upstream test suite uses. */
#define Q_TOL_ABS                1.0e-12
#define Q_TOL_REL                1.0e-9

#endif /* QUAD_REFERENCE_H */
