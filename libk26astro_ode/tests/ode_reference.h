/* ode_reference.h — analytical reference values for the smoke test.
 *
 * Hex-literal IEEE-754 doubles. References pinned from analytical
 * solutions or published canonical-benchmark results so the gate
 * is bit-stable.
 *
 * For numerical ODE solvers like LSODA the gate is necessarily
 * relative-tolerance (the solver itself returns approximate values);
 * absolute hex-comparison only valid for analytical-closed-form
 * problems integrated past one period. */
#ifndef ODE_REFERENCE_H
#define ODE_REFERENCE_H

#include <stdint.h>
#include <math.h>

static inline double o_ref_hex_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* Harmonic oscillator (omega=1): y = [x, v], y(0) = [1, 0].
 * Analytical: y(t) = [cos(t), -sin(t)].
 * At t = pi:   y = [-1, 0]
 * At t = 2pi:  y = [ 1, 0] (period)
 *
 * Hex(cos(pi)) = -1.0 = 0xBFF0000000000000
 * Hex(sin(pi)) ≈ 1.22464679914735e-16 — round-off, not exact zero. */
#define O_REF_HARMONIC_X_AT_PI  (-1.0)
#define O_REF_HARMONIC_V_AT_PI   (0.0)

/* Robertson chemistry benchmark — published canonical reference at
 * various t values (Hindmarsh & Petzold 1980). y(0) = [1, 0, 0],
 * default rate constants k1=0.04, k2=1e4, k3=3e7.
 *
 * At t = 4e-1:  y ≈ [0.985172, 3.38637e-5, 0.0147939]
 * At t = 4e+0:  y ≈ [0.905516, 2.24043e-5, 0.0944614]
 * Hex pinning here would invite tracking-error drift; the test gates
 * to relative 1e-4 against these published values. */
#define O_REF_ROBERTSON_T0_4_Y0  0.9851721
#define O_REF_ROBERTSON_T0_4_Y2  0.0147939

/* Tolerance constants. LSODA produces approximate solutions; tests
 * compare to relative 1e-6 against analytical values (default rtol
 * 1e-9 + atol 1e-12 gives sub-ppm error in our cases). */
#define O_TOL_REL_ANALYTICAL  1.0e-6
#define O_TOL_REL_BENCHMARK   1.0e-4

#endif /* ODE_REFERENCE_H */
