# Vendored: jacobwilliams/quadpack

Upstream: https://github.com/jacobwilliams/quadpack
Pinned commit: 702abfd5f0acbdb51439695334347a4b3c0dc87a (master, 2024-01-27)
Vendored: 2026-05-22
License: BSD-3-Clause (see LICENSE in this directory)

Modernized Fortran 2008 free-form refactor of the original QUADPACK
library (Piessens, de Doncker, Kahaner 1983). Single self-contained
module, no SLATEC / LINPACK / BLAS dependencies, no GOTOs, includes
the October 2021 SciPy bug fixes.

K26 vendors **double precision only** — the upstream provides a
precision-dispatch scheme (real32 / real64 / real128) but K26's
deterministic-binary contract is double precision throughout.

## What's vendored

- `quadpack_generic.F90` — the monolithic module body (8559 LOC); all
  QUADPACK procedures (DQAG, DQAGS, DQAGI, DQAGP, DQAWS, DQAWO, DQAWC,
  DQAWF, plus Gauss-Kronrod kernels, Chebyshev moments, etc.). Wrapped
  in `#ifndef MOD_INCLUDE` so it can be either compiled as its own
  module or `#include`d by a precision-specific wrapper.
- `quadpack_double.F90` — 7-line wrapper that sets `wp => real64`,
  defines `MOD_INCLUDE=1`, and includes `quadpack_generic.F90` to
  produce module `quadpack_double`. K26's ISO_C_BINDING wrapper
  (`src/k26astro_quad_iface.f90`) uses this module.

The upstream's other precision wrappers (`quadpack_single.F90`,
`quadpack_quad.F90`) and the umbrella module (`quadpack.F90`) are NOT
vendored — K26 doesn't expose single or quad precision quadrature.

## What's NOT vendored

- The FPM manifest (`fpm.toml`) — K26 uses its own Makefile pattern
  per the Fortran-link pattern doc.
- The upstream test suite — K26 ships its own
  `tests/test_quad_smoke.c` that exercises the C ABI surface.
- The Ford documentation source — public docs at jacobwilliams's site
  are authoritative; K26 doesn't regenerate them.

## How to update

1. Clone upstream at a new pinned tag/commit.
2. Diff `src/quadpack_generic.F90` + `src/quadpack_double.F90` against
   the in-tree copies.
3. Apply non-K26-breaking changes; re-run `make -C libk26astro_quad
   test` to confirm the C ABI surface still binds.
4. Update this file's pinned-commit + vendored-date lines.

## License compliance

BSD-3-Clause requires the LICENSE file to ship alongside the source.
`Makefile install` copies it into
`/usr/share/licenses/libk26astro_quad/quadpack/LICENSE`.
