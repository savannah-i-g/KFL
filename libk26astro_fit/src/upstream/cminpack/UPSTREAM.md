# Vendored: devernay/cminpack

Upstream: https://github.com/devernay/cminpack
Pinned commit: 17dab75c6160d2ee42a3c95ea55e94738d7e559d (master, 2024-09-13)
Vendored: 2026-05-23
License: BSD-3-Clause (University of Chicago / Argonne; see LICENSE)
Origin: cminpack is a pure-C rewrite of MINPACK 1 (Argonne National
        Laboratory, More/Garbow/Hillstrom 1980). Provides the
        Levenberg-Marquardt nonlinear least-squares + Powell hybrid
        nonlinear-equation solvers.

K26 vendors the **double-precision Levenberg-Marquardt path only**,
which is sufficient for the v0.3 consumer set (orbit determination
from observation residuals, drag/SRP coefficient fitting,
spectroscopic calibration).

## What's vendored

- `cminpack.h` — public header (383 LOC)
- `cminpackP.h` — internal helpers + macros (70 LOC)
- `dpmpar.c` — machine precision constants
- `enorm.c` — Euclidean norm
- `fdjac2.c` — forward-difference Jacobian for lmdif
- `lmdif.c` — Levenberg-Marquardt (numerical Jacobian, full)
- `lmdif1.c` — convenience wrapper around lmdif
- `lmder.c` — Levenberg-Marquardt (user-supplied Jacobian, full)
- `lmder1.c` — convenience wrapper around lmder
- `lmpar.c` — Levenberg-Marquardt parameter solver
- `qrfac.c` — QR factorization
- `qrsolv.c` — QR triangular solve

Total: ~3.1 kLOC vendored C. K26 builds them as part of the
libk26astro_fit static archive.

## What's NOT vendored

- `hybrd*.c`, `hybrj*.c`, `hybrd1.c` — Powell hybrid nonlinear-
  equation solvers (square systems). Add when a consumer drives them.
- `lmstr*.c` — minpack's storage-saving lmder variant for very large
  m problems. K26's first surface handles modest m (< few thousand);
  add when scale demands.
- `chkder.c` — Jacobian-check helper. Useful for analytical-Jacobian
  development; add when lmder is exercised.
- `covar*.c` — covariance matrix from final Jacobian. Add when
  uncertainty propagation is needed.
- The `*_.c` Fortran-compatibility shims, BLAS/LAPACK bindings,
  CUDA variants, CMake/autoconf glue, tests, docs.

## Building

Pure C, no Fortran. Standard CFLAGS:

- `-std=c11 -O2 -ffp-contract=off -fexcess-precision=standard`
  (K26 determinism contract)
- Default macros — vendored sources compile against the default
  configuration (no `__cminpack_double__` define needed; the header
  picks `double` automatically when no precision macro is set).
- Include `-Isrc/upstream/cminpack` so the .c files find `cminpack.h`
  and `cminpackP.h`.

## License compliance

BSD-3-Clause (University of Chicago / Argonne). `make install` copies
LICENSE into `/usr/share/licenses/libk26astro_fit/cminpack/LICENSE`.

## How to update

1. Clone upstream at a new pinned tag/commit.
2. Diff the per-file vendored set above against in-tree copies.
3. Apply non-K26-breaking changes; rebuild + re-run
   `make -C libk26astro_fit test`.
4. Update this file's pinned-commit + vendored-date lines.
