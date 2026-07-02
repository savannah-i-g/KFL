# Vendored: jacobwilliams/odepack (M_odepack)

Upstream: https://github.com/jacobwilliams/odepack
Pinned commit: 94bbe0d5a9dd007dc9e07788ec9c8316417619c8 (master, 2023-04-17)
Vendored: 2026-05-22
License: Public Domain (LLNL release; see LICENSE)
Origin: ODEPACK (Hindmarsh 1983, LLNL); the jacobwilliams fork is a
        modern Fortran refactor of the original Netlib F77 tree.

ODEPACK is a collection of Fortran solvers for the initial-value
problem y' = f(t, y). K26 uses **DLSODA** specifically — the
stiff/non-stiff automatic-switching variant. The other solvers
(DLSODE, DLSODAR, DLSODES, DLSODI, DLSODIS, DLSODKR, DLSODPK,
DLSOIBT) are vendored but not exposed in the C ABI for now; they
can be wrapped later if a consumer drives one of them.

## What's vendored

The entire `src/` tree from upstream:

- `M_odepack.f90` — the wrapper module (~378 LOC at top level); uses
  `#include` to pull in all subroutines from `.inc` files.
- `M_da1/*.inc` — auxiliary routines (DEWSET, DINTDY, DSTODE, DPREP,
  DPREPJ, etc.). ~50 files.
- `M_main/*.inc` — the LSODA/LSODE/etc. entry points. 9 files
  (dlsoda.inc ~2 kLOC; full set ~22 kLOC).
- `M_matrix/*.inc` — BLAS-1 / LINPACK helpers (DAXPY, DDOT, DCOPY,
  DGEFA, etc.) and XERRWD error reporting. ~10 files.

Total: ~35 kLOC vendored Fortran. K26 builds it as a single
compilation unit (M_odepack.f90 with all the #includes resolving
relative to its directory).

## What's NOT vendored

- The fpm.toml manifest, upstream Makefile, examples (`example/`),
  and tests (`test/`) — K26 uses its own Makefile pattern per
  the Fortran-link pattern doc and its own smoke tests
  (`tests/test_ode_smoke.c`).
- The Ford documentation source and `ford.md` — public docs live at
  the upstream site.

## Building

K26's `libk26astro_ode/Makefile` compiles `M_odepack.f90` with
`-cpp` (so the `#include` directives resolve via the C preprocessor;
`.inc` files are searched relative to `M_odepack.f90`'s directory).
Required flags (in addition to the K26 standard FFLAGS):

- `-fallow-argument-mismatch` — upstream passes arrays where some
  callees expect scalars (legacy LSODA pattern); gfortran 14+ would
  hard-error without this flag.
- `-cpp` — enables the C preprocessor for #include directives.

## Public-domain certification

The LSODA family is in the public domain (LLNL release; see LICENSE,
"The ODEPACK package has been declared to be in the Public Domain.").
No license-compatibility constraints downstream.

## How to update

1. Clone upstream at a new pinned commit.
2. Diff `src/M_odepack.f90` and `src/M_*/*.inc` against in-tree.
3. Apply non-K26-breaking changes; rebuild + re-run
   `make -C libk26astro_ode test`.
4. Update this file's pinned-commit + vendored-date lines.

If the upstream restructures (e.g., adopts modules per entry point),
the K26-side Makefile may need to track those changes; see the
quadpack pattern in the Fortran-link pattern doc
for the build-rule template.
