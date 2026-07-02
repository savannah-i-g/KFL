# Vendored: NOAA NCEI IGRF-14 reference

Upstream: http://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f
Pinned: as published 2024-11 by IAGA Working Group V-MOD (released by
        NOAA NCEI); IGRF-14 is the 14th generation, 13th revision.
        Coefficient table valid 1900.0 – 2030.0 (definitive
        1945.0 – 2020.0; predictive secular variation 2025.0 – 2030.0).
Vendored: 2026-05-23
License: Public Domain (US Government work; see LICENSE)

The Fortran reference is the canonical IGRF-14 implementation
maintained by BGS (Susan Macmillan, William Brown, Ciaran Beggan)
and released by NOAA NCEI.

## What's vendored

- `igrf14syn.f` — the synthesis subroutine plus two degree-minute
  conversion helpers (DMDDEC, DDECDM). 754 LOC, fixed-form F77.
  Stripped from the upstream `igrf14.f` which also contains a
  PROGRAM IGRF14 interactive driver — that PROGRAM block would
  conflict with `main` at link time, so we keep only the
  subroutines:

  - `igrf14syn(isv, date, itype, alt, colat, elong, x, y, z, f)` —
    the actual magnetic-field synthesis routine. F77, embeds the
    IGRF-14 coefficient table inline.
  - `DMDDEC(I, M, X)` — degrees+minutes → decimal degrees helper.
  - `DDECDM(X, I, M)` — inverse helper.

  K26 calls only `igrf14syn` from the ISO_C_BINDING wrapper
  (`src/k26astro_geomag_iface.f90`); DMDDEC/DDECDM are unused but
  kept for upstream-diff cleanliness.

## What's NOT vendored

- The PROGRAM IGRF14 interactive driver (lines 1–447 of upstream
  igrf14.f). Stripped because it provides `main` which would
  conflict with consumer apps' `main`.
- The optional `igrf14coeffs.txt` external file. The Fortran
  subroutine has the coefficient table embedded in DATA statements;
  no external file is needed at runtime.

## Building

F77 fixed-form (`.f` extension). Standard K26 Fortran flags + the
gfortran F77 dialect handling:

- Use `gfortran` with default `-std=gnu` (gfortran recognizes
  fixed-form from `.f` extension automatically).
- `-fno-automatic` makes local SAVE-style allocations explicit
  (matches legacy F77 expectations).
- `-fdefault-double-8 -fdefault-real-8` are NOT used — the upstream
  uses `IMPLICIT DOUBLE PRECISION (A-H,O-Z)` so reals default to
  double already.

## Coefficient validity range

- 1900.0 – 1945.0 : non-definitive (extrapolation from historical
  observations)
- 1945.0 – 2020.0 : definitive (5-year snapshots)
- 2020.0 – 2025.0 : definitive at 2020.0, predictive thereafter
- 2025.0 – 2030.0 : predictive secular variation only

Outside [1900.0, 2035.0] the subroutine prints a warning and
returns f = 1.0e8, x = y = z = 0. K26's wrapper translates this
into a non-zero return code so consumers can detect and handle
out-of-range queries.

## How to update (IGRF-15, ~2030)

1. Download new upstream from
   http://www.ngdc.noaa.gov/IAGA/vmod/igrfNN.f (NN = generation).
2. Strip PROGRAM IGRFNN block (replace `sed '448,$p'` cutoff with
   the new file's PROGRAM-end line).
3. Rename the call site in `src/k26astro_geomag_iface.f90` from
   `igrf14syn` to `igrfNNsyn`.
4. Rebuild and re-run smoke tests against the new reference points
   from the BGS validation tables.
5. Update this file's vendored-date + generation lines.
