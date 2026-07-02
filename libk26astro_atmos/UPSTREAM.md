# libk26astro_atmos — vendored upstream

## NRLMSISE-00(NRL FORTRAN direct vendoring)

The K26 atmos library vendors the original NRL FORTRAN distribution
of NRLMSISE-00 as its high-fidelity thermospheric model. The v0.4
ship landed the API surface (`k26astro_atmos_density_nrlmsise`) and
the wrapper (`src/atmos_nrlmsise.c`); the wrapper currently returns
`K26ASTRO_ATMOS_E_NOT_IMPLEMENTED` until the vendored upstream and
ISO_C_BINDING shim land.

### License gate

NRLMSISE-00, Picone et al. 2002, US Naval Research Laboratory.
US-government work — public domain. Passes the
Theme A license gate (BSD-2/3, MIT, Apache-2, CC0, public-domain).

The v0.4 UPSTREAM.md originally targeted Dominik Brodowski's C port,
but the v0.5 license eyes-on review (see memory
`project_nrlmsise_license_finding.md`) found the brodo repository
at `git.linta.de/~brodo/nrlmsise-00.git` lacks an explicit license
grant — only an FSF warranty disclaimer is present, no LICENSE /
COPYING file. Without an explicit BSD / MIT / PD / Apache grant the
derivative C cannot pass Theme A. The original NRL FORTRAN has no
such ambiguity, and the Fortran-link pattern
(QUADPACK / ODEPACK / IGRF — see
`feedback_fortran_lib_pattern.md`) is the natural fit.

### Architecture (Shape A — pure-C ABI via ISO_C_BINDING)

```
KFL surface
    │
    ▼
k26astro_atmos_density_nrlmsise (C, public API in include/)
    │
    ▼
k26astro_atmos_nrlmsise_call    (C, bind(C) entry in iface .f90)
    │
    ▼
GTD7 / GTD7D                    (NRL FORTRAN, vendored under src/upstream/)
```

`src/k26astro_atmos_nrlmsise_iface.f90` provides the ISO_C_BINDING
shim with `bind(C)` entries; `src/atmos_nrlmsise.c` becomes a thin
unit-conversion wrapper that calls those entries.

### Vendoring procedure (executed v0.5)

1. **Fetched NRL FORTRAN distribution** —
   `map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSISE-00/` via manual
   download (path requires browser; curl times out on automated
   fetch from this network). Single Fortran file
   `NRLMSISE-00.FOR` (113 KB), plus a journal-supplementary
   README (`NRLMSISE-00_2002JA009430-readme.txt`, 3.8 KB).

   **SHA-256 of the original combined file** (before split):
   `cce0420e90781c256bc6705c4cc8056b054812d6308adccb7081bf09af0d44cb`

2. **License audit** — Picone, Hedin, Drob (NRL Hulburt Center)
   + Aikin (NASA GSFC). All federal-employee authors per the
   journal-supplementary README — US-government work, PD by
   17 USC §105. No "FOR OFFICIAL USE ONLY" or distribution-
   restricted markers present in the source headers or README.

3. **Split the upstream into model + test driver** — the original
   `NRLMSISE-00.FOR` (2552 lines) bundled both the model
   subroutines/BLOCK DATA (lines 1-2435) AND a standalone test
   driver program (lines 2438-2551). The test driver has no
   `PROGRAM` statement, so in FORTRAN 77 it becomes the main entry
   point — which conflicts with linking the .a into anything that
   has its own `main()`. The fix: a clean split into

     - `NRLMSISE-00_model.FOR`       (109 KB; model only, ships in .a)
     - `NRLMSISE-00_test_driver.FOR` (4.2 KB; standalone driver)

   Both halves are byte-identical to the original (sed-extracted by
   line ranges, no edits). SHA-256s:

     - `NRLMSISE-00_model.FOR`:
       `3683831dae2496acb51013664763a0afde6d37eb07080ab8ef4a26cd95b931c4`
     - `NRLMSISE-00_test_driver.FOR`:
       `0e2fcba39e8359c8cc717ea198b5af382aadf92aea37d1bb4b79eae71ae0ff15`
     - `README.txt` (journal supp, brought along as documentation):
       `3a6e3d6ff89985d13fdbff4c348d9161bf10f50fcbba8f6d09b501b2625058dc`

4. **ISO_C_BINDING shim** at
   `libk26astro_atmos/src/k26astro_atmos_nrlmsise_iface.f90`. Three
   bind(C) entries: `k26astro_atmos_nrlmsise_init_call` (sets METERS
   = .TRUE. once), `k26astro_atmos_nrlmsise_gtd7_call` (the standard
   density+temperature call), and `k26astro_atmos_nrlmsise_gtd7d_call`
   (the drag-effective-density variant including anomalous oxygen).
   Mirrors the libk26astro_quad pattern but simpler — no integrand
   callback machinery.

5. **Makefile**:
     - `FC = gfortran` (force `=`, not `?=` — make's f77 default
       would lose).
     - `FFLAGS += -fdefault-real-8 -fdefault-double-8` — the NRL
       FORTRAN has no IMPLICIT or REAL*8 declarations; promoting
       default REAL to 8 bytes matches c_double in the ISO_C_BINDING
       shim. Standard convention for NRL FORTRAN codes.
     - `UPSTREAM_FFLAGS += -fdec` — enables the DEC extensions
       bundle, including the legacy Hollerith-to-INTEGER DATA
       initialiser pattern that BLOCK DATA GTD7BK uses for
       ISDATE/ISTIME/NAME at lines 1670-1671 of the upstream.
     - `UPSTREAM_FFLAGS += -std=legacy -ffixed-form -Wno-*` —
       silence the harmless conversion/range/argument-mismatch
       warnings the upstream emits under modern gfortran.
     - A consumer that links this library adds `-lgfortran` so the
       Fortran runtime resolves.

6. **Wrapper** at `src/atmos_nrlmsise.c` — one-shot METERS init then
   plain GTD7 call. SI output (kg/m³ + /m³) flows back into the
   K26AstroAtmosNrlmsiseDensities struct 1:1; no unit conversion
   needed because METERS=true handles it inside the upstream.

7. **Test gate** at `tests/test_nrlmsise.c` — three reference points
   at 200/400/1000 km mid-lat-quiet (replicating the NRL test driver's
   case-1 input shape). Reference values captured from a one-time
   run of the wrapper itself (since the wrapper is just an ABI shim,
   K26-wrapper output = NRL-FORTRAN output) and pinned to 17-digit
   `%.17g` doubles per `feedback_no_magic_numbers_future_updates.md`.
   `CAPTURE_REFS=1` env var dumps fresh values if the build flags
   intentionally change (e.g. losing -fdefault-real-8).

8. **Cross-check (informational, future)**: the v0.2 barometric
   profile from `k26astro_atmos_density_at` should agree with
   NRLMSISE-00 to within ~30% at 0-50 km altitude. Not gated.

### Build summary (verified 2026-05-23)

```
$ make -C libk26astro_atmos clean && make -C libk26astro_atmos
$ make -C libk26astro_atmos test
test_inscatter: OK
test_nrlmsise: OK (NRL FORTRAN wrapper agrees with pinned reference
                   values at 200/400/1000 km)
test_refraction: OK
```

### Why NOT brodo C

See `project_nrlmsise_license_finding.md`. Short version: brodo's
`DOCUMENTATION` file carries an FSF warranty disclaimer
("WITHOUT ANY WARRANTY") without naming a license. No LICENSE,
COPYING, or README files in the repo. The disclaimer alone does
not constitute a permission grant under default copyright law.
The provenance policy requires an explicit BSD / MIT / Apache / CC0
or PD grant; ambiguous-license vendoring is out of scope for the
project.
