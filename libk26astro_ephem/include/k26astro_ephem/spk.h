/* k26astro_ephem/spk.h — NAIF SPK (Spacecraft and Planet Kernel)
 * binary reader, focused subset.
 *
 * Subset covered (v0.1):
 *   - DAF/SPK file architecture (NAIF Double-precision Array File)
 *   - Type 2 segments only: Chebyshev polynomials for Cartesian
 *     positions (the standard storage for DE-series planet
 *     ephemerides like DE441)
 *   - Little-endian files only (the LOCFMT identifier must read
 *     "LTL-IEEE"); modern x86_64 / aarch64 distributions are
 *     little-endian everywhere this matters
 *   - Reference frame: J2000/ICRF (frame id 1) — DE441's native
 *
 * Not implemented (deferred until a use case lands):
 *   - Type 3 (Chebyshev positions + velocities) — straightforward
 *     extension; same record shape, twice the coefficients
 *   - Type 8/9/13 (interpolation-based)
 *   - Big-endian SPK files (require byte-swap on every read)
 *   - Frame kernels (FK), planetary constants kernels (PCK), CK
 *     orientation kernels — these live in different file types and
 *     don't intersect with SPK reading
 *
 * Loading is mmap-based: pages are demand-faulted as queries hit
 * them, no copy to RAM, and closing the kernel is a single munmap.
 */
#ifndef K26ASTRO_EPHEM_SPK_H
#define K26ASTRO_EPHEM_SPK_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Opaque kernel handle --------------------------------------- */
typedef struct K26AstroSpk K26AstroSpk;

/* ---- Segment descriptor ---------------------------------------- */
/* One entry per Type-2 SPK segment found in the kernel. Public so
 * the higher-level query layer in ephem.c can iterate the catalogue
 * directly. */
typedef struct {
    int      target_body;     /* NAIF body id (e.g. 399 = Earth) */
    int      center_body;     /* NAIF id of the reference centre */
    int      frame_id;        /* 1 = J2000/ICRF (the only one we serve) */
    double   start_et;        /* TDB seconds past J2000 */
    double   end_et;
    /* Type-2 layout: each record covers `interval_seconds` of ET and
     * holds Chebyshev coefficients for (x, y, z). */
    double   init_et;         /* epoch of first record */
    double   interval_seconds;
    int      record_size_doubles;  /* RSIZE — total doubles per record */
    int      n_records;
    /* Coefficient layout per record:
     *   [0] = MID epoch (record-midpoint, TDB seconds past J2000)
     *   [1] = RADIUS (half-interval, == interval_seconds / 2)
     *   [2 .. 2+order]                       = x coeffs (order+1 of them)
     *   [2+order+1 .. 2+2*(order+1)]          = y coeffs
     *   [2+2*(order+1)+1 .. 2+3*(order+1)]    = z coeffs
     * Note: SPK stores in C0..C_n order, ready for Clenshaw without
     * a reverse. */
    int      coeffs_per_axis;     /* (RSIZE - 2) / 3 */
    /* Pointer into the mmap'd file at the first record. */
    const double *records;
} K26AstroSpkSegment;

/* ---- Load / close ----------------------------------------------- */
K26AstroSpk *k26astro_spk_open (const char *path);
void         k26astro_spk_close(K26AstroSpk *spk);

/* Total Type-2 segments parsed. */
int                       k26astro_spk_n_segments(const K26AstroSpk *spk);
const K26AstroSpkSegment *k26astro_spk_segment   (const K26AstroSpk *spk, int idx);

/* Find a segment by NAIF target body whose [start_et, end_et]
 * contains `et`. Returns NULL if no segment matches. Linear scan;
 * sufficient for the ~10-segment v0.1 inner-planet kernel. */
const K26AstroSpkSegment *
k26astro_spk_find_segment(const K26AstroSpk *spk, int target_body, double et);

/* ---- Position evaluation --------------------------------------- */
/* Evaluate (x, y, z) for `target_body` at TDB-seconds-past-J2000
 * `et`. Returns 0 on success, non-zero if the segment isn't
 * available (no kernel coverage). On success, out_xyz is written in
 * the kernel's native reference frame + length units (km for
 * DE-series). */
int k26astro_spk_pos(const K26AstroSpk *spk,
                     int target_body, double et,
                     double out_xyz[3]);

/* Evaluate position + velocity. Velocity is the analytic derivative
 * of the position Chebyshev series, in length-units / second. */
int k26astro_spk_pos_vel(const K26AstroSpk *spk,
                         int target_body, double et,
                         double out_xyz[3], double out_xyz_dot[3]);

/* ---- Synthetic-fixture writer (for tests) ---------------------- */
/* Write a minimal valid Type-2 SPK file. Used by the test suite to
 * generate fixtures without depending on a real DE441 kernel in
 * tree. Not intended for production ephemerides — analytic test
 * series only. Returns 0 on success, non-zero on I/O failure or
 * invalid input. */
typedef struct {
    int     target_body;
    int     center_body;
    double  start_et;
    double  end_et;
    double  interval_seconds;
    /* Caller supplies coefficient blocks for each interval, packed
     * as [MID, RADIUS, x_coeffs.., y_coeffs.., z_coeffs..]. The
     * helper writes them verbatim; the writer doesn't sample any
     * function on the caller's behalf. */
    const double *records;
    int     n_records;
    int     coeffs_per_axis;
} K26AstroSpkWriteSegment;

int k26astro_spk_write_synthetic(const char *path,
                                  const K26AstroSpkWriteSegment *segs,
                                  int n_segs);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_EPHEM_SPK_H */
