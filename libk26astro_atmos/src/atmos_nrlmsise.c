/* atmos_nrlmsise.c — NRLMSISE-00 wrapper.
 *
 * Wraps the (vendored) NRL FORTRAN distribution of NRLMSISE-00 at
 * src/upstream/nrlmsise00/NRLMSISE-00.FOR (Picone, Hedin, Drob,
 * Aikin 2002; US-government public-domain). The vendoring procedure
 * + license rationale (we use the NRL FORTRAN directly rather than
 * brodo's C port, which lacks an explicit license grant) is at
 * libk26astro_atmos/UPSTREAM.md.
 *
 * The Fortran side exposes three bind(C) entries via the
 * `k26astro_atmos_nrlmsise_iface` module in
 * src/k26astro_atmos_nrlmsise_iface.f90:
 *
 *   k26astro_atmos_nrlmsise_init_call()
 *      Sets the upstream METERS global to .TRUE. so subsequent
 *      density values are in SI (per-m³ + kg/m³) rather than CGS.
 *
 *   k26astro_atmos_nrlmsise_gtd7_call(iyd, sec, alt, glat, glong,
 *      stl, f107a, f107, ap[7], mass, out_d[9], out_t[2])
 *      Forwards to the upstream GTD7 subroutine. Returns the 9-
 *      element density array + 2-element temperature array.
 *
 *   k26astro_atmos_nrlmsise_gtd7d_call(... same signature)
 *      GTD7D variant (anomalous-O included in total-mass output).
 *      Not currently exposed at the K26 API surface; reserved for
 *      future drag-specific paths.
 *
 * Thread safety: the upstream uses COMMON blocks (CSW, DATIME,
 * METERS state) for switch tracking. These are process-global; a
 * concurrent caller can read while another writes. Callers that
 * need parallel NRLMSISE evaluations must serialise externally.
 * For typical K26 use (single propagator + single observer per
 * frame) the single-threaded contract is fine.
 *
 * Unit conversion:
 *   Upstream returns (with METERS=true):
 *     D(1..5,7..9) — number densities (m^-3)
 *     D(6)         — total mass density (kg/m³)
 *     T(1..2)      — temperatures (K)
 *   K26 K26AstroAtmosNrlmsiseDensities surface is already in SI;
 *   no further conversion needed once METERS is set.
 */

#include "k26astro_atmos/atmos.h"

#include <stddef.h>

/* Fortran bind(C) entries — declared in
 * src/k26astro_atmos_nrlmsise_iface.f90. */
extern void k26astro_atmos_nrlmsise_init_call(void);
extern void k26astro_atmos_nrlmsise_gtd7_call(
    int iyd,
    double sec, double alt, double glat, double glong, double stl,
    double f107a, double f107,
    const double *ap,         /* 7 doubles */
    int mass,
    double *out_d,            /* 9 doubles */
    double *out_t);           /* 2 doubles */

/* Pack (year, day_of_year) into the upstream's YYDDD format. The
 * upstream documents that year is ignored by the current model, so
 * we don't bother validating the year value — only the day-of-year
 * matters for the diurnal phase. */
static int pack_iyd_(int year, int doy)
{
    int yy = year % 100;
    if (yy < 0) yy += 100;
    return yy * 1000 + doy;
}

/* One-shot METERS=true init. Idempotent; the upstream's METERS
 * subroutine is safe to call repeatedly. We gate behind a static
 * flag to skip the bind(C) call on hot paths. */
static int g_init_done_ = 0;
static void ensure_init_(void)
{
    if (g_init_done_) return;
    k26astro_atmos_nrlmsise_init_call();
    g_init_done_ = 1;
}

int k26astro_atmos_density_nrlmsise(
    const K26AstroAtmosNrlmsiseInputs *in,
    K26AstroAtmosNrlmsiseDensities    *out)
{
    if (!in || !out) return K26ASTRO_ATMOS_E_NOT_IMPLEMENTED;

    ensure_init_();

    double d_out[9] = {0};
    double t_out[2] = {0};
    int iyd = pack_iyd_(in->year, in->day_of_year);

    /* MASS=48: all species + temperature (the standard "give me
     * everything" mode). */
    k26astro_atmos_nrlmsise_gtd7_call(
        iyd, in->ut_seconds, in->alt_km,
        in->latitude_deg, in->longitude_deg,
        in->local_solar_time_hr,
        in->f107_81day_avg, in->f107_daily,
        in->ap, 48,
        d_out, t_out);

    /* Unpack — units are already SI (m^-3 + kg/m³) thanks to
     * METERS(.TRUE.) set by ensure_init_(). The upstream's D
     * ordering matches our surface 1:1: He, O, N2, O2, Ar,
     * total_mass, H, N, anomalous_O. Anomalous oxygen (D[8]) is
     * not surfaced at the K26 API; consumers wanting it should
     * call GTD7D directly (when that path is added). */
    out->n_he_m3     = d_out[0];
    out->n_o_m3      = d_out[1];
    out->n_n2_m3     = d_out[2];
    out->n_o2_m3     = d_out[3];
    out->n_ar_m3     = d_out[4];
    out->total_kg_m3 = d_out[5];
    out->n_h_m3      = d_out[6];
    out->n_n_m3      = d_out[7];
    /* d_out[8] = anomalous_O — not surfaced. */
    out->t_exo_k     = t_out[0];
    out->t_alt_k     = t_out[1];

    return K26ASTRO_ATMOS_OK;
}
