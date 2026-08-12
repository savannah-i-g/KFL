/* jammer.c — jammer opaque + radar-EW equation.
 *
 * Implements the standard self-protection jamming geometry:
 * jammer is co-located with the target (jammer host = victim
 * radar target). Stand-off jamming is handled by the same
 * J/S formula with separate jammer-to-radar and target-to-radar
 * ranges; v0.1 exposes self-protection only. */

#include "k26astro_softkill/jammer.h"

#include "k26astro_core/consts.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct K26AstroJammer {
    K26AstroPayload     base;
    K26AstroJammerMode mode;
    double             p_j_W;
    double             g_j_db;
    double             freq_Hz;
    double             bandwidth_Hz;
    double             snr_threshold;
};

static double db_to_linear_(double db) { return pow(10.0, db / 10.0); }

K26AstroJammer *k26astro_jammer_new(
    K26AstroJammerMode mode,
    double             p_j_W,
    double             g_j_db,
    double             freq_Hz,
    double             bandwidth_Hz,
    double             snr_threshold)
{
    if (mode != K26ASTRO_JAMMER_MODE_NOISE &&
        mode != K26ASTRO_JAMMER_MODE_COVER_PULSE &&
        mode != K26ASTRO_JAMMER_MODE_DECEPTION) return NULL;
    if (p_j_W <= 0.0 || freq_Hz <= 0.0 || bandwidth_Hz <= 0.0) return NULL;
    if (snr_threshold <= 0.0) return NULL;

    K26AstroJammer *j = (K26AstroJammer *)calloc(1, sizeof(*j));
    if (!j) return NULL;
    j->base.kind             = K26ASTRO_DEFENSE_KIND_JAMMER;
    j->base.owner            = NULL;
    j->base.owner_slot       = -1;
    j->base.generation       = 1;
    j->base.on_owner_destroy = NULL;
    j->mode                  = mode;
    j->p_j_W                 = p_j_W;
    j->g_j_db                = g_j_db;
    j->freq_Hz               = freq_Hz;
    j->bandwidth_Hz          = bandwidth_Hz;
    j->snr_threshold         = snr_threshold;
    return j;
}

void k26astro_jammer_destroy(K26AstroJammer *j)
{
    if (!j) return;
    k26astro_payload_unlink_(&j->base);
    free(j);
}

int k26astro_jammer_bind(K26AstroJammer *j, struct K26AstroVehicle *v)
{
    if (!j) return -1;
    return k26astro_payload_bind(&j->base, v);
}

struct K26AstroPayload *k26astro_jammer_as_payload(K26AstroJammer *j)
{
    return j ? &j->base : NULL;
}

K26AstroJammerMode k26astro_jammer_mode(const K26AstroJammer *j)
{
    return j ? j->mode : K26ASTRO_JAMMER_MODE_NOISE;
}

/* J/S = (P_j · G_j · 4π · R²) / (P_t · G_t · σ)
 *
 * Skolnik 2008 §15.3 self-protection geometry. Derived from
 *
 *   J = P_j · G_j · G_r · λ² / ((4π)² · R²)             one-way
 *   S = σ · P_t · G_t · G_r · λ² / ((4π)³ · R⁴)         two-way
 *   J/S = (P_j · G_j · 4π · R²) / (P_t · G_t · σ)
 *
 * The jamming power arrives at the victim receiver at 1/R²
 * (one-way) whereas the skin echo arrives at 1/R⁴ (two-way), so
 * J/S grows as R² — the jammer wins at long range. Larger target
 * RCS produces a stronger skin echo and reduces J/S; σ is in the
 * denominator.
 */
double k26astro_jammer_js_ratio(const K26AstroJammer *j,
                                 double p_r_W,
                                 double g_r_t_db,
                                 double rcs_target_m2,
                                 double range_m)
{
    if (!j || p_r_W <= 0.0 || range_m <= 0.0) return 0.0;
    if (rcs_target_m2 <= 0.0) return 0.0;

    const double G_j = db_to_linear_(j->g_j_db);
    const double G_t = db_to_linear_(g_r_t_db);
    const double numer = j->p_j_W * G_j * 4.0 * K26A_PI * range_m * range_m;
    const double denom = p_r_W * G_t * rcs_target_m2;
    if (denom <= 0.0) return 0.0;
    return numer / denom;
}

double k26astro_jammer_burn_through_range(const K26AstroJammer *j,
                                           double p_r_W,
                                           double g_r_t_db,
                                           double rcs_target_m2)
{
    if (!j || p_r_W <= 0.0 || rcs_target_m2 <= 0.0) return 0.0;
    /* Solve J/S = snr_threshold for R:
     *   R² = (snr_threshold · P_t · G_t · σ) / (P_j · G_j · 4π)
     *   R  = √(numer / denom). */
    const double G_j = db_to_linear_(j->g_j_db);
    const double G_t = db_to_linear_(g_r_t_db);
    const double numer = j->snr_threshold * p_r_W * G_t * rcs_target_m2;
    const double denom = j->p_j_W * G_j * 4.0 * K26A_PI;
    if (denom <= 0.0) return 0.0;
    return sqrt(numer / denom);
}

double k26astro_jammer_self_signature_W(const K26AstroJammer *j)
{
    return j ? j->p_j_W : 0.0;
}
