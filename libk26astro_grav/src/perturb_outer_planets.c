/* perturb_outer_planets.c — Ephemeris-driven outer-planet perturbation.
 *
 * Adds GM·(r_outer - r_inner) / |r_outer - r_inner|³ acceleration
 * from the outer planets (Jupiter, Saturn; optionally Uranus,
 * Neptune) onto each integrated body. The outer planets themselves
 * are not in the integrated `view->bodies` list; their positions
 * come from the DE441 ephemeris at the integrator's current sim
 * epoch (state->t). They are "on rails" in the sense that their
 * trajectories are external truth, not numerically integrated by
 * the K26 force pipeline.
 *
 * Use case: inner-solar-system 100-year integrations. Without
 * outer-planet perturbations, Jupiter's secular forcing alone
 * accumulates ~1e10 to 1e11 m drift on the inner planets over 100
 * years. With Jupiter + Saturn on rails, the residual drops by an
 * order of magnitude toward the ~100-1000 km regime expected from
 * a realistic inner-planet integration.
 *
 * Construction: caller allocates K26AstroOuterPlanetsCtx, fills in
 * the ephemeris handle + NAIF ids + GMs, then registers via
 * k26astro_grav_register_perturb(state, k26astro_perturb_outer_planets, ctx).
 * The ctx is owned by the caller (state doesn't free it on
 * destroy). */

#include "k26astro_grav/perturb.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/perturb_outer_planets.h"
#include "k26astro_core/pos.h"
#include "k26astro_ephem/ephem.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>

void k26astro_perturb_outer_planets(const K26AstroGravState *state,
                                    const K26AstroGravView  *view,
                                    K26V3 *accel_out, void *ctx)
{
    if (!state || !view || !accel_out || !ctx) return;
    const K26AstroOuterPlanetsCtx *opc =
        (const K26AstroOuterPlanetsCtx *)ctx;
    if (!opc->ephem || opc->n_outer <= 0) return;

    /* Cache outer-planet positions at this integrator step. The
     * ephem query is per-step (called once per gravity evaluation),
     * not per-body; the same outer-planet snapshot applies to all
     * inner bodies' acceleration. */
    K26V3   r_outer[K26ASTRO_MAX_OUTER_PLANETS];
    double  gm_outer[K26ASTRO_MAX_OUTER_PLANETS];
    int     n_active = 0;
    for (int op = 0; op < opc->n_outer; ++op) {
        K26AstroStateXV s = k26astro_ephem_body_state(opc->ephem,
                                                       opc->naif_ids[op],
                                                       &state->t);
        K26V3 r = k26astro_pos_to_m_approx(&s.pos);
        /* Sentinel: ephem returns zero state on a query miss (kernel
         * doesn't cover the time, or NAIF id is unknown). Skip the
         * planet rather than applying nonsense acceleration. */
        if (r.x == 0.0 && r.y == 0.0 && r.z == 0.0 &&
            s.vel.x == 0.0 && s.vel.y == 0.0 && s.vel.z == 0.0) {
            continue;
        }
        r_outer[n_active]  = r;
        gm_outer[n_active] = opc->gms[op];
        n_active++;
    }

    /* Heliocentric tidal-form perturbation:
     *   a_eff_on_i = GM_outer · [(r_outer - r_i)/|r_outer - r_i|³
     *                            - r_outer/|r_outer|³]
     *
     * The first term is the direct pull of the outer planet on body i;
     * the second is the indirect term: the outer planet pulls the
     * central body (Sun), which shifts the heliocentric frame and
     * appears as a fictitious force on every other body. For bodies
     * near the Sun (~1 AU) being perturbed by Jupiter (~5 AU), the
     * direct and indirect terms partially cancel, leaving only the
     * differential ("tidal") pull. Without the indirect subtraction,
     * the perturbation overcounts by ~10x for inner planets.
     *
     * Indirect-term r_outer is measured from the Sun (body 0) since
     * the WH heliocentric frame is Sun-centred. */
    K26V3 r_sun = { 0.0, 0.0, 0.0 };
    if (view->n > 0) r_sun = k26astro_pos_to_m_approx(&view->bodies[0].pos);

    /* Precompute the indirect term per outer planet:
     *   a_indirect_op = GM_outer · (r_outer - r_sun) / |r_outer - r_sun|³ */
    K26V3 a_indirect[K26ASTRO_MAX_OUTER_PLANETS];
    for (int op = 0; op < n_active; ++op) {
        const double dx = r_outer[op].x - r_sun.x;
        const double dy = r_outer[op].y - r_sun.y;
        const double dz = r_outer[op].z - r_sun.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r  = sqrt(r2);
        const double pref = (r2 > 0.0) ? (gm_outer[op] / (r2 * r)) : 0.0;
        a_indirect[op].x = pref * dx;
        a_indirect[op].y = pref * dy;
        a_indirect[op].z = pref * dz;
    }

    /* Apply tidal-form perturbation to every NON-CENTRAL body. The
     * central body (idx 0) stays at the heliocentric origin in WH;
     * no perturbation applied to it; the heliocentric formulation
     * already accounts for the central body's motion via the
     * indirect term subtraction. */
    for (int i = 1; i < view->n; ++i) {
        K26V3 r_inner = k26astro_pos_to_m_approx(&view->bodies[i].pos);
        K26V3 a_total = { 0.0, 0.0, 0.0 };
        for (int op = 0; op < n_active; ++op) {
            const double dx = r_outer[op].x - r_inner.x;
            const double dy = r_outer[op].y - r_inner.y;
            const double dz = r_outer[op].z - r_inner.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 <= 0.0) continue;
            const double r = sqrt(r2);
            const double pref = gm_outer[op] / (r2 * r);
            /* Direct - indirect. */
            a_total.x += pref * dx - a_indirect[op].x;
            a_total.y += pref * dy - a_indirect[op].y;
            a_total.z += pref * dz - a_indirect[op].z;
        }
        accel_out[i].x += a_total.x;
        accel_out[i].y += a_total.y;
        accel_out[i].z += a_total.z;
    }
}

int k26astro_grav_enable_outer_planets(K26AstroGravState *state,
                                       K26AstroOuterPlanetsCtx *ctx)
{
    if (!state || !ctx) return -1;
    if (!ctx->ephem || ctx->n_outer <= 0 ||
        ctx->n_outer > K26ASTRO_MAX_OUTER_PLANETS) {
        return -1;
    }
    return k26astro_grav_register_perturb(state,
                                          k26astro_perturb_outer_planets,
                                          ctx);
}

/* IAU 2015 / DE441 nominal GM values for the four outer-planet
 * barycentres (m³/s²). Folkner et al. 2014 §3 published values,
 * hex-pinned for cross-platform determinism.
 *
 * Jupiter  barycentre GM = 1.26686534e17 m³/s²
 *   IEEE-754 bits: 0x434C19DE3A800000
 * Saturn   barycentre GM = 3.7931187e16 m³/s²
 *   IEEE-754 bits: 0x4340D7D4FC400000
 * Uranus   barycentre GM = 5.793939e15  m³/s²
 *   IEEE-754 bits: 0x4334958E8BF1FE00
 * Neptune  barycentre GM = 6.836529e15  m³/s²
 *   IEEE-754 bits: 0x433849C9728AAA00 */
static double k_hex_d_(uint64_t bits)
{
    union { double d; uint64_t u; } cvt;
    cvt.u = bits;
    return cvt.d;
}

/* Singleton ctx for the convenience default-registrar. One per
 * process; multi-world programs that need distinct outer-planet
 * configs (different ephem files, etc.) must use the explicit
 * k26astro_grav_enable_outer_planets() with their own
 * K26AstroOuterPlanetsCtx allocations.
 *
 * Lifetime: process-scoped. OS reclaims on exit; no destroy hook
 * needed because the singleton outlives any grav_state. */
static K26AstroOuterPlanetsCtx k_default_outer_ctx_;

int k26astro_grav_enable_outer_planets_default(K26AstroGravState *state,
                                                K26AstroEphem    *ephem)
{
    if (!state || !ephem) return -1;
    k_default_outer_ctx_.ephem    = ephem;
    k_default_outer_ctx_.n_outer  = 4;
    k_default_outer_ctx_.naif_ids[0] = 5;  /* Jupiter barycentre */
    k_default_outer_ctx_.naif_ids[1] = 6;  /* Saturn barycentre  */
    k_default_outer_ctx_.naif_ids[2] = 7;  /* Uranus barycentre  */
    k_default_outer_ctx_.naif_ids[3] = 8;  /* Neptune barycentre */
    k_default_outer_ctx_.gms[0] = k_hex_d_(0x434C19DE3A800000ULL);
    k_default_outer_ctx_.gms[1] = k_hex_d_(0x4340D7D4FC400000ULL);
    k_default_outer_ctx_.gms[2] = k_hex_d_(0x4334958E8BF1FE00ULL);
    k_default_outer_ctx_.gms[3] = k_hex_d_(0x433849C9728AAA00ULL);
    return k26astro_grav_enable_outer_planets(state, &k_default_outer_ctx_);
}
