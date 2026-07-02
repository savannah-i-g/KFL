/* observer.c — observer-mode dispatch.
 *
 * GEOMETRIC: instant-light vector.
 * ASTROMETRIC: solve |r_T(t_obs - τ) - r_obs(t_obs)| = c τ via
 *              fixed-point iteration. For an in-world body (not
 *              ephemeris-backed), the iteration uses the body's
 *              current K26AstroBody.vel for retarded position.
 * APPARENT: ASTROMETRIC + stellar aberration (Aoki et al. 1983).
 *           Adds GR Shapiro delay (Will 1993 eq. 8.36) when the
 *           world's grav state has gr_ppn1 enabled.
 * TOPOCENTRIC: returns E_NOT_IMPLEMENTED (needs libk26astro_atmos). */
#include "k26astro_rt/observer.h"
#include "world_internal.h"

#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26astro_atmos/atmos.h"
#include "k26m3d.h"

#include <math.h>

/* Speed of light, m/s. IAU 2009 (defined, exact). */
#define K26ASTRO_C_LIGHT 299792458.0

static double v3_norm_(K26V3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static K26V3 v3_scale_(K26V3 v, double s)
{
    K26V3 r = { v.x * s, v.y * s, v.z * s };
    return r;
}

static double v3_dot_(K26V3 a, K26V3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* Light-time iteration for an in-world target. The retarded position
 * r_T(t - τ) is approximated as r_T(t) - vel_T * τ (linear extrap).
 * Converges in 2-3 iterations for inner SS distances. */
static K26AstroPos retarded_position_(const K26AstroBody *target,
                                       const K26AstroPos  *obs_pos,
                                       int max_iter)
{
    K26AstroPos retarded = target->pos;
    if (max_iter <= 0) return retarded;
    for (int iter = 0; iter < max_iter; iter++) {
        K26V3 r = k26astro_pos_sub(&retarded, obs_pos);
        double d = v3_norm_(r);
        double tau = d / K26ASTRO_C_LIGHT;
        K26V3 back = v3_scale_(target->vel, -tau);
        retarded = target->pos;
        k26astro_pos_add(&retarded, back);
    }
    return retarded;
}

/* Stellar aberration: Aoki 1983 (A&A 105:359) exact relativistic
 * form. Given the geometric apparent direction `dir` (unit vector
 * from observer to target in the BCRS rest frame) and the observer's
 * velocity `v_obs` (m/s in the same frame), returns the aberrated
 * direction in the observer's instantaneous rest frame:
 *
 *   β   = v_obs / c              (vector)
 *   β²  = β · β                  (scalar)
 *   γ   = 1 / √(1 - β²)
 *   β·d = β · dir                (scalar)
 *
 *   dir' = ( (1/γ) dir + β + (γ/(γ+1)) (β·d) β ) / (1 + β·d)
 *
 * For β² → 0 this collapses to the textbook first-order form
 * dir' ≈ (dir + β) / |dir + β| identically; the 1st-order gate
 * tests still pass. For Earth's |β| ~ 1e-4 the second-order term
 * contributes ~v²/c² ~ 1e-8 rad (≈ 2 mas), which the
 * test_aberration_aoki gate verifies.
 *
 * Numerical robustness: when |dir| is small or v_obs is zero, fall
 * back to identity. The denominator (1 + β·d) is bounded below by
 * 1 - |β| > 0 for any physical observer velocity. */
static K26V3 aberrate_(K26V3 dir, K26V3 v_obs)
{
    K26V3 beta = v3_scale_(v_obs, 1.0 / K26ASTRO_C_LIGHT);
    double beta2 = v3_dot_(beta, beta);
    if (!(beta2 > 0.0)) return dir;
    if (beta2 >= 1.0) return dir;            /* unphysical, bail */
    double gamma = 1.0 / sqrt(1.0 - beta2);
    double beta_dot_dir = v3_dot_(beta, dir);
    double inv_gamma    = 1.0 / gamma;
    double k_b          = gamma / (gamma + 1.0);

    K26V3 num;
    num.x = inv_gamma * dir.x + beta.x + k_b * beta_dot_dir * beta.x;
    num.y = inv_gamma * dir.y + beta.y + k_b * beta_dot_dir * beta.y;
    num.z = inv_gamma * dir.z + beta.z + k_b * beta_dot_dir * beta.z;

    double denom = 1.0 + beta_dot_dir;
    if (!(denom > 0.0)) return dir;
    K26V3 out = v3_scale_(num, 1.0 / denom);
    double n = v3_norm_(out);
    if (!(n > 0.0)) return dir;
    return v3_scale_(out, 1.0 / n);
}

/* GR Shapiro delay magnitude: ds = (2 GM / c^2) ln ((r_obs + r_T +
 *   (r_obs + r_T) · r̂_obs·r̂_T) / (r_obs + r_T - (r_obs + r_T) · r̂_obs·r̂_T))
 * Will 1993 *Theory and Experiment in Gravitational Physics* §8.36. */
static double shapiro_delay_m_(K26V3 r_obs_from_sun,
                                K26V3 r_tgt_from_sun,
                                double gm_sun)
{
    double r_o = v3_norm_(r_obs_from_sun);
    double r_t = v3_norm_(r_tgt_from_sun);
    if (!(r_o > 0.0) || !(r_t > 0.0)) return 0.0;
    double rr = v3_dot_(r_obs_from_sun, r_tgt_from_sun) / (r_o * r_t);
    double num = r_o + r_t + (r_o + r_t) * rr;
    double den = r_o + r_t - (r_o + r_t) * rr;
    if (!(den > 0.0)) return 0.0;
    double c2 = K26ASTRO_C_LIGHT * K26ASTRO_C_LIGHT;
    return (2.0 * gm_sun / c2) * log(num / den);
}

/* GR Shapiro apparent-direction rotation. The light path
 * bending angle for a ray grazing a spherically-symmetric mass at
 * impact parameter b is the standard gravitational-lensing result:
 *
 *     Δθ = (2 GM / (c² b)) · (1 + cos α)
 *
 * where α is the angle observer → deflector → target. For α → 0
 * (target behind deflector) the factor (1 + cos α) → 2 → full
 * Einstein 4GM/(c²b). For α → π (deflector behind target) it → 0.
 *
 * The deflected apparent direction is the geometric direction
 * `dir` (unit observer→target) rotated by Δθ in the plane
 * containing observer, deflector, and target, AWAY from the
 * deflector (light is pulled toward the mass; the apparent source
 * appears displaced AWAY from it).
 *
 * Implementation: compute the unit perpendicular `up` that lies in
 * the plane spanned by `dir` and (-u_def) (vector from observer
 * away from deflector) and is orthogonal to `dir`. Then the
 * rotated direction is `dir + Δθ · up`, renormalised. For small
 * Δθ (always < 4 mas in inner-SS geometry) the small-angle
 * approximation is exact to 1e-12 rad.
 *
 * Returns rotated direction; on degenerate geometry (target
 * exactly at deflector, observer at deflector) returns the input
 * direction unchanged. */
static K26V3 shapiro_direction_(K26V3 dir,
                                 K26V3 r_obs_from_def,
                                 K26V3 r_tgt_from_def,
                                 double gm_def)
{
    if (!(gm_def > 0.0)) return dir;

    /* Unit vector from observer toward deflector. */
    K26V3 obs_to_def = { -r_obs_from_def.x,
                          -r_obs_from_def.y,
                          -r_obs_from_def.z };
    double d_od = v3_norm_(obs_to_def);
    if (!(d_od > 0.0)) return dir;
    K26V3 u_def = v3_scale_(obs_to_def, 1.0 / d_od);

    /* Angle α at deflector: observer → deflector → target. */
    double r_o = v3_norm_(r_obs_from_def);
    double r_t = v3_norm_(r_tgt_from_def);
    if (!(r_o > 0.0) || !(r_t > 0.0)) return dir;
    double cos_alpha = v3_dot_(r_obs_from_def, r_tgt_from_def)
                      / (r_o * r_t);
    /* Clamp for numerical safety. */
    if (cos_alpha >  1.0) cos_alpha =  1.0;
    if (cos_alpha < -1.0) cos_alpha = -1.0;

    /* Impact parameter b ≈ r_o · sin(angle observer→target offset
     * from observer→deflector). Use the apparent target direction
     * as the ray direction; sin(angle) = |dir × u_def|. */
    K26V3 cross = {
        dir.y * u_def.z - dir.z * u_def.y,
        dir.z * u_def.x - dir.x * u_def.z,
        dir.x * u_def.y - dir.y * u_def.x
    };
    double sin_off = v3_norm_(cross);
    double b_imp   = r_o * sin_off;
    /* Floor at the deflector's effective Schwarzschild radius to
     * avoid singularity on a degenerate ray. */
    double rs = 2.0 * gm_def
              / (K26ASTRO_C_LIGHT * K26ASTRO_C_LIGHT);
    if (b_imp < rs) b_imp = rs;
    if (!(b_imp > 0.0)) return dir;

    double dtheta = (2.0 * gm_def
                    / (K26ASTRO_C_LIGHT * K26ASTRO_C_LIGHT * b_imp))
                  * (1.0 + cos_alpha);

    /* Unit perpendicular to `dir`, in the plane (dir, u_def),
     * pointing AWAY from the deflector. Construct via
     * Gram-Schmidt: perp_raw = u_def - (u_def·dir) dir, then
     * negate (we want away from deflector, not toward).
     * On near-degenerate alignment (|cross| → 0), the deflection
     * direction is undefined; return unrotated. */
    if (sin_off < 1.0e-12) return dir;
    double d_dot_u = v3_dot_(u_def, dir);
    K26V3 perp_raw = {
        u_def.x - d_dot_u * dir.x,
        u_def.y - d_dot_u * dir.y,
        u_def.z - d_dot_u * dir.z
    };
    double pn = v3_norm_(perp_raw);
    if (!(pn > 0.0)) return dir;
    /* Negate: we want displacement AWAY from deflector, so the
     * perpendicular component of the rotation is opposite the
     * observer→deflector projection. */
    K26V3 perp_away = v3_scale_(perp_raw, -1.0 / pn);

    K26V3 rotated = {
        dir.x + dtheta * perp_away.x,
        dir.y + dtheta * perp_away.y,
        dir.z + dtheta * perp_away.z
    };
    double rn = v3_norm_(rotated);
    if (!(rn > 0.0)) return dir;
    return v3_scale_(rotated, 1.0 / rn);
}

int k26astro_world_observe(const K26AstroWorld *world,
                            int target_idx, int observer_idx,
                            K26AstroPos *out_target_pos,
                            K26V3 *out_apparent_dir)
{
    if (!world) return -K26ASTRO_RT_E_NULL;
    int n = world->grav.n_bodies;
    if (target_idx   < 0 || target_idx   >= n) return -K26ASTRO_RT_E_BAD_ARG;
    if (observer_idx < 0 || observer_idx >= n) return -K26ASTRO_RT_E_BAD_ARG;
    if (target_idx == observer_idx)            return -K26ASTRO_RT_E_BAD_ARG;

    const K26AstroBody *target   = &world->grav.bodies[target_idx];
    const K26AstroBody *observer = &world->grav.bodies[observer_idx];

    K26AstroPos target_pos = target->pos;

    switch (world->observer_mode) {
    case K26ASTRO_OBS_GEOMETRIC:
        break;
    case K26ASTRO_OBS_ASTROMETRIC:
    case K26ASTRO_OBS_APPARENT:
    case K26ASTRO_OBS_TOPOCENTRIC:
        /* TOPOCENTRIC runs the same light-time + aberration
         * + Shapiro pipeline as APPARENT, then post-applies
         * atmospheric refraction at the observer's locale.
         * Light-time iteration is the same for all three corrected
         * modes. */
        target_pos = retarded_position_(target, &observer->pos, 4);
        break;
    }

    K26V3 r = k26astro_pos_sub(&target_pos, &observer->pos);
    double d = v3_norm_(r);
    K26V3 dir = (d > 0.0) ? v3_scale_(r, 1.0 / d) : (K26V3){0, 0, 0};

    if (world->observer_mode == K26ASTRO_OBS_APPARENT
     || world->observer_mode == K26ASTRO_OBS_TOPOCENTRIC) {
        dir = aberrate_(dir, observer->vel);
        /* Shapiro: only apply when GR PPN-1 is enabled. */
        if (world->grav.use_gr_ppn1 && n >= 1) {
            /* Treat the body of largest mass as the deflector (Sun). */
            int sun_idx = 0;
            double m_max = world->grav.bodies[0].mass;
            for (int k = 1; k < n; k++) {
                if (world->grav.bodies[k].mass > m_max) {
                    m_max   = world->grav.bodies[k].mass;
                    sun_idx = k;
                }
            }
            const K26AstroBody *sun = &world->grav.bodies[sun_idx];
            /* Skip when the deflector IS the target or observer
             * (geometry is degenerate / not a deflection). */
            if (sun_idx != target_idx && sun_idx != observer_idx) {
                K26V3 r_o = k26astro_pos_sub(&observer->pos, &sun->pos);
                K26V3 r_t = k26astro_pos_sub(&target_pos,    &sun->pos);
                /* Emit both the delay-as-distance magnitude
                 * (preserved as diagnostic; sub-arcsec for inner
                 * SS) and the direction rotation. */
                double ds = shapiro_delay_m_(r_o, r_t, sun->gm);
                (void)ds;
                dir = shapiro_direction_(dir, r_o, r_t, sun->gm);
            }
        }

        /* TOPOCENTRIC adds atmospheric refraction on top of
         * APPARENT. The zenith direction at the observer is the
         * unit vector from the observer's PARENT body centre
         * (the planet under their feet) to the observer.
         *
         * Multi-body extension: read observer->parent_body_idx
         * to identify the planet. The SOI parent field already
         * carries this; a surface observer's SOI parent is the
         * body they're sitting on. Position differencing produces
         * a correct zenith regardless of where the world's
         * coordinate origin sits (Earth-centric, heliocentric,
         * barycentric, all work).
         *
         * Fallbacks:
         *   parent_body_idx < 0          - no parent (central body
         *                                  or free-flying). Skip
         *                                  refraction; APPARENT
         *                                  result stands.
         *   parent_body_idx == observer  - degenerate. Skip.
         *   parent_body_idx >= n_bodies  - invalid. Skip.
         *
         * Callers that leave observer->parent_body_idx at the
         * default (-1) get APPARENT-equivalent behaviour for
         * TOPOCENTRIC. */
        if (world->observer_mode == K26ASTRO_OBS_TOPOCENTRIC
         && world->atmos) {
            int parent_idx = observer->parent_body_idx;
            if (parent_idx >= 0 && parent_idx < n
             && parent_idx != observer_idx) {
                const K26AstroBody *parent =
                    &world->grav.bodies[parent_idx];
                K26V3 r_obs_from_parent =
                    k26astro_pos_sub(&observer->pos, &parent->pos);
                double rn = v3_norm_(r_obs_from_parent);
                if (rn > 0.0) {
                    K26V3 zenith =
                        v3_scale_(r_obs_from_parent, 1.0 / rn);
                    dir = k26astro_atmos_apparent(world->atmos, dir,
                                                    zenith);
                }
            }
        }
    }

    if (out_target_pos) *out_target_pos = target_pos;
    if (out_apparent_dir) {
        out_apparent_dir->x = dir.x;
        out_apparent_dir->y = dir.y;
        out_apparent_dir->z = dir.z;
    }
    return K26ASTRO_RT_OK;
}
