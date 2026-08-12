/* laser.c — laser opaque + engage() pipeline. */

#include "k26astro_laser/laser.h"

#include "k26astro_core/consts.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct K26AstroLaser {
    K26AstroPayload base;
    double primary_diam_m;
    double wavelength_m;
    double p_output_W;
    double m_squared;
    double pointing_jitter_rad;
    double rms_wavefront_m;
    double plasma_attn_k;
};

K26AstroLaser *k26astro_laser_new(
    double primary_diam_m,
    double wavelength_nm,
    double p_output_W,
    double m_squared,
    double pointing_jitter_rad,
    double rms_wavefront_m,
    double plasma_attn_k)
{
    if (primary_diam_m <= 0.0 || wavelength_nm <= 0.0) return NULL;
    if (p_output_W <= 0.0) return NULL;
    if (m_squared < 1.0) m_squared = 1.0;
    if (pointing_jitter_rad < 0.0) pointing_jitter_rad = 0.0;
    if (rms_wavefront_m < 0.0) rms_wavefront_m = 0.0;
    if (plasma_attn_k < 0.0) plasma_attn_k = K26ASTRO_LASER_PLASMA_PLUG_DEFAULT_K;

    K26AstroLaser *l = (K26AstroLaser *)calloc(1, sizeof(*l));
    if (!l) return NULL;
    l->base.kind             = K26ASTRO_DEFENSE_KIND_LASER;
    l->base.owner            = NULL;
    l->base.owner_slot       = -1;
    l->base.generation       = 1;
    l->base.on_owner_destroy = NULL;
    l->primary_diam_m        = primary_diam_m;
    l->wavelength_m          = wavelength_nm * 1.0e-9;
    l->p_output_W            = p_output_W;
    l->m_squared             = m_squared;
    l->pointing_jitter_rad   = pointing_jitter_rad;
    l->rms_wavefront_m       = rms_wavefront_m;
    l->plasma_attn_k         = plasma_attn_k;
    return l;
}

void k26astro_laser_destroy(K26AstroLaser *l)
{
    if (!l) return;
    k26astro_payload_unlink_(&l->base);
    free(l);
}

int k26astro_laser_bind(K26AstroLaser *l, struct K26AstroVehicle *v)
{
    if (!l) return -1;
    return k26astro_payload_bind(&l->base, v);
}

struct K26AstroPayload *k26astro_laser_as_payload(K26AstroLaser *l)
{
    return l ? &l->base : NULL;
}

K26AstroLaserAblationEvent k26astro_laser_engage(
    const K26AstroLaser   *emitter,
    K26AstroLaserMaterial  target_material,
    double                 target_area_m2,
    double                 target_reflectivity,
    double                 range_m,
    double                 dwell_s)
{
    K26AstroLaserAblationEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.range_m = range_m;

    if (!emitter || range_m <= 0.0 || dwell_s <= 0.0) return ev;
    if (target_area_m2 <= 0.0) return ev;
    if (target_reflectivity < 0.0) target_reflectivity = 0.0;
    if (target_reflectivity > 1.0) target_reflectivity = 1.0;

    /* 1. Spot geometry. */
    ev.spot_diameter_m = k26astro_airy_spot_diameter(
        range_m, emitter->wavelength_m, emitter->primary_diam_m,
        emitter->m_squared, emitter->pointing_jitter_rad);
    /* Target-radius bound on spot (whichever is smaller bounds the
     * encircled-energy calculation). */
    const double r_target = sqrt(target_area_m2 / K26A_PI);
    ev.encircled_fraction = k26astro_airy_encircled_fraction(
        r_target, range_m, emitter->primary_diam_m, emitter->wavelength_m);

    /* 2. Strehl scaling. */
    const double strehl = k26astro_mirror_strehl(
        emitter->rms_wavefront_m, emitter->wavelength_m);
    ev.p_at_aperture_W = emitter->p_output_W * strehl;

    /* 3. On-target intensity = encircled_fraction × P / target_area. */
    const double p_on_target = ev.p_at_aperture_W * ev.encircled_fraction;
    ev.on_target_intensity_W_per_m2 = p_on_target / target_area_m2;
    ev.fluence_J_per_m2 = ev.on_target_intensity_W_per_m2 * dwell_s;

    /* 4. Plasma-plug attenuation. */
    const double Phi_th = k26astro_phipps_threshold_fluence(
        target_material, emitter->wavelength_m * 1.0e9);
    ev.threshold_fluence_J_per_m2 = Phi_th;
    ev.plasma_ignited = (Phi_th > 0.0) &&
                        (ev.fluence_J_per_m2 >= Phi_th);
    ev.plasma_transmissivity = k26astro_plasma_plug_transmissivity(
        ev.fluence_J_per_m2, Phi_th, emitter->plasma_attn_k);

    /* 5. Coupled energy. */
    const double absorbed = (1.0 - target_reflectivity) *
                            ev.plasma_transmissivity;
    ev.p_coupled_W = absorbed * p_on_target;
    const double E_coupled = ev.p_coupled_W * dwell_s;

    /* 6. Mass loss + impulse. */
    const double Q_star = k26astro_phipps_q_star(target_material);
    ev.mass_loss_kg = (Q_star > 0.0) ? (E_coupled / Q_star) : 0.0;
    /* Below the plasma threshold, mass loss derates by 10×
     * (vapour-pressure regime), matching k26astro_phipps_dwell_mass_loss. */
    if (Phi_th > 0.0 && ev.fluence_J_per_m2 < Phi_th) {
        ev.mass_loss_kg *= 0.1;
    }

    const double C_m = k26astro_phipps_c_m(target_material,
        emitter->wavelength_m * 1.0e9, ev.fluence_J_per_m2);
    ev.impulse_N_s = C_m * E_coupled;
    ev.effective_dwell_s = dwell_s;

    return ev;
}
