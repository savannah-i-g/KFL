/* k26astro_laser/laser.h — directed-energy laser.
 *
 * Composes against libk26astro_vehicle's PAYLOAD composition slot
 * via the K26AstroPayload shared base. The concrete struct is
 * private to src/laser.c; consumers see only the opaque
 * K26AstroLaser typedef plus the new / destroy / engage entry
 * points.
 *
 * The engage() function combines the four submodules — Airy
 * diffraction geometry, Phipps coupling, plasma-plug attenuation,
 * mirror Strehl — into a single ablation-event predictor that
 * the mission driver invokes per tick during an engage decision.
 *
 * Determinism: laser ablation is fully deterministic (no
 * stochastic noise on the predicted impulse / mass loss); the
 * laser does not consume an RNG.
 *
 * Stated limitations:
 *   - Plasma-plug attenuation factor is a caller-supplied
 *     first-order parameter; the empirical anchor is weak.
 *   - Mirror wavefront state is caller-supplied; thermal
 *     distortion under sustained flux is not modeled internally.
 *   - Pointing jitter is treated as Gaussian-RMS spreading the
 *     spot; biased mount offset is not modeled.
 *   - The encircled-energy fraction is the ideal diffraction-
 *     limited Airy curve; M^2 and jitter widen the reported spot
 *     diameter but do not degrade the encircled fraction, so
 *     energy coupling is optimistic for low-quality beams.
 *   - Atmospheric path attenuation is out of scope (vacuum
 *     engagement-bubble regime); atmospheric callers should
 *     pair this lib with libk26astro_atmos.
 *
 * References:
 *   Born & Wolf 2002    Principles of Optics, 7th ed.
 *   Phipps 1989         Impulse coupling in vacuum (KrF/HF/CO2)
 *   Phipps 2010         Laser-ablation propulsion review
 */
#ifndef K26ASTRO_LASER_LASER_H
#define K26ASTRO_LASER_LASER_H

#include "k26astro_laser/laser_consts.h"
#include "k26astro_laser/airy.h"
#include "k26astro_laser/phipps_coupling.h"
#include "k26astro_laser/plasma_plug.h"
#include "k26astro_laser/mirror_strehl.h"

#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque laser. */
typedef struct K26AstroLaser K26AstroLaser;

/* Construct a laser.
 *
 *   primary_diam_m      primary aperture diameter
 *   wavelength_nm       operating wavelength (nm)
 *   p_output_W          continuous output power (W)
 *   m_squared           beam quality factor (>= 1)
 *   pointing_jitter_rad RMS pointing jitter
 *   rms_wavefront_m     mirror wavefront RMS error
 *   plasma_attn_k       plasma-plug attenuation strength
 *
 * Returns NULL on invalid inputs or allocation failure. */
K26AstroLaser *k26astro_laser_new(
    double primary_diam_m,
    double wavelength_nm,
    double p_output_W,
    double m_squared,
    double pointing_jitter_rad,
    double rms_wavefront_m,
    double plasma_attn_k);

/* Destroy. Unlinks from the bound vehicle (if any) and frees the
 * concrete struct. Safe on NULL. */
void k26astro_laser_destroy(K26AstroLaser *l);

/* Bind to a vehicle's PAYLOAD composition slot. Returns 0 on
 * success, -1 on NULL inputs or slot append failure. */
int k26astro_laser_bind(K26AstroLaser *l, struct K26AstroVehicle *v);

/* Downcast to the shared K26AstroPayload base. The base is the
 * laser's first member; the cast is a pointer reinterpret. */
struct K26AstroPayload *k26astro_laser_as_payload(K26AstroLaser *l);

/* ---- Ablation event ------------------------------------------ */

typedef struct {
    /* Geometry */
    double range_m;
    double spot_diameter_m;
    double encircled_fraction;

    /* Power flow */
    double p_at_aperture_W;        /* after Strehl + jitter */
    double on_target_intensity_W_per_m2;
    double fluence_J_per_m2;       /* intensity × dwell */
    double threshold_fluence_J_per_m2; /* Phi_th(material) at 1064 nm */
    double plasma_transmissivity;  /* in [0, 1] */
    double p_coupled_W;

    /* Output */
    double mass_loss_kg;
    double impulse_N_s;            /* C_m · E_coupled */
    double effective_dwell_s;

    /* Regime flag. True iff fluence at target >= Phi_th(material) at
     * 1064 nm; false signals vapour-regime, mass_loss already derated
     * 10× internally (Phipps 2010 §III.B vapour-pressure coupling). */
    int    plasma_ignited;
} K26AstroLaserAblationEvent;

/* Engage a target. Computes the full ablation chain in one call:
 *   1. Spot geometry at range R (Airy + M^2 + jitter)
 *   2. Strehl scaling of the on-aperture power
 *   3. On-target power = encircled_fraction · p_at_aperture;
 *      intensity = on-target power / target_area
 *   4. Plasma-plug transmissivity above threshold
 *   5. Coupled energy E_coupled = (1 - reflectivity) · T_plasma ·
 *      on-target power · τ
 *   6. Mass loss E_coupled / Q*; impulse C_m · E_coupled
 *
 * The encircled_fraction reports the fraction of total beam
 * energy that falls inside the target's projected area; the
 * intensity calculation uses this fraction (rather than assuming
 * uniform spot illumination). */
K26AstroLaserAblationEvent k26astro_laser_engage(
    const K26AstroLaser   *emitter,
    K26AstroLaserMaterial  target_material,
    double                 target_area_m2,
    double                 target_reflectivity,
    double                 range_m,
    double                 dwell_s);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_LASER_LASER_H */
