/* phipps_coupling.h — laser-target coupling (Phipps).
 *
 * Two regimes:
 *
 *   1. Below plasma-ignition fluence (Φ < Φ_th): vapour-pressure
 *      coupling. C_m is roughly an order of magnitude smaller
 *      than the plasma regime; ablation is slow.
 *
 *   2. Above threshold (Φ >= Φ_th): plasma-mediated coupling.
 *      The ablation plume forms an inverse-bremsstrahlung
 *      absorption layer; C_m scales as ~(I λ √τ)^(1/4) per
 *      Phipps 2010 eq. (3) but the lib uses tabulated nominal
 *      values rather than the scaling relation by default —
 *      reflecting the measurement vs model discipline of the
 *      published data.
 *
 * The lib exposes both regimes through a single API: the
 * `c_m_coupling` function returns the appropriate C_m for the
 * caller-supplied (material, wavelength, fluence). Mass loss is
 * derived from the coupled energy and the material's specific
 * ablation energy Q*. */
#ifndef K26ASTRO_LASER_PHIPPS_COUPLING_H
#define K26ASTRO_LASER_PHIPPS_COUPLING_H

#include "k26astro_laser/laser_consts.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Momentum-coupling coefficient C_m in N/W for the given
 * material, beam wavelength, and incident fluence. Caller may
 * pass material = K26ASTRO_LASER_MAT_NONE to use a conservative
 * generic-metal default. */
double k26astro_phipps_c_m(K26AstroLaserMaterial material,
                           double wavelength_nm,
                           double fluence_J_per_m2);

/* Plasma-ignition fluence threshold (J/m^2). Below this fluence,
 * vapour-pressure coupling applies; above, plasma coupling. */
double k26astro_phipps_threshold_fluence(K26AstroLaserMaterial material,
                                          double wavelength_nm);

/* Specific ablation energy Q* (J/kg). Mass loss for energy E
 * coupled into the target is E / Q*. */
double k26astro_phipps_q_star(K26AstroLaserMaterial material);

/* Convenience: predicted mass loss (kg) from dwell time τ at
 * uniform intensity I (W/m^2) onto a spot of area A (m^2). The
 * coupled energy is E = (1 - reflectivity) · I · A · τ; mass
 * loss is E / Q* once Φ = I · τ exceeds the plasma threshold,
 * otherwise the vapour-pressure regime applies with a 10x
 * smaller mass yield. */
double k26astro_phipps_dwell_mass_loss(K26AstroLaserMaterial material,
                                        double wavelength_nm,
                                        double intensity_W_per_m2,
                                        double spot_area_m2,
                                        double dwell_s,
                                        double reflectivity);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_LASER_PHIPPS_COUPLING_H */
