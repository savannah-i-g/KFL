/* phipps_coupling.c — laser-target momentum + energy coupling.
 *
 * Default values are nominal Phipps 1988 + 2010 table entries.
 * Wavelength scaling on C_m follows Phipps 2010 Eq. 19
 * (page 616): C_m ∝ (Iλ√τ)^(-1/4), so the wavelength dependence
 * alone is C_m(λ) ∝ λ^(-1/4) at fixed (I, τ). Mission fit tables
 * override per-engagement.
 *
 * Plasma-ignition fluence threshold Φ_th is tabulated at 1064 nm.
 * Phipps 2010 §III.B (page 616) gives I·τ^(1/2) = constant for
 * the plasma-ignition intensity; no closed-form wavelength
 * dependence is supplied. Per-wavelength thresholds must come
 * from caller-supplied tables. */

#include "k26astro_laser/phipps_coupling.h"
#include "k26astro_laser/laser_consts.h"
#include "k26astro_laser/plasma_plug.h"

#include <math.h>

static double material_cm_default_(K26AstroLaserMaterial m)
{
    switch (m) {
    case K26ASTRO_LASER_MAT_ALUMINUM:
        return K26ASTRO_LASER_CM_DEFAULT_AL_N_PER_W;
    case K26ASTRO_LASER_MAT_STEEL:
        return K26ASTRO_LASER_CM_DEFAULT_STEEL_N_PER_W;
    case K26ASTRO_LASER_MAT_TITANIUM:
        return K26ASTRO_LASER_CM_DEFAULT_TITANIUM_N_PER_W;
    case K26ASTRO_LASER_MAT_COPPER:
        return K26ASTRO_LASER_CM_DEFAULT_COPPER_N_PER_W;
    case K26ASTRO_LASER_MAT_COMPOSITE:
        return K26ASTRO_LASER_CM_DEFAULT_COMPOSITE_N_PER_W;
    case K26ASTRO_LASER_MAT_FUSED_SILICA:
        return K26ASTRO_LASER_CM_DEFAULT_FUSED_SILICA_N_PER_W;
    default:
        /* Conservative generic-metal anchor. */
        return K26ASTRO_LASER_CM_DEFAULT_STEEL_N_PER_W;
    }
}

static double wavelength_scale_(double wavelength_nm)
{
    /* Phipps 2010 Eq. 19 (page 616): C_m = 1.84e-4 · Ψ^(9/16) ·
     * A^(1/8) / (Iλ√τ)^(1/4). Holding (I, τ) fixed and varying λ
     * alone, the scaling is C_m(λ) / C_m(1064) = (1064 / λ)^(1/4).
     * Reference wavelength is the 1064 nm Nd:YAG table value. */
    if (wavelength_nm <= 0.0) return 1.0;
    return pow(1064.0 / wavelength_nm, 0.25);
}

double k26astro_phipps_c_m(K26AstroLaserMaterial material,
                           double wavelength_nm,
                           double fluence_J_per_m2)
{
    const double C_m_plasma = material_cm_default_(material) *
                              wavelength_scale_(wavelength_nm);
    const double Phi_th = k26astro_phipps_threshold_fluence(
        material, wavelength_nm);
    if (Phi_th <= 0.0) return C_m_plasma;
    if (fluence_J_per_m2 < Phi_th) {
        /* Sub-threshold vapour-pressure regime: ~0.1× plasma C_m. */
        return 0.1 * C_m_plasma;
    }
    return C_m_plasma;
}

double k26astro_phipps_threshold_fluence(K26AstroLaserMaterial material,
                                          double wavelength_nm)
{
    /* Tabulated thresholds are at 1064 nm. Phipps 2010 §III.B (page
     * 616) gives I·τ^(1/2) = constant for plasma-ignition intensity
     * (pulse-width scaling); no wavelength dependence appears in the
     * cited paragraph. The lib therefore returns the 1064 nm value
     * for all wavelengths and leaves wavelength-specific thresholds
     * to caller-supplied tables. */
    if (wavelength_nm <= 0.0) return 0.0;
    switch (material) {
    case K26ASTRO_LASER_MAT_ALUMINUM:
        return K26ASTRO_LASER_THRESHOLD_AL_J_PER_M2_AT_1064NM;
    case K26ASTRO_LASER_MAT_STEEL:
        return K26ASTRO_LASER_THRESHOLD_STEEL_J_PER_M2_AT_1064NM;
    case K26ASTRO_LASER_MAT_TITANIUM:
        return K26ASTRO_LASER_THRESHOLD_TITANIUM_J_PER_M2_AT_1064NM;
    case K26ASTRO_LASER_MAT_COPPER:
        return K26ASTRO_LASER_THRESHOLD_COPPER_J_PER_M2_AT_1064NM;
    case K26ASTRO_LASER_MAT_COMPOSITE:
        return K26ASTRO_LASER_THRESHOLD_COMPOSITE_J_PER_M2_AT_1064NM;
    case K26ASTRO_LASER_MAT_FUSED_SILICA:
        return K26ASTRO_LASER_THRESHOLD_FUSED_SILICA_J_PER_M2_AT_1064NM;
    default:
        return K26ASTRO_LASER_THRESHOLD_STEEL_J_PER_M2_AT_1064NM;
    }
}

double k26astro_phipps_q_star(K26AstroLaserMaterial material)
{
    switch (material) {
    case K26ASTRO_LASER_MAT_ALUMINUM:
        return K26ASTRO_LASER_Q_STAR_AL_J_PER_KG;
    case K26ASTRO_LASER_MAT_STEEL:
        return K26ASTRO_LASER_Q_STAR_STEEL_J_PER_KG;
    case K26ASTRO_LASER_MAT_TITANIUM:
        return K26ASTRO_LASER_Q_STAR_TITANIUM_J_PER_KG;
    case K26ASTRO_LASER_MAT_COPPER:
        return K26ASTRO_LASER_Q_STAR_COPPER_J_PER_KG;
    case K26ASTRO_LASER_MAT_COMPOSITE:
        return K26ASTRO_LASER_Q_STAR_COMPOSITE_J_PER_KG;
    case K26ASTRO_LASER_MAT_FUSED_SILICA:
        return K26ASTRO_LASER_Q_STAR_FUSED_SILICA_J_PER_KG;
    default:
        return K26ASTRO_LASER_Q_STAR_STEEL_J_PER_KG;
    }
}

double k26astro_phipps_dwell_mass_loss(K26AstroLaserMaterial material,
                                        double wavelength_nm,
                                        double intensity_W_per_m2,
                                        double spot_area_m2,
                                        double dwell_s,
                                        double reflectivity)
{
    if (intensity_W_per_m2 <= 0.0 || spot_area_m2 <= 0.0) return 0.0;
    if (dwell_s <= 0.0) return 0.0;
    if (reflectivity < 0.0) reflectivity = 0.0;
    if (reflectivity > 1.0) reflectivity = 1.0;

    const double absorbed_frac = 1.0 - reflectivity;
    const double power_in_W = intensity_W_per_m2 * spot_area_m2;
    const double E_coupled = absorbed_frac * power_in_W * dwell_s;
    const double Q_star = k26astro_phipps_q_star(material);
    if (Q_star <= 0.0) return 0.0;

    /* Check whether the cumulative fluence exceeds the plasma
     * threshold; if not, derate by 10× for the vapour-pressure
     * regime. */
    const double fluence = intensity_W_per_m2 * dwell_s;
    const double Phi_th = k26astro_phipps_threshold_fluence(
        material, wavelength_nm);
    double mass = E_coupled / Q_star;
    if (Phi_th > 0.0 && fluence < Phi_th) mass *= 0.1;
    return mass;
}
