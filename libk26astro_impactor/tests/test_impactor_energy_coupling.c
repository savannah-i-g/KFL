/* test_impactor_energy_coupling.c — impactor energy-coupling API gate.
 *
 * Validates k26astro_impactor_energy_delivered_j across three regimes:
 *
 *   1. Bumper-stopped Whipple — only surface_spallation_eff
 *      (~10%) of impact KE delivered.
 *   2. Whipple bumper breached, inner wall present — bumper_partial_eff
 *      (~50%) delivered.
 *   3. Full penetration of monolithic target — full_penetration_eff
 *      (~90%) delivered.
 *
 * Also asserts:
 *   - Energy bounded by 0 ≤ E ≤ impact_KE.
 *   - Partial-monolithic-penetration scales quadratically with
 *     depth ratio.
 *   - NULL inputs return 0.
 *
 * Rationale: without this API a driver would assign delivered
 * energy as raw ½mv² regardless of target shielding; with it the
 * delivered energy scales with the penetration outcome. */
#include "k26astro_impactor/energy_coupling.h"
#include "k26astro_impactor/impactor.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx_rel(double got, double want, double rel_tol)
{
    const double denom = (fabs(want) > 1.0) ? fabs(want) : 1.0;
    return fabs(got - want) / denom <= rel_tol;
}

int main(void)
{
    /* 50 kg tungsten impactor (ρ=19,250 kg/m³) at 7 km/s impact. */
    const double m_proj = 50.0;
    const double rho_w  = 19250.0;
    const double d_proj = 0.108;  /* ~10 cm diameter sphere */
    K26AstroImpactor *impactor = k26astro_impactor_new(
        K26ASTRO_IMPACTOR_PATTERN_SINGLE, m_proj, rho_w, d_proj, 1, 0.0);
    assert(impactor != NULL);

    const double impact_KE = 0.5 * m_proj * 7000.0 * 7000.0;  /* 1.225 GJ */
    fprintf(stderr,
            "test_impactor_energy_coupling: impact_KE = %.3e J\n",
            impact_KE);

    /* ---- Regime 1: Bumper-stopped Whipple --------------------- *
     *
     * Thick bumper (3 mm Al), small projectile (10 mm) → Whipple
     * critical diameter exceeds projectile diameter → bumper
     * stops it. Only surface spallation delivers energy inward. */
    K26AstroImpactEvent ev_stopped = {
        .impact_velocity_m_per_s   = 7000.0,
        .cos_impact_angle          = 1.0,
        .critical_diameter_m       = 0.150,   /* > d_proj */
        .penetrates                = 0,
        .monolithic_penetration_m  = 0.0,
    };
    K26AstroTargetStructureSpec target_whipple_thick = {
        .outer_thickness_m         = 3.0e-3,
        .outer_yield_stress_ksi    = 35.0,
        .bumper_density_kg_per_m3  = 2780.0,
        .bumper_gap_m              = 0.10,
        .inner_thickness_m         = 5.0e-3,
        .inner_yield_stress_ksi    = 50.0,
        /* monolithic params zero */
        /* coupling effs zero → defaults applied */
    };
    double E_stopped = k26astro_impactor_energy_delivered_j(
        impactor, &ev_stopped, &target_whipple_thick);
    /* Default surface_spallation_eff = 0.10. */
    assert(approx_rel(E_stopped, 0.10 * impact_KE, 1.0e-12));
    fprintf(stderr,
            "test_impactor_energy_coupling: bumper-stopped Whipple — "
            "E_delivered = %.3e J (%.1f%% of impact_KE)\n",
            E_stopped, 100.0 * E_stopped / impact_KE);

    /* ---- Regime 2: Whipple bumper breached, inner wall present  *
     *
     * Larger projectile, thinner bumper → Whipple critical
     * diameter LESS than projectile diameter → bumper breached.
     * Inner wall present → bumper_partial_eff (default 0.50). */
    K26AstroImpactEvent ev_breach = {
        .impact_velocity_m_per_s   = 7000.0,
        .cos_impact_angle          = 1.0,
        .critical_diameter_m       = 0.020,   /* < d_proj */
        .penetrates                = 1,
        .monolithic_penetration_m  = 0.0,
    };
    double E_breach = k26astro_impactor_energy_delivered_j(
        impactor, &ev_breach, &target_whipple_thick);
    assert(approx_rel(E_breach, 0.50 * impact_KE, 1.0e-12));
    fprintf(stderr,
            "test_impactor_energy_coupling: Whipple breach + inner wall — "
            "E_delivered = %.3e J (%.1f%% of impact_KE)\n",
            E_breach, 100.0 * E_breach / impact_KE);

    /* ---- Regime 3: Full monolithic penetration --------------- *
     *
     * Monolithic plate, projectile depth exceeds wall thickness. */
    K26AstroImpactEvent ev_full_mono = {
        .impact_velocity_m_per_s   = 7000.0,
        .cos_impact_angle          = 1.0,
        .critical_diameter_m       = 0.0,
        .penetrates                = 0,
        .monolithic_penetration_m  = 0.50,
    };
    K26AstroTargetStructureSpec target_mono_thin = {
        /* Whipple disabled (zero outer_thickness). */
        .monolithic_thickness_m    = 0.05,    /* < penetration */
        .monolithic_brinell_hardness    = 150.0,
        .monolithic_density_kg_per_m3   = 7800.0,
        .monolithic_speed_of_sound_m_per_s = 5000.0,
    };
    double E_full_mono = k26astro_impactor_energy_delivered_j(
        impactor, &ev_full_mono, &target_mono_thin);
    assert(approx_rel(E_full_mono, 0.90 * impact_KE, 1.0e-12));
    fprintf(stderr,
            "test_impactor_energy_coupling: full monolithic penetration — "
            "E_delivered = %.3e J (%.1f%% of impact_KE)\n",
            E_full_mono, 100.0 * E_full_mono / impact_KE);

    /* ---- Regime 4: Partial monolithic penetration ------------ *
     *
     * Quadratic falloff with depth ratio. Penetration = 0.5 ×
     * wall thickness → delivered = 0.25 × 0.90 × impact_KE. */
    K26AstroImpactEvent ev_half_mono = {
        .impact_velocity_m_per_s   = 7000.0,
        .cos_impact_angle          = 1.0,
        .critical_diameter_m       = 0.0,
        .penetrates                = 0,
        .monolithic_penetration_m  = 0.25,   /* = half wall */
    };
    K26AstroTargetStructureSpec target_mono_thick = {
        .monolithic_thickness_m    = 0.50,
        .monolithic_brinell_hardness    = 150.0,
        .monolithic_density_kg_per_m3   = 7800.0,
        .monolithic_speed_of_sound_m_per_s = 5000.0,
    };
    double E_half = k26astro_impactor_energy_delivered_j(
        impactor, &ev_half_mono, &target_mono_thick);
    assert(approx_rel(E_half, 0.25 * 0.90 * impact_KE, 1.0e-12));

    /* ---- Bounds: 0 ≤ E ≤ impact_KE ---------------------------- */
    assert(E_stopped >= 0.0 && E_stopped <= impact_KE);
    assert(E_breach >= 0.0 && E_breach <= impact_KE);
    assert(E_full_mono >= 0.0 && E_full_mono <= impact_KE);
    assert(E_half >= 0.0 && E_half <= impact_KE);

    /* ---- NULL safety ----------------------------------------- */
    assert(k26astro_impactor_energy_delivered_j(NULL, &ev_breach,
                                            &target_whipple_thick) == 0.0);
    assert(k26astro_impactor_energy_delivered_j(impactor, NULL,
                                            &target_whipple_thick) == 0.0);
    assert(k26astro_impactor_energy_delivered_j(impactor, &ev_breach,
                                            NULL) == 0.0);

    /* ---- Custom coupling efficiencies ------------------------ */
    K26AstroTargetStructureSpec target_custom_eff = {
        .outer_thickness_m         = 3.0e-3,
        .bumper_gap_m              = 0.10,
        .inner_thickness_m         = 5.0e-3,
        .surface_spallation_eff    = 0.05,  /* tougher structure */
        .bumper_partial_eff        = 0.30,
        .full_penetration_eff      = 0.95,
    };
    double E_custom_stopped = k26astro_impactor_energy_delivered_j(
        impactor, &ev_stopped, &target_custom_eff);
    assert(approx_rel(E_custom_stopped, 0.05 * impact_KE, 1.0e-12));

    k26astro_impactor_destroy(impactor);
    fprintf(stderr, "test_impactor_energy_coupling: OK\n");
    return 0;
}
