/* test_impactor_physics.c — closed-form gates for Whipple ballistic-
 * limit, Modified Cour-Palais monolithic, swarm dispatch geometry,
 * and terminal-phase math. */

#include "k26astro_impactor/impactor.h"
#include "k26astro_core/consts.h"
#include "k26compute.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx_(double got, double want, double rel_tol)
{
    const double denom = (fabs(want) > 1.0) ? fabs(want) : 1.0;
    return fabs(got - want) / denom <= rel_tol;
}

int main(void)
{
    /* ---- Whipple ballistic-limit (Christiansen 1993 Eq. 11) --- *
     *
     * Anchor case: t_w = 1 mm, ρ_b = 2700 (Al bumper), ρ_p = 2700,
     * V = 7 km/s, cos = 1, S = 0.1 m, σ_w = 70 ksi (Al 6061-T6).
     *
     * Native cgs/ksi form (page 150):
     *   d_c[cm] = 3.918 × (0.1)^{2/3} × 2.7^{-1/9} × 2.7^{-1/3}
     *             × 7^{-2/3} × 10^{1/3} × 1
     *           = 3.918 × 0.2154 × 0.8956 × 0.7181 × 0.2733 × 2.154
     *           ≈ 0.320 cm = 3.20e-3 m
     *
     * SI form using K_CHRISTIANSEN_SI ≈ 8441 gives the same value. */
    const double tw = 1e-3;
    const double rho_b = 2700.0;
    const double rho_p = 2700.0;
    const double V = 7000.0;
    const double S = 0.1;
    const double sigma_ksi = 70.0;
    const double d_c = k26astro_whipple_critical_diameter(
        tw, rho_b, rho_p, V, 1.0, S, sigma_ksi);
    const double expected =
        K26ASTRO_IMPACTOR_K_CHRISTIANSEN_SI
        * pow(tw, 2.0 / 3.0)
        * pow(rho_b, -1.0 / 9.0)
        * pow(rho_p, -1.0 / 3.0)
        * pow(V, -2.0 / 3.0)
        * pow(S, 1.0 / 3.0)
        * pow(sigma_ksi / 70.0, 1.0 / 3.0);
    assert(approx_(d_c, expected, 1.0e-12));
    /* Anchor: Christiansen 1993 Eq. 11 closed-form at canonical inputs. */
    assert(approx_(d_c, 3.20e-3, 0.05));

    /* Projectile diameter larger than d_c → perforates. */
    assert(k26astro_whipple_penetrates(
        tw, rho_b, /* d_proj */ 5e-3, rho_p, V, 1.0, S, sigma_ksi) == 1);
    /* Projectile diameter smaller than d_c → no perforation. */
    assert(k26astro_whipple_penetrates(
        tw, rho_b, /* d_proj */ 1e-6, rho_p, V, 1.0, S, sigma_ksi) == 0);

    /* Oblique impact (cos θ < 1) raises d_c → harder to perforate. */
    const double d_c_oblique = k26astro_whipple_critical_diameter(
        tw, rho_b, rho_p, V, /* cos */ 0.5, S, sigma_ksi);
    assert(d_c_oblique > d_c);

    /* ---- Monolithic penetration (Hayashida-Robinson §2.5) ----- *
     *
     * Modified Cour-Palais form (NASA TM-103565 page 4):
     *   p[cm] = 5.24 · d[cm]^(19/18) · BH^(-0.25) · (ρ_p/ρ_t)^0.5
     *           · (V_n / C)^(2/3)
     *
     * Anchor case: m = 0.01 kg steel sphere → d = (6m/(π·ρ_p))^(1/3)
     * ≈ 1.345e-2 m = 1.345 cm; BH = 95 (Al 6061-T6); ρ_p = 7850;
     * ρ_t = 2700; V = 7 km/s; cos = 1; C_Al = 5100 m/s.
     *
     *   ρ_p/ρ_t = 2.907 → (ρ_p/ρ_t)^0.5 = 1.705
     *   V_n/C   = 7/5.1 = 1.373 → (V_n/C)^(2/3) = 1.235
     *   d^(19/18) = 1.345^1.0556 = 1.368
     *   BH^(-0.25) = 95^(-0.25) = 0.320
     *   p = 5.24 · 1.368 · 0.320 · 1.705 · 1.235 ≈ 4.83 cm
     *     ≈ 0.0483 m */
    const double d_proj = 1.345e-2;
    const double BH = 95.0;
    const double rho_t = 2700.0;
    const double C_target = 5100.0;
    const double P = k26astro_monolithic_penetration_depth(
        d_proj, BH, 7850.0, rho_t, V, 1.0, C_target);
    /* Recoded form for ulp comparison. */
    const double V_n_km = V * 1.0e-3;
    const double C_km   = C_target * 1.0e-3;
    const double P_expected =
        5.24
        * pow(d_proj * 1.0e2, 19.0 / 18.0)
        * pow(BH, -0.25)
        * pow(7850.0 / rho_t, 0.5)
        * pow(V_n_km / C_km, 2.0 / 3.0)
        * 1.0e-2;
    assert(approx_(P, P_expected, 1.0e-12));
    /* Anchor: Hayashida-Robinson §2.5 closed-form at canonical inputs. */
    assert(approx_(P, 0.0483, 0.05));

    /* Doubling V multiplies P by 2^{2/3} ≈ 1.587 (the V_n/C
     * exponent in the canonical form is 2/3). */
    const double P_2V = k26astro_monolithic_penetration_depth(
        d_proj, BH, 7850.0, rho_t, 14000.0, 1.0, C_target);
    assert(approx_(P_2V / P, pow(2.0, 2.0 / 3.0), 1.0e-12));

    /* ---- Swarm dispatch geometry --------------------------- */
    /* Per-unit mass: 50 / 100 = 0.5 kg. */
    assert(k26astro_swarm_per_unit_mass(100, 50.0) == 0.5);
    /* Swarm footprint at 10 km, 10 mrad half-angle:
     *   r = 10000 × sin(0.01) ≈ 99.998 m
     *   A = π r^2 ≈ 31415 m^2 */
    const double A_fp = k26astro_swarm_footprint_m2(1.0e4, 0.01);
    const double r_fp = 1.0e4 * sin(0.01);
    const double A_fp_expected = K26A_PI * r_fp * r_fp;
    assert(approx_(A_fp, A_fp_expected, 1.0e-12));

    /* Sampling: an aligned cone (axis = +z, narrow half-angle)
     * produces directions with z-component near 1. Sample over
     * 1000 draws and check that all are within the cone. */
    K26CRng rng;
    k26c_rng_init(&rng, 42);
    K26V3 axis = { 0.0, 0.0, 1.0 };
    const double half = 0.05;
    const double cos_half = cos(half);
    for (int i = 0; i < 1000; i++) {
        K26V3 d = k26astro_swarm_sample_direction(axis, half, &rng);
        /* Unit vector. */
        const double len = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
        assert(fabs(len - 1.0) < 1.0e-9);
        /* Inside the cone (dot product with axis >= cos(half_angle)). */
        assert(d.z >= cos_half - 1.0e-12);
    }

    /* NULL rng returns the cone axis. */
    K26V3 axis_only = k26astro_swarm_sample_direction(axis, half, NULL);
    assert(axis_only.x == 0.0 && axis_only.y == 0.0 && axis_only.z == 1.0);

    /* Bit-identical reproduction with same seed. */
    K26CRng rng_a, rng_b;
    k26c_rng_init(&rng_a, 12345);
    k26c_rng_init(&rng_b, 12345);
    for (int i = 0; i < 100; i++) {
        K26V3 da = k26astro_swarm_sample_direction(axis, 0.02, &rng_a);
        K26V3 db = k26astro_swarm_sample_direction(axis, 0.02, &rng_b);
        assert(da.x == db.x && da.y == db.y && da.z == db.z);
    }

    /* ---- Terminal-phase geometry ---------------------------- *
     *
     * Two bodies on a head-on collision: projectile at (-1e6, 0, 0)
     * heading +x at 10 km/s; target at (1e6, 0, 0) heading -x at
     * 5 km/s. Closing speed = 15 km/s. */
    K26V3 r_p = { -1.0e6, 0.0, 0.0 };
    K26V3 v_p = {  1.0e4, 0.0, 0.0 };
    K26V3 r_t = {  1.0e6, 0.0, 0.0 };
    K26V3 v_t = { -5.0e3, 0.0, 0.0 };
    const double v_close = k26astro_impactor_closing_speed(v_p, v_t);
    assert(approx_(v_close, 1.5e4, 1.0e-12));

    /* Time to closest approach: |r_p - r_t| = 2e6 m, closing at
     * 1.5e4 m/s → t ≈ 133.33 s. */
    const double tca = k26astro_impactor_time_to_closest_approach(
        r_p, v_p, r_t, v_t);
    assert(approx_(tca, 2.0e6 / 1.5e4, 1.0e-12));

    /* Closest-approach distance for the head-on collision: 0. */
    const double cad = k26astro_impactor_closest_approach_distance(
        r_p, v_p, r_t, v_t);
    assert(cad < 1.0e-6);

    /* Off-axis offset case: projectile passes target with miss
     * distance 1000 m. */
    K26V3 r_p2 = { -1.0e6, 1000.0, 0.0 };
    K26V3 v_p2 = {  1.0e4, 0.0,    0.0 };
    K26V3 r_t2 = {  1.0e6, 0.0,    0.0 };
    K26V3 v_t2 = { -5.0e3, 0.0,    0.0 };
    const double miss = k26astro_impactor_closest_approach_distance(
        r_p2, v_p2, r_t2, v_t2);
    assert(approx_(miss, 1000.0, 1.0e-9));

    /* Impact cos-angle: -closing-velocity = +x; surface normal
     * = +x; cos = 1 (head-on). */
    K26V3 closing_unit = { 1.0, 0.0, 0.0 };
    K26V3 surface_norm = { -1.0, 0.0, 0.0 };  /* facing projectile */
    const double cos_imp = k26astro_impactor_impact_cos_angle(
        closing_unit, surface_norm);
    assert(approx_(cos_imp, 1.0, 1.0e-12));

    /* ---- Integrated impactor_analyse_impact ---------------------- */
    K26AstroImpactor *k = k26astro_impactor_new(
        K26ASTRO_IMPACTOR_PATTERN_SINGLE,
        /* m */ 0.01,
        /* ρ_p */ K26ASTRO_IMPACTOR_RHO_TUNGSTEN_KG_PER_M3,
        /* d */  0.01,
        0, 0.0);
    K26AstroImpactEvent ev = k26astro_impactor_analyse_impact(
        k,
        /* V */ 1.0e4,
        /* cos */ 1.0,
        /* t_w */ 0.001,
        /* ρ_b */ K26ASTRO_IMPACTOR_BUMPER_DENSITY_AL_KG_PER_M3,
        /* S */   0.1,
        /* σ_ksi */ K26ASTRO_IMPACTOR_WALL_YIELD_AL6061T6_KSI,
        /* BH */  95.0,
        /* ρ_target */ K26ASTRO_IMPACTOR_RHO_ALUMINUM_KG_PER_M3,
        /* C_target */ 5100.0);
    assert(ev.critical_diameter_m > 0.0);
    assert(ev.monolithic_penetration_m > 0.0);
    /* A 10 mm tungsten projectile at 10 km/s easily perforates a
     * 1 mm aluminium Whipple shield. */
    assert(ev.penetrates == 1);
    k26astro_impactor_destroy(k);

    printf("test_impactor_physics: OK\n");
    return 0;
}
