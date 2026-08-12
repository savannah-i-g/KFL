/* impactor_consts.h — Christiansen + Hayashida-Robinson defaults.
 *
 * Hypervelocity-impact scaling for Whipple-shielded and monolithic
 * targets is well-measured in the 1-50 km/s regime. The lib uses
 * the Christiansen 1993 IJIE Eq. 11 ballistic-limit equation
 * (NNO regime, V_n ≥ 7 km/sec) for Whipple shields. Monolithic
 * penetration is the Modified Cour-Palais form reported in
 * Hayashida & Robinson 1991 §2.5.
 *
 * References:
 *   Christiansen 1993          IJIE 14:145-156, Eq. 11
 *   Hayashida & Robinson 1991  NASA TM-103565, §2.5 Modified Cour-Palais
 *
 * Default densities are nominal SI values for common shield /
 * projectile materials. Mission fits override per-vehicle.
 */
#ifndef K26ASTRO_IMPACTOR_CONSTS_H
#define K26ASTRO_IMPACTOR_CONSTS_H

/* Common material densities (kg / m^3). */
#define K26ASTRO_IMPACTOR_RHO_ALUMINUM_KG_PER_M3      2700.0
#define K26ASTRO_IMPACTOR_RHO_STEEL_KG_PER_M3         7850.0
#define K26ASTRO_IMPACTOR_RHO_TITANIUM_KG_PER_M3      4506.0
#define K26ASTRO_IMPACTOR_RHO_TUNGSTEN_KG_PER_M3      19300.0
#define K26ASTRO_IMPACTOR_RHO_KEVLAR_KG_PER_M3        1440.0
#define K26ASTRO_IMPACTOR_RHO_NEXTEL_KG_PER_M3        2700.0

/* Christiansen Eq. 11 constant in SI units. Native form (cgs/ksi):
 *   d_c[cm] = 3.918 · t_w[cm]^(2/3) · ρ_b[g/cm³]^(-1/9) ·
 *             ρ_p[g/cm³]^(-1/3) · V[km/s]^(-2/3) · S[cm]^(1/3) ·
 *             (σ_w[ksi] / 70)^(1/3)
 *
 * Converting to SI (m, kg/m³, m/s) folds 10^(10/3) ≈ 2154.4 into K:
 *   K_SI = 3.918 · 10^(10/3) ≈ 8441
 *
 * The (σ_w / 70) factor stays dimensionless because σ_w is taken
 * in ksi throughout — there is no SI σ_w form in Eq. 11. */
#define K26ASTRO_IMPACTOR_K_CHRISTIANSEN_SI           8441.0

/* Bumper density defaults (kg/m³) — outer-layer material in a
 * Whipple shield. Mission fits override per-target. */
#define K26ASTRO_IMPACTOR_BUMPER_DENSITY_AL_KG_PER_M3   2700.0
#define K26ASTRO_IMPACTOR_BUMPER_DENSITY_TI_KG_PER_M3   4506.0

/* Rear-wall yield-stress defaults (ksi). Christiansen's
 * normalisation pivot is 70 ksi (Al 6061-T6); other materials
 * scale via (σ/70)^(1/3). */
#define K26ASTRO_IMPACTOR_WALL_YIELD_AL6061T6_KSI       70.0
#define K26ASTRO_IMPACTOR_WALL_YIELD_AL2024T3_KSI       65.0
#define K26ASTRO_IMPACTOR_WALL_YIELD_TI6AL4V_KSI        130.0
#define K26ASTRO_IMPACTOR_WALL_YIELD_STEEL_MILD_KSI     50.0

/* Spacing between front bumper and rear wall (m), nominal.
 * Christiansen 2009 ballistic-limit equation assumes the bumper
 * fragments the projectile and the expanded debris cloud
 * spreads over this distance. Default 0.1 m. */
#define K26ASTRO_IMPACTOR_DEFAULT_BUMPER_SPACING_M    0.10

/* Default swarm half-angle (rad). The dispersion cone
 * is set by the projectile-ejection geometry of the carrier
 * vehicle; mission fits override. */
#define K26ASTRO_IMPACTOR_DEFAULT_SWARM_HALF_ANGLE_RAD  0.01

/* Default number of swarm units. Mission fits override. */
#define K26ASTRO_IMPACTOR_DEFAULT_SWARM_COUNT     100

#endif /* K26ASTRO_IMPACTOR_CONSTS_H */
