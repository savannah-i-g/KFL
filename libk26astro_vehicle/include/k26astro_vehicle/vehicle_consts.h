/* k26astro_vehicle/vehicle_consts.h — vehicle-tier mass / inertia
 * canonical constants and constructor defaults.
 *
 * Constructor defaults are consumed by k26astro_vehicle_new. The
 * reference-fixture section is the home for canonical spacecraft
 * mass / inertia values (e.g. Falcon 9, Apollo Command Module,
 * Cassini) used by example programs and by reference-fixture tests.
 *
 * Discipline:
 *
 *   - Hex-literal IEEE-754 form alongside the decimal form for any
 *     value where bit-exact reproducibility across CPU vendors is
 *     load-bearing (cross-vendor reference comparison). Decimal
 *     forms are sufficient for magnitudes that round-trip through
 *     strtod on both glibc and musl within 1 ULP — see the
 *     ias15_coeffs.c rationale in libk26astro_grav.
 *   - Each fixture cites its source measurement (manufacturer
 *     technical report, NASA mission press kit, post-flight
 *     reconstruction paper). No undocumented numbers.
 *
 * References:
 *
 *   - NASA-STD-1000, "Mass Properties Control for Space Systems."
 *   - GSFC-STD-1000RevH, "Rules for the Design, Development,
 *     Verification, and Operation of Flight Systems."
 *   - AIAA S-120A-2015, "Mass Properties Control for Space Systems." */
#ifndef K26ASTRO_VEHICLE_CONSTS_H
#define K26ASTRO_VEHICLE_CONSTS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Constructor defaults ------------------------------------- */

/* Initial NASA basic mass on a fresh vehicle. Zero so that a caller
 * who never sets a dry mass produces a clearly-broken vehicle
 * (mass queries return 0) rather than a silently-wrong default. */
#define K26ASTRO_VEHICLE_DEFAULT_BASIC_MASS_KG  0.0

/* Initial AIAA S-120A-2015 mass growth allowance. Zero so that
 * programs not tracking flight-class mass budgets see
 * predicted_mass == basic_mass. */
#define K26ASTRO_VEHICLE_DEFAULT_MGA_KG         0.0

/* ---- Reference fixtures --------------------------------------- *
 *
 * Canonical spacecraft constants belong here. Format:
 *
 *   #define K26ASTRO_VEHICLE_<ship>_BASIC_MASS_KG  549054.0  // 0x4110B5C1A0000000
 *   #define K26ASTRO_VEHICLE_<ship>_IXX_KG_M2      ...
 *
 * Each entry cites its measurement source in an adjacent comment. */

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_VEHICLE_CONSTS_H */
