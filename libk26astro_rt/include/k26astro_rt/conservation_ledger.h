/* libk26astro_rt — closed-system conservation ledger.
 *
 * The world's body / vehicle subsystem does not track the energy +
 * momentum carried out of the surviving-vehicle population by
 * mass-ejection and radiative emission events and propulsion
 * exhaust. Without those terms the trace-side conservation
 * diagnostic (sum over surviving vehicles) is a vanity metric:
 * arbitrary mass + KE may leave the system as solid-mass ejecta,
 * radiated energy, or reaction-mass exhaust and the surviving-bodies
 * bookkeeping still reads as bit-perfect.
 *
 * The ledger is a pure accumulator: callers (typically a scenario
 * driver) report each emission at the moment it occurs; the
 * ledger maintains running totals (E_ejecta, p_ejecta_xyz,
 * L_z_ejecta, E_radiated, E_reaction_mass, p_reaction_xyz,
 * E_chem_source) that the trace writer combines with the
 * surviving-bodies subtotals to get the closed-system invariant
 * for a conservation gate.
 *
 * The ledger holds no references to the world or to any vehicle
 * handle. It is a passive bookkeeping object — register-only,
 * never reads back into the simulation state. Determinism is
 * trivial: deterministic-add into a fixed-shape struct.
 *
 * Notation: kinetic-energy KE = ½ m v² uses |v|²; momentum is
 * three-vector m·v; angular momentum L_z = m (r_x v_y − r_y v_x).
 *
 * ── FRAME CONVENTION (load-bearing) ──
 *
 * Every velocity vector passed into the record_* APIs must be in
 * the same inertial frame as the surviving-vehicles bookkeeping
 * that the driver sums against the ledger. The closed-system
 * invariant
 *
 *     E_closed = sum_i ½ m_i v_i²        (surviving vehicles)
 *              + E_ejecta                (ledger)
 *              + E_radiated              (ledger)
 *              + E_reaction_mass         (ledger)
 *              - E_chem_source           (ledger, see below)
 *
 * is a meaningful invariant only when every term is in the same
 * frame. Passing rocket-frame exhaust velocity into a ledger that
 * also receives inertial-frame surviving-bodies KE leaves a
 * per-burn drift on the order of m_p · v_vehicle · v_e — for
 * Orion-class burns at 30 km/s vehicle velocity and 30 km/s
 * exhaust, that is ~10¹⁴ J per burn per vessel, easily 10¹⁵ J
 * across a six-vessel squadron's full schedule.
 *
 * ── E_reaction_mass_J interpretation trap ──
 *
 * "E_reaction_mass_J" is unambiguous as a *ledger field* — it is
 * the inertial-frame kinetic energy of ejected propellant, i.e.
 *
 *     E_reaction_mass_J += ½ m_p |v_p_inertial|²
 *
 * where v_p_inertial is the propellant's velocity in the inertial
 * frame after ejection.
 *
 * Three different physical quantities sometimes get the same label
 * in writeups and casual prose; do not conflate them:
 *
 *   (a) Exhaust-frame propellant KE   = ½ m_p v_e²
 *       What the legacy _record_reaction_mass computes IF the caller
 *       passes the rocket-frame exhaust velocity. NOT what the
 *       closed-system invariant needs.
 *   (b) Inertial-frame propellant KE  = ½ m_p |V_vehicle − v_e d̂|²
 *       What this ledger's E_reaction_mass_J means.
 *       This is the canonical interpretation.
 *   (c) Vehicle ΔKE                   = ½ M_v (V_after² − V_before²)
 *       A characterisation of the burn's effect on the vehicle.
 *       Useful for writeup prose but NOT what the ledger records.
 *
 * For a Tsiolkovsky burn at Δv/v_e ≈ 0.3-0.5 these three values
 * differ by factors of 3-10. Conflating them silently in a
 * conservation ledger produces drift on the order of the
 * difference.
 *
 * Prefer the high-level k26astro_rt_ledger_record_burn helper for
 * new code — it takes the burn's physical inputs (vehicle inertial
 * velocity, exhaust speed, exhaust direction) and computes both
 * E_reaction_mass_J (interpretation b) and the chemistry source
 * term atomically, with no opportunity for frame confusion at the
 * call site.
 *
 * ── E_chem_source ──
 *
 * When a propulsion burn converts chemical (or nuclear) potential
 * energy into kinetic energy, the surviving-vehicle KE + the
 * ejected propellant KE together exceed their pre-burn sum by the
 * chemistry energy released. Callers record that chemistry input
 * (typically ½ m_p v_e², the propellant's rocket-frame kinetic
 * energy at the nozzle exit) via record_chemical_source. The
 * closed-system invariant subtracts E_chem_source from the
 * post-burn sum to recover the pre-burn baseline.
 *
 * E_chem_source has no momentum field by design. Chemistry release
 * in the rocket frame produces equal-and-opposite momentum on the
 * vehicle and exhaust; the rocket equation captures that
 * redistribution in the surviving-vehicles + p_reaction
 * bookkeeping. There is no net momentum input from chemistry in
 * any frame, so no separate p_chem_source_* field is needed.
 *
 * ── Documented limitations ──
 *
 * - Radiated energy bookkeeping is one-sided. record_radiated
 *   credits the receiver with absorbed photons but does not debit
 *   the emitter vehicle's KE budget. For a 1 MW × 5 s emitter at
 *   55 emissions the missing emitter-side term is ~275 MJ,
 *   negligible against engagement-class drift; for sustained
 *   high-power radiative emission this becomes a real bookkeeping
 *   gap. Caller can balance manually via record_dispatch with a
 *   photon-energy entry if needed.
 * - L_z_ejecta tracks angular momentum about the z-axis only
 *   (matching the surviving-bodies conservation.csv L_z column).
 *   3-D angular momentum bookkeeping is a future extension. */
#ifndef K26ASTRO_RT_CONSERVATION_LEDGER_H
#define K26ASTRO_RT_CONSERVATION_LEDGER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle. Allocated by ledger_new, freed by ledger_destroy. */
typedef struct K26AstroRtConservationLedger K26AstroRtConservationLedger;

/* Plain-old-data totals snapshot. All values are running sums since
 * ledger_new. */
typedef struct {
    /* Ejecta (solid-mass emissions): mass × velocity contributions
     * from every recorded emission since ledger_new. */
    double E_ejecta_J;
    double p_ejecta_x_kgmps;
    double p_ejecta_y_kgmps;
    double p_ejecta_z_kgmps;

    /* Angular momentum about the z-axis carried away by ejecta.
     * Computed as mass × (origin_x v_y − origin_y v_x) at the moment
     * of emission. Closes the L_z column in the trace-side
     * conservation.csv when ejected masses are removed from the
     * propagator. */
    double L_z_ejecta_kg_m2ps;

    /* Radiated energy: cumulative joules emitted radiatively. No
     * momentum: photon momentum at the relevant scales is negligible
     * (1 MJ / c ≈ 3.3e-3 kg·m/s). */
    double E_radiated_J;

    /* Reaction mass (propulsion burns): inertial-frame kinetic
     * energy and momentum of ejected propellant. See FRAME
     * CONVENTION at the top of this file. */
    double E_reaction_mass_J;
    double p_reaction_x_kgmps;
    double p_reaction_y_kgmps;
    double p_reaction_z_kgmps;

    /* Chemistry energy released into the system by propulsion
     * burns. Subtracted from the closed-system sum to recover the
     * pre-burn baseline (the burn converted chemistry → KE, so the
     * post-burn KE exceeds the pre-burn KE by E_chem). Typically
     * ½ m_p v_e² (propellant rocket-frame KE at nozzle exit) per
     * burn. Zero in pure ballistic regimes (no burns). */
    double E_chem_source_J;

    /* Bookkeeping counters: number of recorded events of each kind. */
    size_t dispatch_count;
    size_t radiated_count;
    size_t reaction_mass_count;
    size_t chem_source_count;
} K26AstroRtConservationTotals;

/* Lifecycle ------------------------------------------------------- */

/* Allocate a zero-initialised ledger. Returns NULL on allocation
 * failure. */
K26AstroRtConservationLedger *k26astro_rt_ledger_new(void);

/* Free. Safe on NULL. */
void k26astro_rt_ledger_destroy(K26AstroRtConservationLedger *ledger);

/* Record slots ---------------------------------------------------- */

/* Record a solid-mass emission (any mass leaving the surviving-
 * vehicle population — ejecta, a released projectile, a jettisoned
 * component). Accumulates:
 *   E_ejecta_J       += ½ × mass_kg × (vx² + vy² + vz²)
 *   p_ejecta_*       += mass_kg × v_ejecta_*
 *   L_z_ejecta_kg_m2ps += mass_kg × (origin_x v_y − origin_y v_x)
 *
 * `origin_*` is the dispatch position in the same inertial frame
 * as the surviving-bodies subtotal. The angular-momentum term
 * closes the propagator's L_z balance when the projectile is
 * removed from integration (the propagator no longer sees its
 * contribution; the ledger carries it).
 *
 * Safe on NULL ledger (no-op). Negative mass / NaN velocity are
 * accepted as-is — the caller is responsible for emission validity. */
void k26astro_rt_ledger_record_dispatch(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_ejecta_x, double v_ejecta_y, double v_ejecta_z,
    double origin_x,   double origin_y,   double origin_z);

/* Record radiated energy in joules (absorbed-at-receiver or
 * coupled-to-receiver; the choice is the caller's, but should be
 * consistent across the run to compare against an integrator-drift
 * gate). */
void k26astro_rt_ledger_record_radiated(
    K26AstroRtConservationLedger *ledger,
    double energy_J);

#ifdef __GNUC__
#define K26_DEPRECATED_FRAME_CONVENTION \
    __attribute__((deprecated( \
        "Use k26astro_rt_ledger_record_burn (recommended) or " \
        "k26astro_rt_ledger_record_reaction_mass_inertial. " \
        "This legacy entry takes any-frame velocity and produces a " \
        "meaningful closed-system invariant only when fed the " \
        "propellant's inertial-frame velocity. See the FRAME " \
        "CONVENTION note in conservation_ledger.h.")))
#else
#define K26_DEPRECATED_FRAME_CONVENTION
#endif

/* Record a propulsion-burn reaction-mass dispatch (legacy entry).
 *
 * DEPRECATED: this entry accepts any-frame velocity and the math
 * is correct only when the caller passes the propellant's
 * inertial-frame velocity. Passing rocket-frame velocity here
 * leaves a closed-system drift on the order of m_p · v_vehicle · v_e
 * per burn. New code should use either:
 *
 *   - k26astro_rt_ledger_record_burn (ergonomic; takes vehicle
 *     velocity and rocket-frame exhaust speed/direction, computes
 *     both inertial-frame reaction mass and chemistry source
 *     atomically).
 *   - k26astro_rt_ledger_record_reaction_mass_inertial (low-level;
 *     takes the inertial-frame velocity directly).
 *
 * The legacy entry is retained for backward compatibility with
 * tests written before the convention was made explicit. */
K26_DEPRECATED_FRAME_CONVENTION
void k26astro_rt_ledger_record_reaction_mass(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_exhaust_x, double v_exhaust_y, double v_exhaust_z);

/* Record a propulsion-burn reaction-mass dispatch — inertial-frame
 * variant. Same accumulator math as the legacy entry, with the
 * frame requirement encoded in the function name.
 *
 * Accumulates:
 *   E_reaction_mass_J += ½ × mass_kg × |v_inertial|²
 *   p_reaction_*      += mass_kg × v_inertial_*
 *
 * v_inertial is the propellant's velocity in the same inertial
 * frame as the surviving-vehicles ke sum. For a vehicle at
 * inertial velocity v_v ejecting propellant at rocket-frame
 * exhaust speed v_e along direction -d (d = vehicle thrust
 * direction unit vector), the inertial exhaust velocity is
 * v_v - v_e * d. Most callers should prefer record_burn, which
 * does this computation internally. */
void k26astro_rt_ledger_record_reaction_mass_inertial(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_inertial_x, double v_inertial_y, double v_inertial_z);

/* Record chemistry energy released into the system by a
 * propulsion burn. Accumulates:
 *   E_chem_source_J += energy_J
 *
 * Typically energy_J = ½ m_p v_e² — the propellant's rocket-frame
 * kinetic energy at the nozzle exit, which equals the chemistry
 * energy released by the burn under the ideal-rocket assumption
 * that all chemistry energy converts to exhaust KE in the rocket
 * frame. For non-ideal engines (heat loss, ionisation, etc.) the
 * caller can supply a fraction-of-ideal value.
 *
 * The closed-system invariant subtracts E_chem_source from the
 * post-burn sum to recover the pre-burn baseline. */
void k26astro_rt_ledger_record_chemical_source(
    K26AstroRtConservationLedger *ledger,
    double energy_J);

/* Ergonomic helper: record a propulsion burn atomically from
 * physical inputs.
 *
 * Inputs:
 *   mass_kg           Propellant mass ejected (kg).
 *   v_vehicle_x/y/z   Vehicle velocity at burn epoch in the
 *                     inertial frame (same frame as the
 *                     surviving-vehicles ke sum).
 *   exhaust_speed_rocket_frame
 *                     Exhaust speed in the rocket frame (m/s),
 *                     i.e. v_e = Isp × g0 for a chemical engine.
 *   exhaust_dir_x/y/z Unit vector pointing in the exhaust
 *                     direction (opposite to the vehicle's
 *                     thrust direction). The function does not
 *                     re-normalise; the caller should pass a
 *                     unit vector.
 *
 * The helper computes:
 *
 *   v_p_inertial = v_vehicle + exhaust_speed_rocket_frame * exhaust_dir
 *   E_reaction_mass_J  += ½ m_p |v_p_inertial|²        (interpretation b)
 *   p_reaction_*       += m_p × v_p_inertial_*
 *   E_chem_source_J    += ½ m_p × exhaust_speed_rocket_frame²
 *
 * Both the reaction-mass and chemistry-source slots are updated in
 * a single call — there is no foot-gun where the caller records
 * one but forgets the other, and no opportunity to pass the wrong
 * frame at the call site. The exhaust_dir convention is the
 * physical exhaust direction in the inertial frame (typically
 * opposite to the vehicle thrust direction; for a vehicle thrusting
 * along +x to gain +x velocity, exhaust_dir = (-1, 0, 0)). */
void k26astro_rt_ledger_record_burn(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_vehicle_x, double v_vehicle_y, double v_vehicle_z,
    double exhaust_speed_rocket_frame,
    double exhaust_dir_x, double exhaust_dir_y, double exhaust_dir_z);

/* Accessors ------------------------------------------------------- */

/* Copy the running totals into `out`. Safe on NULL ledger (writes
 * zero-filled totals to `out`). Safe on NULL out (no-op). */
void k26astro_rt_ledger_totals(
    const K26AstroRtConservationLedger *ledger,
    K26AstroRtConservationTotals *out);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_CONSERVATION_LEDGER_H */
