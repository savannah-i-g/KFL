/* conservation_ledger.c — closed-system conservation accumulator.
 *
 * Pure POD accumulator; deterministic by construction. See
 * include/k26astro_rt/conservation_ledger.h for the contract. */
#include "k26astro_rt/conservation_ledger.h"

#include <stdlib.h>
#include <string.h>

struct K26AstroRtConservationLedger {
    K26AstroRtConservationTotals t;
};

K26AstroRtConservationLedger *k26astro_rt_ledger_new(void)
{
    K26AstroRtConservationLedger *ledger =
        (K26AstroRtConservationLedger *)calloc(1, sizeof(*ledger));
    return ledger; /* NULL on OOM */
}

void k26astro_rt_ledger_destroy(K26AstroRtConservationLedger *ledger)
{
    free(ledger);
}

void k26astro_rt_ledger_record_dispatch(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_ejecta_x, double v_ejecta_y, double v_ejecta_z,
    double origin_x,   double origin_y,   double origin_z)
{
    if (!ledger) return;
    const double v2 = v_ejecta_x * v_ejecta_x +
                      v_ejecta_y * v_ejecta_y +
                      v_ejecta_z * v_ejecta_z;
    ledger->t.E_ejecta_J        += 0.5 * mass_kg * v2;
    ledger->t.p_ejecta_x_kgmps  += mass_kg * v_ejecta_x;
    ledger->t.p_ejecta_y_kgmps  += mass_kg * v_ejecta_y;
    ledger->t.p_ejecta_z_kgmps  += mass_kg * v_ejecta_z;
    /* Angular momentum about z: L_z = m × (x v_y − y v_x).
     * Closes the propagator's L_z column when the projectile is
     * removed from the surviving-bodies population. */
    ledger->t.L_z_ejecta_kg_m2ps +=
        mass_kg * (origin_x * v_ejecta_y - origin_y * v_ejecta_x);
    ledger->t.dispatch_count    += 1;
    (void)origin_z; /* z-only L tracking; full 3-D is a future ext */
}

void k26astro_rt_ledger_record_radiated(
    K26AstroRtConservationLedger *ledger,
    double energy_J)
{
    if (!ledger) return;
    ledger->t.E_radiated_J  += energy_J;
    ledger->t.radiated_count += 1;
}

void k26astro_rt_ledger_record_reaction_mass(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_exhaust_x, double v_exhaust_y, double v_exhaust_z)
{
    if (!ledger) return;
    const double v2 = v_exhaust_x * v_exhaust_x +
                      v_exhaust_y * v_exhaust_y +
                      v_exhaust_z * v_exhaust_z;
    ledger->t.E_reaction_mass_J     += 0.5 * mass_kg * v2;
    ledger->t.p_reaction_x_kgmps    += mass_kg * v_exhaust_x;
    ledger->t.p_reaction_y_kgmps    += mass_kg * v_exhaust_y;
    ledger->t.p_reaction_z_kgmps    += mass_kg * v_exhaust_z;
    ledger->t.reaction_mass_count   += 1;
}

void k26astro_rt_ledger_record_reaction_mass_inertial(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_inertial_x, double v_inertial_y, double v_inertial_z)
{
    /* The accumulator math is identical to the legacy entry; what
     * the _inertial suffix changes is the contract — this entry
     * promises the caller is passing inertial-frame velocity. We
     * inline the accumulator here (rather than delegating to the
     * deprecated legacy function) so this entry's call sites do
     * not trip the legacy function's deprecation warning. */
    if (!ledger) return;
    const double v2 = v_inertial_x * v_inertial_x +
                      v_inertial_y * v_inertial_y +
                      v_inertial_z * v_inertial_z;
    ledger->t.E_reaction_mass_J     += 0.5 * mass_kg * v2;
    ledger->t.p_reaction_x_kgmps    += mass_kg * v_inertial_x;
    ledger->t.p_reaction_y_kgmps    += mass_kg * v_inertial_y;
    ledger->t.p_reaction_z_kgmps    += mass_kg * v_inertial_z;
    ledger->t.reaction_mass_count   += 1;
}

void k26astro_rt_ledger_record_chemical_source(
    K26AstroRtConservationLedger *ledger,
    double energy_J)
{
    if (!ledger) return;
    ledger->t.E_chem_source_J  += energy_J;
    ledger->t.chem_source_count += 1;
}

void k26astro_rt_ledger_record_burn(
    K26AstroRtConservationLedger *ledger,
    double mass_kg,
    double v_vehicle_x, double v_vehicle_y, double v_vehicle_z,
    double exhaust_speed_rocket_frame,
    double exhaust_dir_x, double exhaust_dir_y, double exhaust_dir_z)
{
    if (!ledger) return;
    /* v_p_inertial = v_vehicle + v_e × exhaust_dir
     * (exhaust_dir is the physical exhaust direction in inertial
     * frame, typically opposite to vehicle thrust direction). */
    const double v_p_x = v_vehicle_x +
                         exhaust_speed_rocket_frame * exhaust_dir_x;
    const double v_p_y = v_vehicle_y +
                         exhaust_speed_rocket_frame * exhaust_dir_y;
    const double v_p_z = v_vehicle_z +
                         exhaust_speed_rocket_frame * exhaust_dir_z;
    k26astro_rt_ledger_record_reaction_mass_inertial(
        ledger, mass_kg, v_p_x, v_p_y, v_p_z);
    /* Chemistry source: ½ m_p v_e² (rocket-frame propellant KE,
     * equals chemistry energy released under ideal-rocket
     * assumption). */
    const double E_chem = 0.5 * mass_kg *
        exhaust_speed_rocket_frame * exhaust_speed_rocket_frame;
    k26astro_rt_ledger_record_chemical_source(ledger, E_chem);
}

void k26astro_rt_ledger_totals(
    const K26AstroRtConservationLedger *ledger,
    K26AstroRtConservationTotals *out)
{
    if (!out) return;
    if (!ledger) {
        memset(out, 0, sizeof(*out));
        return;
    }
    *out = ledger->t;
}
