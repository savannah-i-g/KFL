/* mass_step.c — k26astro_vehicle_commit_mass_step.
 *
 * Called by the integrator at the close of each gravity substep.
 * Applies any accumulated dot_m * dt to the vehicle's basic mass and
 * propagates the new mass to the bound body (which the gravity
 * integrator uses for GM at the next substep). The accumulator is
 * cleared so the next substep starts from zero.
 *
 * Engine-cluster thrust callbacks add to the accumulator via
 * k26astro_vehicle_mass_accum_add with a negative dot_m for mass
 * leaving the vehicle. Vehicles with no registered thrust callback
 * see a zero accumulator on every substep and this commit is
 * effectively a no-op.
 *
 * Inertia re-evaluation: if a stage event fired mid-substep, the
 * stage-event handler in stage_event.c has already updated the Ext
 * inertia tensor (and its inverse) before this function runs. No
 * additional work needed here.
 *
 * COM re-evaluation: ties to the propellant-distribution coupling
 * across registered tanks; that coupling is not currently modelled.
 * com_offset is treated as static here. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "vehicle_internal.h"

void k26astro_vehicle_commit_mass_step(K26AstroVehicle *v, double dt)
{
    if (!v) return;

    double dm = v->mass_accum * dt;
    if (dm != 0.0) {
        double new_mass = v->basic_mass_kg + dm;
        if (new_mass < 0.0) new_mass = 0.0;
        v->basic_mass_kg = new_mass;
        if (v->body) {
            k26astro_body_set_mass(v->body, new_mass);
        }
    }
    v->mass_accum = 0.0;
}
