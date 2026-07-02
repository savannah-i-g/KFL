/* vehicle.c — K26AstroVehicle struct definition, constructor /
 * destructor, scalar setters, and queries. Slot bookkeeping is in
 * composition.c; stage events in stage_event.c; lifecycle windows
 * in lifecycle.c; the per-substep mass commit in mass_step.c. */

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/vehicle_consts.h"
#include "k26astro_body/attitude.h"
#include "k26astro_core/epoch.h"
#include "vehicle_internal.h"

#include <stdlib.h>
#include <string.h>

/* ---- Helpers --------------------------------------------------- */

static K26M3 identity_m3_(void)
{
    K26M3 I;
    memset(&I, 0, sizeof(I));
    I.m[0][0] = I.m[1][1] = I.m[2][2] = 1.0;
    return I;
}

/* ---- new / destroy -------------------------------------------- */

K26AstroVehicle *k26astro_vehicle_new(void)
{
    K26AstroVehicle *v = (K26AstroVehicle *)calloc(1, sizeof(*v));
    if (!v) return NULL;
    v->basic_mass_kg = K26ASTRO_VEHICLE_DEFAULT_BASIC_MASS_KG;
    v->mga_kg        = K26ASTRO_VEHICLE_DEFAULT_MGA_KG;
    /* Identity inertia keeps the Ext path's inverse well-defined
     * before the caller sets a real tensor. */
    k26astro_attitude_init_ext(&v->attitude_ext, identity_m3_());
    v->slot_attitude.live = 0;
    return v;
}

void k26astro_vehicle_destroy(K26AstroVehicle *v)
{
    if (!v) return;

    /* Notify every live subsystem slot so subsystem-side
     * back-references blank out before the vehicle storage
     * disappears. */
    k26astro_vehicle_notify_all_slots_(v);

    /* Bump generation so any subsystem cached snapshot can detect
     * the dead owner via k26astro_vehicle_generation. */
    v->generation++;

    /* Free per-event handler contexts that were heap-allocated by
     * k26astro_vehicle_register_stage_events_with. */
    if (v->event_ctxs) {
        for (int i = 0; i < v->n_ctxs; i++) {
            free(v->event_ctxs[i]);
        }
        free(v->event_ctxs);
    }

    free(v->slots_engine);
    free(v->slots_tank);
    free(v->slots_thermal);
    free(v->slots_rcs);
    free(v->slots_eps);
    free(v->slots_sensor);
    free(v->slots_payload);
    free(v->events);
    free(v->windows);

    k26astro_attitude_destroy_ext(&v->attitude_ext);
    /* body is non-owning. */
    free(v);
}

/* ---- bind_body ------------------------------------------------- */

void k26astro_vehicle_bind_body(K26AstroVehicle *v, K26AstroBody *b)
{
    if (!v) return;
    v->body = b;
}

/* ---- Mass setters --------------------------------------------- */

void k26astro_vehicle_set_dry_mass(K26AstroVehicle *v, double m)
{
    if (!v) return;
    v->basic_mass_kg = (m < 0.0) ? 0.0 : m;
}

void k26astro_vehicle_set_mga(K26AstroVehicle *v, double m)
{
    if (!v) return;
    v->mga_kg = (m < 0.0) ? 0.0 : m;
}

/* ---- Inertia setters ------------------------------------------ */

void k26astro_vehicle_set_inertia_diag(K26AstroVehicle *v,
                                       double ixx, double iyy, double izz)
{
    if (!v) return;
    if (ixx < 0.0 || iyy < 0.0 || izz < 0.0) return;
    K26M3 I;
    memset(&I, 0, sizeof(I));
    I.m[0][0] = ixx;
    I.m[1][1] = iyy;
    I.m[2][2] = izz;
    k26astro_attitude_update_inertia(&v->attitude_ext, I);
}

void k26astro_vehicle_set_inertia_full(K26AstroVehicle *v, K26M3 inertia)
{
    if (!v) return;
    k26astro_attitude_update_inertia(&v->attitude_ext, inertia);
}

/* ---- COM setter ----------------------------------------------- */

void k26astro_vehicle_set_com_offset(K26AstroVehicle *v,
                                     double x, double y, double z)
{
    if (!v) return;
    v->com_offset.x = x;
    v->com_offset.y = y;
    v->com_offset.z = z;
}

/* ---- Queries -------------------------------------------------- */

double k26astro_vehicle_mass_now(const K26AstroVehicle *v)
{
    if (!v) return 0.0;
    return v->basic_mass_kg;
}

double k26astro_vehicle_predicted_mass(const K26AstroVehicle *v)
{
    if (!v) return 0.0;
    return v->basic_mass_kg + v->mga_kg;
}

double k26astro_vehicle_mev_mass(const K26AstroVehicle *v)
{
    return k26astro_vehicle_predicted_mass(v);
}

double k26astro_vehicle_mass_at(const K26AstroVehicle *v, K26AstroEpoch t)
{
    if (!v) return 0.0;
    double m = v->basic_mass_kg;
    for (int i = 0; i < v->n_events; i++) {
        double dt_s = k26astro_epoch_diff_seconds(&t, &v->events[i].epoch);
        if (dt_s >= 0.0) {
            m -= v->events[i].mass_drop_kg;
            if (m < 0.0) m = 0.0;
        } else {
            /* events[] is sorted ascending — no further events fire. */
            break;
        }
    }
    return m;
}

K26V3 k26astro_vehicle_com_at(const K26AstroVehicle *v, K26AstroEpoch t)
{
    (void)t;
    K26V3 z = { 0.0, 0.0, 0.0 };
    if (!v) return z;
    return v->com_offset;
}

K26M3 k26astro_vehicle_inertia_at(const K26AstroVehicle *v, K26AstroEpoch t)
{
    if (!v) return identity_m3_();
    K26M3 I = v->attitude_ext.inertia;
    for (int i = 0; i < v->n_events; i++) {
        double dt_s = k26astro_epoch_diff_seconds(&t, &v->events[i].epoch);
        if (dt_s >= 0.0) {
            I.m[0][0] = v->events[i].ixx_after;
            I.m[1][1] = v->events[i].iyy_after;
            I.m[2][2] = v->events[i].izz_after;
        } else {
            break;
        }
    }
    return I;
}

/* ---- Internal accessors --------------------------------------- */

K26AstroAttitudeStateExt *k26astro_vehicle_attitude_ext(K26AstroVehicle *v)
{
    return v ? &v->attitude_ext : NULL;
}

K26AstroBody *k26astro_vehicle_body(K26AstroVehicle *v)
{
    return v ? v->body : NULL;
}

double k26astro_vehicle_mass_accum_get(const K26AstroVehicle *v)
{
    return v ? v->mass_accum : 0.0;
}

void k26astro_vehicle_mass_accum_add(K26AstroVehicle *v, double dot_m)
{
    if (v) v->mass_accum += dot_m;
}

void k26astro_vehicle_mass_accum_clear(K26AstroVehicle *v)
{
    if (v) v->mass_accum = 0.0;
}

uint64_t k26astro_vehicle_generation(const K26AstroVehicle *v)
{
    return v ? v->generation : 0;
}

K26V3 k26astro_vehicle_last_non_grav_accel_inertial(const K26AstroVehicle *v)
{
    K26V3 z = { 0.0, 0.0, 0.0 };
    if (!v) return z;
    return v->last_non_grav_accel_inertial;
}

void k26astro_vehicle_set_last_non_grav_accel_inertial(K26AstroVehicle *v,
                                                       K26V3 a)
{
    if (!v) return;
    v->last_non_grav_accel_inertial = a;
}
