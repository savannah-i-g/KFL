/* test_infostate_smoke.c — construction, target push, latest
 * inspection, vehicle-slot binding, and destructor protocol.
 *
 * Validates that a fresh infostate accepts pushes, exposes the
 * latest sample correctly, integrates with the vehicle composition
 * slot, and tears down cleanly through both the infostate-first
 * and vehicle-first destruction orders. */

#include "k26astro_infostate/infostate.h"
#include "k26astro_infostate/infostate_consts.h"

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"
#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static K26AstroEpoch epoch_j2000_plus_s_(double dt_s)
{
    K26AstroEpoch e = k26astro_epoch_j2000_tt();
    k26astro_epoch_add_seconds(&e, dt_s);
    return e;
}

static K26AstroVehicle *make_vehicle_with_body_(K26AstroBody *b)
{
    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);
    k26astro_vehicle_bind_body(v, b);
    return v;
}

int main(void)
{
    /* ---- Construction with NULL observer fails. -------------- */
    assert(k26astro_infostate_new(NULL, 64) == NULL);

    /* ---- Construction with observer missing a body fails. ---- */
    K26AstroVehicle *v_no_body = k26astro_vehicle_new();
    assert(k26astro_infostate_new(v_no_body, 64) == NULL);
    k26astro_vehicle_destroy(v_no_body);

    /* ---- Happy-path construction + capacity report. ---------- */
    K26AstroBody obs_body;
    memset(&obs_body, 0, sizeof(obs_body));
    K26AstroVehicle *observer = make_vehicle_with_body_(&obs_body);

    K26AstroInfostate *s = k26astro_infostate_new(observer, 8);
    assert(s != NULL);
    assert(k26astro_infostate_observer(s) == observer);
    assert(k26astro_infostate_history_capacity(s) == 8);

    /* Minimum capacity clamps to 2. */
    K26AstroInfostate *s_small = k26astro_infostate_new(observer, 0);
    assert(s_small != NULL);
    assert(k26astro_infostate_history_capacity(s_small) ==
           K26ASTRO_INFOSTATE_MIN_HISTORY_CAPACITY);
    k26astro_infostate_destroy(s_small);

    /* ---- Target push + latest read --------------------------- */
    K26AstroBody tgt_body;
    memset(&tgt_body, 0, sizeof(tgt_body));
    K26AstroVehicle *target = make_vehicle_with_body_(&tgt_body);

    assert(k26astro_infostate_history_length(s, target) == 0);

    K26V3 r1 = { 1000.0, 0.0, 0.0 };
    K26V3 v1 = { 0.0,    0.0, 0.0 };
    k26astro_infostate_target_push(s, target, epoch_j2000_plus_s_(0.0), r1, v1);
    assert(k26astro_infostate_history_length(s, target) == 1);

    K26AstroInfostateObservation latest =
        k26astro_infostate_latest(s, target);
    assert(latest.valid == 1);
    assert(latest.position.x == 1000.0);

    /* Same-epoch push replaces in place. */
    K26V3 r1b = { 1500.0, 0.0, 0.0 };
    k26astro_infostate_target_push(s, target, epoch_j2000_plus_s_(0.0), r1b, v1);
    assert(k26astro_infostate_history_length(s, target) == 1);
    latest = k26astro_infostate_latest(s, target);
    assert(latest.position.x == 1500.0);

    /* Older-than-latest push silently drops. */
    K26V3 r0 = { 9999.0, 0.0, 0.0 };
    k26astro_infostate_target_push(s, target, epoch_j2000_plus_s_(-1.0), r0, v1);
    assert(k26astro_infostate_history_length(s, target) == 1);
    latest = k26astro_infostate_latest(s, target);
    assert(latest.position.x == 1500.0);

    /* Forward push extends the ring. */
    K26V3 r2 = { 2500.0, 0.0, 0.0 };
    k26astro_infostate_target_push(s, target, epoch_j2000_plus_s_(1.0), r2, v1);
    assert(k26astro_infostate_history_length(s, target) == 2);

    /* Ring overflow: push past capacity drops the oldest. */
    for (int i = 0; i < 20; i++) {
        K26V3 r = { 1000.0 * (double)(i + 3), 0.0, 0.0 };
        k26astro_infostate_target_push(s, target,
            epoch_j2000_plus_s_((double)(i + 2)), r, v1);
    }
    assert(k26astro_infostate_history_length(s, target) ==
           k26astro_infostate_history_capacity(s));

    /* ---- Vehicle-slot integration ----------------------------- */
    assert(k26astro_infostate_attach(s) == 0);

    /* Tearing down the infostate first should pull its own slot
     * off the vehicle. Vehicle subsequently destroyed with no
     * dangling reference. */
    k26astro_infostate_destroy(s);
    k26astro_vehicle_destroy(observer);

    /* ---- Tear-down second order: vehicle first. -------------- */
    K26AstroBody obs_body2;
    memset(&obs_body2, 0, sizeof(obs_body2));
    K26AstroVehicle *o2 = make_vehicle_with_body_(&obs_body2);
    K26AstroInfostate *s2 = k26astro_infostate_new(o2, 4);
    assert(k26astro_infostate_attach(s2) == 0);

    /* Push a target so the observe() call below reaches the
     * observer-NULL guard rather than the empty-track early-out. */
    k26astro_infostate_target_push(s2, target,
        epoch_j2000_plus_s_(0.0), r1, v1);

    /* Vehicle destroy fires the payload base's on_owner_destroy
     * callback which blanks s2->observer. Calling destroy on s2
     * afterwards must not dereference the freed vehicle. */
    k26astro_vehicle_destroy(o2);
    assert(k26astro_infostate_observer(s2) == NULL);

    /* observe() against a stale infostate must return valid=0, not
     * silently fall back to the world-origin observer (HIGH-3). */
    K26AstroInfostateObservation stale_obs =
        k26astro_infostate_observe(s2, target,
            epoch_j2000_plus_s_(0.0),
            K26ASTRO_INFOSTATE_MODALITY_NONE);
    assert(stale_obs.valid == 0);

    k26astro_infostate_destroy(s2);

    /* ---- NULL safety ----------------------------------------- */
    k26astro_infostate_destroy(NULL);
    k26astro_infostate_target_push(NULL, NULL,
        k26astro_epoch_j2000_tt(),
        (K26V3){ 0, 0, 0 }, (K26V3){ 0, 0, 0 });
    assert(k26astro_infostate_history_length(NULL, NULL) == 0);
    K26AstroInfostateObservation null_obs =
        k26astro_infostate_latest(NULL, NULL);
    assert(null_obs.valid == 0);

    printf("test_infostate_smoke: OK\n");
    return 0;
}
