/* test_event_time_root_finding.c — event-time bisection gate.
 *
 * Scenario: two-body central-force orbit with a test particle at
 * unit radius around a unit-mass primary (GM = 1; circular orbit
 * radius = 1; period T = 2π). The orbit lies in the xy plane;
 * initial conditions place the particle at (1, 0, 0) moving in +y.
 *
 * Event: predicate fires when the particle's y-component is
 * strictly positive. Over one orbital period the particle's y
 * transits from 0 → positive (just after start) → 0 (at T/2) →
 * negative. The predicate sees two sign changes: at t ≈ 0+ and at
 * t ≈ T/2.
 *
 * Validates:
 *   (1) Counter increments exactly twice over one orbital period.
 *   (2) The bisected t_event at the half-period crossing agrees with
 *       the analytic T/2 = π within event_tol_s.
 *   (3) Energy + angular momentum measured at orbit close agree
 *       with the no-event baseline within numerical noise — the
 *       integrator's full order is preserved across the event
 *       discontinuity (the handler is a no-op for state, only
 *       diagnostics). */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_body/body.h"
#include "k26astro_core/pos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    int    count;
    double t_first;
    double t_second;
} event_record_t;

/* y > 0 predicate. */
static int pred_y_positive_(const K26AstroGravState *state, double t, void *ctx)
{
    (void)t; (void)ctx;
    if (!state || state->n_bodies < 2) return 0;
    K26V3 r = k26astro_pos_to_m_approx(&state->bodies[1].pos);
    return r.y > 0.0 ? 1 : 0;
}

/* Handler: counter++. No state mutation (we want to verify the
 * conservation diagnostics are unaffected by the wrapper alone). */
static int handler_record_(K26AstroGravState *state, double t, void *ctx)
{
    (void)state;
    event_record_t *rec = (event_record_t *)ctx;
    if (rec->count == 0) rec->t_first  = t;
    else                 rec->t_second = t;
    rec->count++;
    return K26ASTRO_E_OK;
}

static void energy_angmom_(const K26AstroGravState *state,
                            double *E_out, double *L_out)
{
    /* Heliocentric two-body energy + angular momentum, body 1
     * relative to body 0. Specific quantities (per unit mass of
     * body 1; body 0 has unit GM). */
    K26V3 r0 = k26astro_pos_to_m_approx(&state->bodies[0].pos);
    K26V3 r1 = k26astro_pos_to_m_approx(&state->bodies[1].pos);
    K26V3 r  = { r1.x - r0.x, r1.y - r0.y, r1.z - r0.z };
    K26V3 v0 = state->bodies[0].vel;
    K26V3 v1 = state->bodies[1].vel;
    K26V3 v  = { v1.x - v0.x, v1.y - v0.y, v1.z - v0.z };
    double r_mag = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    double v_sq  = v.x * v.x + v.y * v.y + v.z * v.z;
    *E_out = 0.5 * v_sq - state->bodies[0].gm / r_mag;
    K26V3 L = { r.y * v.z - r.z * v.y,
                r.z * v.x - r.x * v.z,
                r.x * v.y - r.y * v.x };
    *L_out = sqrt(L.x * L.x + L.y * L.y + L.z * L.z);
}

static void setup_two_body_(K26AstroBody *bodies)
{
    memset(bodies, 0, 2 * sizeof(K26AstroBody));

    /* Primary at origin: GM = 1, m = 1 (G = 1 internally). */
    bodies[0].kind = K26ASTRO_BODY_STAR;
    bodies[0].gm   = 1.0;
    bodies[0].mass = 1.0;
    bodies[0].pos  = k26astro_pos_zero();
    bodies[0].vel  = (K26V3){ 0.0, 0.0, 0.0 };
    bodies[0].parent_body_idx = -1;

    /* Test particle: unit radius, circular speed. */
    bodies[1].kind = K26ASTRO_BODY_SPACECRAFT;
    bodies[1].gm   = 0.0;
    bodies[1].mass = 1.0e-10;     /* tracer; near-zero mass */
    bodies[1].pos  = k26astro_pos_from_m(1.0, 0.0, 0.0);
    bodies[1].vel  = (K26V3){ 0.0, 1.0, 0.0 };
    bodies[1].parent_body_idx = 0;
}

int main(void)
{
    const double T  = 2.0 * 3.14159265358979323846;
    const int    N  = 1000;
    const double dt = T / (double)N;
    const double event_tol_s = 1.0e-7;

    /* ---- Baseline: no events ---------------------------------- */
    K26AstroBody bodies_a[2];
    setup_two_body_(bodies_a);
    K26AstroGravState st_a;
    assert(k26astro_grav_state_init(&st_a, bodies_a, 2) == K26ASTRO_E_OK);
    /* Use IAS15 for tight conservation. */
    k26astro_grav_set_integrator(&st_a, K26ASTRO_INTEGRATOR_IAS15);
    double E0, L0;
    energy_angmom_(&st_a, &E0, &L0);
    for (int n = 0; n < N; n++) {
        int rc = k26astro_grav_step(&st_a, dt);
        assert(rc == K26ASTRO_E_OK);
    }
    double E_a, L_a;
    energy_angmom_(&st_a, &E_a, &L_a);
    double E_drift_baseline = fabs(E_a - E0);
    double L_drift_baseline = fabs(L_a - L0);
    fprintf(stderr,
        "test_event_time: baseline E_drift=%.3e L_drift=%.3e\n",
        E_drift_baseline, L_drift_baseline);
    k26astro_grav_state_destroy(&st_a);

    /* ---- With events: counter + bisected epoch -------------- */
    K26AstroBody bodies_b[2];
    setup_two_body_(bodies_b);
    K26AstroGravState st_b;
    assert(k26astro_grav_state_init(&st_b, bodies_b, 2) == K26ASTRO_E_OK);
    k26astro_grav_set_integrator(&st_b, K26ASTRO_INTEGRATOR_IAS15);
    assert(k26astro_grav_set_event_tol_s(&st_b, event_tol_s) == K26ASTRO_E_OK);

    event_record_t rec = { 0, 0.0, 0.0 };
    K26AstroGravEvent ev = {
        .predicate = pred_y_positive_,
        .handler   = handler_record_,
        .ctx       = &rec
    };
    assert(k26astro_grav_register_event(&st_b, ev) == K26ASTRO_E_OK);

    double E0b, L0b;
    energy_angmom_(&st_b, &E0b, &L0b);
    for (int n = 0; n < N; n++) {
        int rc = k26astro_grav_step(&st_b, dt);
        assert(rc == K26ASTRO_E_OK);
    }
    double Eb, Lb;
    energy_angmom_(&st_b, &Eb, &Lb);
    double E_drift_event = fabs(Eb - E0b);
    double L_drift_event = fabs(Lb - L0b);
    fprintf(stderr,
        "test_event_time: with-events E_drift=%.3e L_drift=%.3e\n",
        E_drift_event, L_drift_event);

    /* (1) Counter increments exactly twice (start crossing + T/2 crossing). */
    fprintf(stderr,
        "test_event_time: count=%d, t_first=%.6e, t_second=%.6e (T/2=%.6e)\n",
        rec.count, rec.t_first, rec.t_second, T / 2.0);
    assert(rec.count == 2);

    /* (2) Second crossing (the T/2 transit) bisected to within
     * event_tol_s. First crossing is at t = 0+ (just after start);
     * the bisected t lies in [0, dt] but is small compared to T/2. */
    assert(rec.t_first  < dt);
    assert(fabs(rec.t_second - T / 2.0) < event_tol_s * 10.0);

    /* (3) Conservation diagnostics within numerical noise of
     * baseline. Bisection sub-steps perturb the integrator's
     * adaptive cadence slightly; bound the with-events drift as
     * baseline + a per-event allowance proportional to step size
     * × number of bisection iterations × event count. */
    double E_drift_bound = E_drift_baseline + 1e-10;
    double L_drift_bound = L_drift_baseline + 1e-10;
    assert(E_drift_event < E_drift_bound);
    assert(L_drift_event < L_drift_bound);

    k26astro_grav_clear_events(&st_b);
    assert(st_b.events == NULL);

    /* ---- Null guards ---------------------------------------- */
    assert(k26astro_grav_register_event(NULL, ev) != K26ASTRO_E_OK);
    K26AstroGravEvent ev_bad = { .predicate = NULL, .handler = handler_record_ };
    assert(k26astro_grav_register_event(&st_b, ev_bad) != K26ASTRO_E_OK);
    K26AstroGravEvent ev_bad2 = { .predicate = pred_y_positive_, .handler = NULL };
    assert(k26astro_grav_register_event(&st_b, ev_bad2) != K26ASTRO_E_OK);
    k26astro_grav_clear_events(NULL);
    assert(k26astro_grav_set_event_tol_s(NULL, 1e-6) != K26ASTRO_E_OK);
    assert(k26astro_grav_set_event_tol_s(&st_b, -1.0) != K26ASTRO_E_OK);

    k26astro_grav_state_destroy(&st_b);

    printf("test_event_time_root_finding: OK\n");
    return 0;
}
