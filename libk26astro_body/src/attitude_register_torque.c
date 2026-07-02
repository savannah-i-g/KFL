/* attitude_register_torque.c — torque-source registry for
 * K26AstroAttitudeState.
 *
 * Mirrors the perturbation-list pattern in libk26astro_grav
 * (perturb_registry.c): lazy heap allocation, doubling realloc from
 * an initial capacity of 4, additive function-pointer table.
 *
 * Each registered function adds its body-frame torque contribution
 * to a caller-supplied accumulator; the registry is dispatched in
 * registration order. k26astro_attitude_step_torque_registry sums
 * the registry into a single body-frame torque and forwards to the
 * existing single-torque step.
 *
 * Reference: Markley & Crassidis, Fundamentals of Spacecraft
 * Attitude Determination and Control (2014), §3.10 — torque
 * superposition. */
#include "k26astro_body/attitude.h"

#include <stdlib.h>

struct K26AstroTorqueList {
    K26AstroTorqueFn *fns;
    void            **ctxs;
    int               count;
    int               capacity;
};

int k26astro_attitude_register_torque(K26AstroAttitudeState *state,
                                      K26AstroTorqueFn fn, void *ctx)
{
    if (!state || !fn) return 1;
    if (!state->torques) {
        state->torques = calloc(1, sizeof(struct K26AstroTorqueList));
        if (!state->torques) return 2;
    }
    struct K26AstroTorqueList *p = state->torques;
    if (p->count >= p->capacity) {
        int new_cap = p->capacity ? p->capacity * 2 : 4;
        K26AstroTorqueFn *new_fns = realloc(p->fns,
            (size_t)new_cap * sizeof(K26AstroTorqueFn));
        void **new_ctxs = realloc(p->ctxs,
            (size_t)new_cap * sizeof(void *));
        if (!new_fns || !new_ctxs) {
            free(new_fns);
            free(new_ctxs);
            return 2;
        }
        p->fns  = new_fns;
        p->ctxs = new_ctxs;
        p->capacity = new_cap;
    }
    p->fns [p->count] = fn;
    p->ctxs[p->count] = ctx;
    p->count++;
    return 0;
}

void k26astro_attitude_clear_torques(K26AstroAttitudeState *state)
{
    if (!state || !state->torques) return;
    struct K26AstroTorqueList *p = state->torques;
    free(p->fns);
    free(p->ctxs);
    free(p);
    state->torques = NULL;
}

void k26astro_attitude_destroy(K26AstroAttitudeState *state)
{
    k26astro_attitude_clear_torques(state);
}

void k26astro_attitude_step_torque_registry(K26AstroAttitudeState *state,
                                            double t, double dt)
{
    if (!state) return;
    K26V3 sum = { 0.0, 0.0, 0.0 };
    if (state->torques) {
        struct K26AstroTorqueList *p = state->torques;
        for (int i = 0; i < p->count; i++) {
            p->fns[i](state, t, &sum, p->ctxs[i]);
        }
    }
    k26astro_attitude_step_torque(state, sum, dt);
}
