/* perturb_registry.c — K26AstroPerturbList registration + dispatch.
 *
 * Also holds weak (no-op) defaults for the three built-in
 * perturbations (J2 / SRP / GR PPN-1); the real impls live in
 * perturb_j2.c / perturb_srp.c / perturb_gr_ppn1.c and override
 * these symbols at link time. The weak defaults exist so the
 * accel_total dispatch in force_direct.c doesn't dangle if a
 * caller links only a subset of perturbations (e.g. for a quick
 * test that doesn't need GR). */
#include "k26astro_grav/perturb.h"

#include <stdlib.h>
#include <string.h>

int k26astro_grav_register_perturb(K26AstroGravState *state,
                                    K26AstroPerturbFn fn,
                                    void *ctx)
{
    if (!state || !fn) return K26ASTRO_E_NULL;
    if (!state->perturbs) {
        state->perturbs = calloc(1, sizeof(K26AstroPerturbList));
        if (!state->perturbs) return K26ASTRO_E_ALLOC;
    }
    K26AstroPerturbList *p = state->perturbs;
    if (p->count >= p->capacity) {
        int new_cap = p->capacity ? p->capacity * 2 : 4;
        K26AstroPerturbFn *new_fns  = realloc(p->fns,
            (size_t)new_cap * sizeof(K26AstroPerturbFn));
        void **new_ctxs            = realloc(p->ctxs,
            (size_t)new_cap * sizeof(void *));
        if (!new_fns || !new_ctxs) {
            free(new_fns);  free(new_ctxs);
            return K26ASTRO_E_ALLOC;
        }
        p->fns  = new_fns;
        p->ctxs = new_ctxs;
        p->capacity = new_cap;
    }
    p->fns [p->count] = fn;
    p->ctxs[p->count] = ctx;
    p->count++;
    return K26ASTRO_E_OK;
}
