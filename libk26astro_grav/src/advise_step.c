/* advise_step.c — integrator-appropriate suggested step size.
 *
 *   WH / Verlet :  T_min / 20
 *   IAS15       :  T_min / 100  (the adaptive controller will refine)
 *   RK*         :  T_min / 40
 *
 * where T_min is the shortest two-body orbital period among the
 * bodies in the state, evaluated against the central body (body[0]).
 *
 * Returns 0 if no orbital period can be inferred (e.g. fewer than 2
 * bodies, or no body is bound to the central). */
#include "k26astro_grav/grav.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/consts.h"

#include <math.h>

double k26astro_grav_advise_step(const K26AstroGravState *state)
{
    if (!state || state->n_bodies < 2) return 0.0;
    const K26AstroBody *central = &state->bodies[0];
    if (central->gm <= 0.0) return 0.0;

    double T_min = INFINITY;
    for (int i = 1; i < state->n_bodies; i++) {
        const K26AstroBody *bi = &state->bodies[i];
        K26V3 r = k26astro_pos_sub(&bi->pos, &central->pos);
        double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        if (r_mag <= 0.0) continue;
        double v2 = (bi->vel.x - central->vel.x)*(bi->vel.x - central->vel.x)
                  + (bi->vel.y - central->vel.y)*(bi->vel.y - central->vel.y)
                  + (bi->vel.z - central->vel.z)*(bi->vel.z - central->vel.z);
        /* Semi-major axis from vis-viva: 1/a = 2/r - v²/μ */
        double inv_a = 2.0 / r_mag - v2 / central->gm;
        if (inv_a <= 0.0) continue;   /* unbound; skip */
        double a = 1.0 / inv_a;
        double T = K26A_TWO_PI * sqrt(a*a*a / central->gm);
        if (T < T_min) T_min = T;
    }

    if (!isfinite(T_min)) return 0.0;

    switch (state->integrator) {
    case K26ASTRO_INTEGRATOR_WH:
    case K26ASTRO_INTEGRATOR_VERLET:
        return T_min / 20.0;
    case K26ASTRO_INTEGRATOR_IAS15:
        return T_min / 100.0;
    case K26ASTRO_INTEGRATOR_RK4:
    case K26ASTRO_INTEGRATOR_RK45:
        return T_min / 40.0;
    case K26ASTRO_INTEGRATOR_MERCURIUS:
        return T_min / 20.0;
    }
    return T_min / 20.0;
}
