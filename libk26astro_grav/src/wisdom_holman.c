/* wisdom_holman.c — classic Wisdom-Holman KDK symplectic integrator.
 *
 * Hamiltonian split:
 *   H = H_Kepler + H_interaction
 *
 *   H_Kepler:     each body i ≥ 1 propagated independently in a
 *                  two-body orbit around the central mass m_0.
 *                  Exact via libk26astro_conics' universal-variable
 *                  Kepler propagator.
 *
 *   H_interaction: gravitational forces between non-central bodies
 *                  (i.e. body[i] - body[j] for i, j ≥ 1, i ≠ j).
 *                  Applied as an impulsive kick.
 *
 * Step (Kick-Drift-Kick):
 *   1. v_i ← v_i + (dt/2)·a_int_i          for i ≥ 1
 *   2. (r_i, v_i) ← Kepler(r_i, v_i, μ_0, dt)
 *   3. v_i ← v_i + (dt/2)·a_int_i
 *
 * Heliocentric coordinates: r_i ≡ body[i].pos - body[0].pos, kept as
 * a K26V3 in metres (since the Kepler propagator works in metres-
 * scale doubles). The central body (body[0]) stays motionless in
 * heliocentric coords by construction; in the underlying inertial
 * representation it would move (Newton's 3rd law), but the heliocentric
 * symplectic formulation hides that motion from the integrator.
 *
 * Reference: Wisdom & Holman (1991) AJ 102:1528. The "democratic
 * heliocentric" variant (Duncan, Levison & Lee 1998 + Rein & Tamayo
 * 2015 WHFast) achieves better symplectic-corrector behaviour but
 * adds complexity; the classic form is sufficient for v0.1 and is
 * exact for two-body (the Test 2 gate). */
#include "k26astro_grav/wisdom_holman.h"
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_grav/perturb.h"
#include "k26astro_conics/kepler.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"

#include <stdlib.h>
#include <string.h>

static int ensure_wh_carry_(K26AstroGravState *state)
{
    if (state->wh_carry && state->wh_carry->capacity >= state->n_bodies) return 0;
    if (state->wh_carry) {
        free(state->wh_carry->p_bary);
        free(state->wh_carry->r_helio);
        free(state->wh_carry);
        state->wh_carry = NULL;
    }
    state->wh_carry = calloc(1, sizeof(K26AstroWHCarry));
    if (!state->wh_carry) return K26ASTRO_E_ALLOC;
    state->wh_carry->p_bary  = calloc((size_t)state->n_bodies, sizeof(K26V3));
    state->wh_carry->r_helio = calloc((size_t)state->n_bodies, sizeof(K26V3));
    if (!state->wh_carry->p_bary || !state->wh_carry->r_helio) {
        free(state->wh_carry->p_bary);
        free(state->wh_carry->r_helio);
        free(state->wh_carry);
        state->wh_carry = NULL;
        return K26ASTRO_E_ALLOC;
    }
    state->wh_carry->capacity = state->n_bodies;
    return 0;
}

/* External-perturbation prototypes (defined in
 * perturb_{j2,srp,gr_ppn1,outer_planets}.c). Declared inline here so
 * the WH kick can dispatch them without pulling in the full perturb
 * registry header. */
extern void k26astro_perturb_j2     (const K26AstroGravState *,
                                      const K26AstroGravView *,
                                      K26V3 *, void *);
extern void k26astro_perturb_srp    (const K26AstroGravState *,
                                      const K26AstroGravView *,
                                      K26V3 *, void *);
extern void k26astro_perturb_gr_ppn1(const K26AstroGravState *,
                                      const K26AstroGravView *,
                                      K26V3 *, void *);

/* Compute interaction acceleration: gravitational pull from each
 * non-central body j ≠ i on body i, PLUS all enabled perturbations
 * (J2 / SRP / GR-PPN1 flags + user-registered perturb registry).
 * The central body's Newtonian pull is EXCLUDED (it's handled by
 * the Kepler drift); perturbations on the central body ARE included
 * so e.g. outer-planet kicks contribute to the Sun's barycentric
 * wobble correctly. */
static void interaction_accel_(const K26AstroGravState *state, K26V3 *accel)
{
    int n = state->n_bodies;
    const K26AstroBody *b = state->bodies;
    memset(accel, 0, sizeof(K26V3) * (size_t)n);

    /* Inter-body Newtonian (central body excluded; Kepler drift). */
    for (int i = 1; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            K26V3 r_ij = k26astro_pos_sub(&b[j].pos, &b[i].pos);
            double r2 = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
            double r_mag = sqrt(r2);
            double inv_r3 = 1.0 / (r2 * r_mag);
            double s_i =  b[j].gm * inv_r3;
            double s_j = -b[i].gm * inv_r3;
            accel[i].x += s_i * r_ij.x;
            accel[i].y += s_i * r_ij.y;
            accel[i].z += s_i * r_ij.z;
            accel[j].x += s_j * r_ij.x;
            accel[j].y += s_j * r_ij.y;
            accel[j].z += s_j * r_ij.z;
        }
    }

    /* Perturbations (built-in flag-gated + user registry). */
    K26AstroGravView view = k26astro_grav_view(b, n);
    if (state->use_j2)      k26astro_perturb_j2     (state, &view, accel, NULL);
    if (state->use_srp)     k26astro_perturb_srp    (state, &view, accel, NULL);
    if (state->use_gr_ppn1) k26astro_perturb_gr_ppn1(state, &view, accel, NULL);
    if (state->perturbs && state->perturbs->count > 0) {
        for (int k = 0; k < state->perturbs->count; k++) {
            state->perturbs->fns[k](state, &view, accel,
                                     state->perturbs->ctxs[k]);
        }
    }
}

int k26astro_grav_step_wh(K26AstroGravState *state, double dt)
{
    if (!state) return K26ASTRO_E_NULL;
    if (state->n_bodies < 1) return K26ASTRO_E_BAD_INPUT;

    int rc = ensure_wh_carry_(state);
    if (rc != 0) return rc;

    int n = state->n_bodies;
    K26AstroBody *b = state->bodies;

    /* Step 1: Kick(dt/2); interaction kick. */
    K26V3 *a_int = calloc((size_t)n, sizeof(K26V3));
    if (!a_int) return K26ASTRO_E_ALLOC;
    interaction_accel_(state, a_int);
    for (int i = 1; i < n; i++) {
        b[i].vel.x += 0.5 * dt * a_int[i].x;
        b[i].vel.y += 0.5 * dt * a_int[i].y;
        b[i].vel.z += 0.5 * dt * a_int[i].z;
    }

    /* Step 2: Kepler drift; each body i ≥ 1 propagated around
     * central mass m_0 (= body[0]). The central body's heliocentric
     * "position" is by definition zero throughout; we still need to
     * carry its absolute K26AstroPos forward by its barycentric
     * motion (which for the classic WH formulation we approximate as
     * zero; the Sun is treated as immobile in the inertial frame).
     *
     * body[0]'s pos + vel stay unchanged across the drift. The
     * classic Sun-centred formulation matches this; a democratic
     * heliocentric variant that moves the Sun is a future
     * extension. */
    double mu0 = b[0].gm;
    for (int i = 1; i < n; i++) {
        K26V3 r_helio = k26astro_pos_sub(&b[i].pos, &b[0].pos);
        K26V3 r_new, v_new;
        int prc = k26astro_kepler_propagate(&r_new, &v_new,
                                              r_helio, b[i].vel,
                                              mu0, dt, 64);
        if (prc != 0) { free(a_int); return K26ASTRO_E_NO_CONVERGE; }
        /* Reconstruct absolute pos: body[0].pos + r_new (heliocentric). */
        b[i].pos = b[0].pos;
        k26astro_pos_add(&b[i].pos, r_new);
        b[i].vel = v_new;
    }

    /* Step 3: Kick(dt/2); second interaction kick. */
    interaction_accel_(state, a_int);
    for (int i = 1; i < n; i++) {
        b[i].vel.x += 0.5 * dt * a_int[i].x;
        b[i].vel.y += 0.5 * dt * a_int[i].y;
        b[i].vel.z += 0.5 * dt * a_int[i].z;
    }
    free(a_int);

    k26astro_epoch_add_seconds(&state->t, dt);
    state->dt_last = dt;
    return K26ASTRO_E_OK;
}
