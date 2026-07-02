/* k26astro_grav/ias15.h — Gauss-Radau IAS15 (Rein & Spiegel 2015).
 *
 * IAS15 is a 15th-order adaptive-step integrator built on Gauss-Radau
 * quadrature nodes. Predictor-corrector iteration converges to
 * near-machine-precision per step (~10⁻¹⁵ relative error). The
 * step-size controller preserves relative precision through close
 * encounters; the algorithm degrades gracefully rather than blowing
 * up at near-singular geometry.
 *
 * Reference: Rein & Spiegel, "IAS15: a fast, adaptive, high-order
 * integrator for gravitational dynamics, accurate to machine
 * precision over a billion orbits", MNRAS 446:1424, 2015.
 *
 * The substep coefficient table is the canonical Gauss-Radau set
 * (8 internal nodes plus boundaries); each coefficient is sourced
 * as a hex-literal double in ias15_coeffs.c for compiler-
 * independent bit-exactness. The predictor-corrector tolerance is
 * 1e-16 relative, matched verbatim against the REBOUND reference.
 */
#ifndef K26ASTRO_GRAV_IAS15_H
#define K26ASTRO_GRAV_IAS15_H

#include "k26astro_grav/grav.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* IAS15 carry-over: the b-coefficients of the Gauss-Radau expansion
 * (predictor seed for the next step), e-coefficients (predictor
 * error tracking), and substep scratch buffers. */
struct K26AstroIAS15Carry {
    /* Per-body G-coefficients across the 7 internal Gauss-Radau nodes. */
    K26V3 *b[7];   /* arrays of length n_bodies, indexed [k][i] */
    K26V3 *e[7];   /* predictor error tracking */
    K26V3 *g[7];   /* Gauss-Radau substep acceleration values */

    K26V3 *at0;    /* acceleration at the start of the step */
    K26V3 *r_sub;  /* position scratch */
    K26V3 *v_sub;  /* velocity scratch */
    K26V3 *a_sub;  /* acceleration scratch */

    double dt_proposed;   /* step-size controller's next proposal */
    int    capacity;
    int    initialised;   /* 0 until first step's b_coeffs settle */
};

/* Advance the body system by dt seconds via one or more IAS15
 * substeps. `dt` is the **requested span**; the Rein-Spiegel 2015
 * §4 step-size controller sub-divides it transparently:
 *
 *   dt_required = dt_try · (epsilon / |b_6|)^(1/7)
 *
 * Each substep is accepted if dt_required ≥ 0.5·dt_try, otherwise
 * the substep is rejected (positions unchanged) and retried with
 * dt_try = dt_required. The controller's next-substep hint persists
 * in state->ias15_dt_hint across calls so the integrator gradually
 * shrinks dt at close approach and expands after.
 *
 * On return, state->dt_last holds the actual final substep taken.
 * Diagnostic counter state->ias15_rejected_steps tracks cumulative
 * rejections. */
int k26astro_grav_step_ias15(K26AstroGravState *state, double dt);

/* Reset the IAS15 carry-over (clear b/e predictor history). Useful
 * after a non-conservative perturbation (e.g. SOI handoff)
 * that invalidates the smooth-trajectory assumption. */
void k26astro_grav_ias15_reset(K26AstroGravState *state);

/* Set the per-step relative error tolerance (epsilon). Default is
 * 1.0e-9 (Rein-Spiegel 2015 §4). Lower means smaller substeps and
 * higher accuracy. Reasonable range: 1e-12 .. 1e-6. */
void k26astro_grav_ias15_set_tol(K26AstroGravState *state, double eps);

/* Read the last completed substep size (seconds). 0 until the
 * first IAS15 step. */
double k26astro_grav_ias15_get_dt_last(const K26AstroGravState *state);

/* Cumulative rejected-substep counter since state init. */
uint32_t k26astro_grav_ias15_rejected_steps(const K26AstroGravState *state);

/* Set a per-call wall-time budget (seconds). 0.0 (default)
 * disables the budget; positive values cause k26astro_grav_step_ias15
 * to return K26ASTRO_E_TIME_BUDGET if a single call exceeds the budget.
 *
 * Deterministic runs MUST keep the budget at 0.0; wall-clock
 * measurement is platform-dependent and the early-exit is NOT a
 * determinism gate. The budget exists as a safety net for chaotic-
 * region stalls (e.g. the Burrau Pythagorean 3-4-5 close-encounter
 * cascade) where the adaptive controller would otherwise spend
 * minutes converging on a step at micrometric tolerance. */
void k26astro_grav_ias15_set_wall_budget(K26AstroGravState *state,
                                          double budget_s);

#ifdef __cplusplus
}
#endif

#endif
