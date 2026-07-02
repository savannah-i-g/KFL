/* ias15_internal.h — private API for IAS15 implementation.
 *
 * Owned by libk26astro_grav; not exposed via include/k26astro_grav/.
 * The split: coefficients + scratch buffers in ias15_coeffs.c,
 * predictor-corrector loop in ias15_predictor.c, top-level driver
 * in ias15.c. */
#ifndef K26ASTRO_GRAV_IAS15_INTERNAL_H
#define K26ASTRO_GRAV_IAS15_INTERNAL_H

#include "k26astro_grav/grav.h"
#include "k26m3d.h"

/* Gauss-Radau 8-node table (h[0] = 0, h[1..7] = Radau nodes). */
extern const double k26_ias15_h[8];

/* Initialise c[], d[], r[] matrices from h[]. Called once at first
 * IAS15 step. Idempotent. */
void k26_ias15_init_matrices(void);

/* Read-back accessors (after init_matrices has been called). */
const double *k26_ias15_c(void);   /* 21 entries, lower-triangle */
const double *k26_ias15_d(void);   /* 21 entries, lower-triangle */
const double *k26_ias15_rr(void);  /* 28 entries, 1/(h_i - h_j) for i>j */

/* IAS15 carry init/free (called by ensure_carry_ helper). */
int  k26_ias15_carry_alloc  (K26AstroIAS15Carry **out, int n);
void k26_ias15_carry_release(K26AstroIAS15Carry *carry);

/* Predictor-corrector inner loop. Iterates until |Δb_6| / |b_6|
 * < tol; returns number of iterations used (1..max_iter) or -1 on
 * non-convergence.
 *
 * `out_max_b6`, if non-NULL, receives the converged max-over-bodies
 * of |b_6|, which the step-size controller uses to set the next
 * substep size. `out_max_a0`, if non-NULL, receives the converged
 * max-over-bodies of |a_0| (start-of-step acceleration magnitude). */
int k26_ias15_pc_iterate(K26AstroGravState *state, double dt,
                          int max_iter, double tol,
                          double *out_max_b6, double *out_max_a0);

#endif
