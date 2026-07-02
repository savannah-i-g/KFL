/* k26astro_grav/wisdom_holman.h — WHFast (Rein & Tamayo 2015).
 *
 * The Wisdom-Holman mixed-variable symplectic integrator splits the
 * Hamiltonian into a dominant Keplerian part (each planet around the
 * central mass) and a perturbation kick from planet-planet
 * interactions:
 *
 *     Kick(dt/2) → Drift(dt via exact Kepler) → Kick(dt/2)
 *
 * The "democratic heliocentric" formulation (Duncan, Levison & Lee
 * 1998) keeps heliocentric coordinates but uses barycentric
 * momenta — this avoids the secular drift of the simpler "Jacobi"
 * variant and is what WHFast (Rein & Tamayo 2015) ships in REBOUND.
 *
 * Reference: Rein & Tamayo, "WHFast: a fast and unbiased
 * implementation of a symplectic Wisdom-Holman integrator for long-
 * term gravitational simulations", MNRAS 452:376, 2015.
 *
 * Energy oscillation bounded at second order in the perturbation
 * amplitude: for inner-solar-system runs ~10⁻⁶ peak-to-peak, no
 * secular drift. The Kepler-drift step uses
 * libk26astro_conics' k26astro_kepler_propagate.
 */
#ifndef K26ASTRO_GRAV_WISDOM_HOLMAN_H
#define K26ASTRO_GRAV_WISDOM_HOLMAN_H

#include "k26astro_grav/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

/* WH carry-over: cached barycentric momenta + Jacobi transform
 * working buffers. Allocated on first WH step; freed by
 * k26astro_grav_state_destroy. */
struct K26AstroWHCarry {
    K26V3 *p_bary;   /* barycentric momenta (n_bodies) */
    K26V3 *r_helio;  /* heliocentric position scratch */
    int    capacity;
};

/* Advance the body system by dt via one WH KDK substep. The first
 * body (state->bodies[0]) is the central mass; perturber bodies are
 * indices 1..n-1. */
int k26astro_grav_step_wh(K26AstroGravState *state, double dt);

#ifdef __cplusplus
}
#endif

#endif
