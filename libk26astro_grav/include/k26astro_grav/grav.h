/* k26astro_grav/grav.h — N-body integrator state + dispatch.
 *
 * libk26astro_grav owns the *numerical* integration of a body system.
 * Conic-style propagation lives in libk26astro_conics; the runtime
 * orchestration (close-encounter handoff, observer light-time,
 * snapshot format) lives in libk26astro_rt. This lib is the
 * numerical workhorse between them.
 *
 * K26AstroGravState carries:
 *   - the body array under integration
 *   - the chosen integrator (WH / IAS15 / Verlet / RK4 / RK45)
 *   - perturbation flags (J2, SRP, GR PPN-1) + user-registered perturbs
 *   - softening parameter (Plummer ε for direct N²)
 *   - integrator-specific carry-over (WH Jacobi cache, IAS15 G-coeffs)
 *   - current sim epoch
 *
 * libk26astro_rt's K26AstroWorld composes K26AstroGravState as a
 * member field (the KFL `world` opaque type's "grav" sub-state). The
 * KFL builtins in this lib's .kflbi bind to K26AstroGravState*; the
 * rt layer adds a thin forwarding shim so user KFL code sees them as
 * world-level operations.
 *
 * Determinism: in portable mode, force summation switches to
 * compensated (libk26astro_core/sum.h) when n_bodies ≥
 * K26ASTRO_COMPENSATED_SUM_THRESHOLD. FPU state pinned at
 * k26astro_grav_state_init (fpu_pin.c). Coefficient tables sourced
 * as hex-literal doubles (ias15_coeffs.c).
 */
#ifndef K26ASTRO_GRAV_H
#define K26ASTRO_GRAV_H

#include "k26astro_body/body.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/sum.h"
#include "k26m3d.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Integrator selection --------------------------------------- */

typedef enum {
    K26ASTRO_INTEGRATOR_WH        = 1,   /* Wisdom-Holman (default) */
    K26ASTRO_INTEGRATOR_IAS15     = 2,   /* Gauss-Radau adaptive */
    K26ASTRO_INTEGRATOR_VERLET    = 3,   /* DKD Velocity Verlet */
    K26ASTRO_INTEGRATOR_RK4       = 4,   /* fixed-step classical RK4 */
    K26ASTRO_INTEGRATOR_RK45      = 5,   /* Dormand-Prince via libk26compute */
    K26ASTRO_INTEGRATOR_MERCURIUS = 6    /* hybrid symplectic; orchestration in libk26astro_rt */
} K26AstroIntegrator;

/* ---- Perturbation context --------------------------------------- */

/* Forward declarations — registry shapes in perturb.h / forces.h. */
typedef struct K26AstroPerturbList    K26AstroPerturbList;
typedef struct K26AstroEventList      K26AstroEventList;
typedef struct K26AstroIAS15Carry     K26AstroIAS15Carry;
typedef struct K26AstroWHCarry        K26AstroWHCarry;
typedef struct K26AstroMercuriusContext K26AstroMercuriusContext;

/* ---- The state struct ------------------------------------------- */

typedef struct K26AstroGravState {
    K26AstroBody  *bodies;
    int            n_bodies;
    K26AstroIntegrator integrator;
    double         softening;            /* Plummer ε, metres */
    K26AstroEpoch  t;                    /* current sim epoch */
    double         dt_last;              /* last step taken (s) */

    uint8_t        use_j2;
    uint8_t        use_srp;
    uint8_t        use_gr_ppn1;

    K26AstroPerturbList *perturbs;       /* heap-owned registry */
    K26AstroEventList   *events;         /* heap-owned event registry */

    /* Event-time root-finding snapshot. Lazy-allocated by the
     * advance-with-events wrapper; grown if n_bodies exceeds capacity.
     * Reused across event-bearing steps within one state's lifetime. */
    K26AstroBody        *event_snapshot;
    int                  event_snapshot_cap;

    /* Tolerance for the event-time bisection (seconds). Default
     * 1e-6 s; clamp to >= 0. */
    double               event_tol_s;

    K26AstroIAS15Carry  *ias15_carry;    /* heap; lazy init on first IAS15 step */
    K26AstroWHCarry     *wh_carry;       /* heap; lazy init on first WH step */

    /* IAS15 step-size controller (Rein-Spiegel 2015 §4). */
    double   ias15_tol;            /* epsilon; per-step relative error budget */
    double   ias15_dt_last;        /* actual final substep taken */
    double   ias15_dt_hint;        /* controller-suggested next substep (carry) */
    uint32_t ias15_rejected_steps; /* diagnostic counter (cumulative) */

    /* IAS15 last-substep PC convergence diagnostics (Rein-Spiegel
     * 2015 §4 eq. 9). Refreshed at every accepted substep; left
     * at the rejected attempt's values if the most recent substep
     * was a reject (the controller's response to the reject is
     * captured by ias15_dt_hint). The defaults (0 iterations,
     * 0.0 eps_b) hold before the first IAS15 substep runs.
     *
     * Drivers wanting to emit a per-tick integrator diagnostic
     * read these AFTER calling k26astro_grav_step. A typical
     * healthy run shows pc_iter ∈ [2, 6] outside close encounters
     * and rises toward 12 inside them; 12 across every step is
     * the cap-and-quit pattern, indicating either too-tight
     * ias15_tol or genuinely stiff dynamics. */
    int      ias15_last_pc_iterations;  /* PC iter count on most recent attempt */
    double   ias15_last_eps_b_achieved; /* max_b6 / max_a0 on most recent attempt */

    /* IAS15 explicit snapshot/rollback. On reject the
     * b-coefficient carry is reused as predictor seed for the
     * smaller dt, but the body array is also byte-restored from
     * this snapshot so rejected-step state is bit-identical to
     * pre-step. Lazy-allocated on first ias15 substep; reused
     * thereafter. Reallocated if n_bodies grows past cap. */
    K26AstroBody *ias15_snapshot;  /* heap; size = ias15_snapshot_cap */
    int           ias15_snapshot_cap;

    /* IAS15 wall-time budget. Safety net for chaotic-region
     * stalls; default is 0.0 (no budget). When set,
     * the controller checks CLOCK_MONOTONIC at each substep entry
     * and returns K26ASTRO_E_TIME_BUDGET if the cumulative wall
     * time on this call exceeds the budget. NOT a determinism
     * gate; wall-clock measurement is platform-dependent. Use
     * only for test/debugging safety nets. */
    double   ias15_wall_budget_s;

    /* MERCURIUS orchestration context (non-owning pointer). Set by
     * the rt-layer orbit-step driver for the duration of a single
     * WH-far or IAS15-near substep, then cleared. NULL = standard
     * full-force integration (no MERCURIUS active). */
    const K26AstroMercuriusContext *mercurius;
} K26AstroGravState;

/* ---- Error codes ------------------------------------------------ */

#define K26ASTRO_E_OK            0
#define K26ASTRO_E_NULL          1
#define K26ASTRO_E_BAD_INPUT     2
#define K26ASTRO_E_NOT_WIRED     3   /* integrator selected but not yet implemented */
#define K26ASTRO_E_NO_CONVERGE   4
#define K26ASTRO_E_ALLOC         5
#define K26ASTRO_E_TIME_BUDGET   6   /* IAS15 wall-budget exceeded */

/* ---- Lifecycle -------------------------------------------------- */

/* Initialise grav state over a caller-owned body array. The state
 * keeps a pointer to `bodies` — caller must keep it alive for the
 * state's lifetime. Defaults: integrator=WH, softening=0,
 * perturbations off. FPU mode pinned to FE_TONEAREST with DAZ/FTZ
 * explicitly cleared. */
int  k26astro_grav_state_init(K26AstroGravState *state,
                               K26AstroBody *bodies,
                               int n_bodies);

/* Release perturbation registry, IAS15 carry, WH carry. Does NOT
 * free the body array (caller-owned). */
void k26astro_grav_state_destroy(K26AstroGravState *state);

/* ---- Step + dispatch ------------------------------------------- */

/* Advance the body system by dt seconds using the selected integrator.
 * Returns K26ASTRO_E_OK or one of the K26ASTRO_E_* codes. */
int k26astro_grav_step(K26AstroGravState *state, double dt);

/* Returns an integrator-appropriate suggested step:
 *   WH / Verlet :  T_min / 20
 *   IAS15       :  T_min / 100 (adaptive will refine)
 *   RK*         :  T_min / 40
 * where T_min is the shortest two-body orbital period among the
 * bodies in the state. Returns 0 if no orbital period can be
 * inferred (e.g. fewer than 2 bodies or all unbound). */
double k26astro_grav_advise_step(const K26AstroGravState *state);

/* ---- Integrator selection + tuning ------------------------------ */

int k26astro_grav_set_integrator(K26AstroGravState *state,
                                  K26AstroIntegrator which);
int k26astro_grav_set_softening (K26AstroGravState *state,
                                  double softening_m);

int k26astro_grav_enable_j2     (K26AstroGravState *state, int enable);
int k26astro_grav_enable_srp    (K26AstroGravState *state, int enable);
int k26astro_grav_enable_gr_ppn1(K26AstroGravState *state, int enable);

/* ---- IAS15 substep diagnostics --------------------------------- *
 *
 * Drivers that emit a per-tick integrator-diagnostic row read
 * these accessors after every k26astro_grav_step call. All four
 * report state from the most recent substep attempt; NULL state
 * yields zero-valued returns.
 *
 * pc_iterations: PC iteration count (1..IAS15_MAX_PC_ITER). Always
 *   at the cap signals controller failing to converge to ias15_tol
 *   on every step — read alongside eps_b_achieved to disambiguate
 *   "tolerance too tight" from "dynamics demand it."
 * dt_last_taken: the actual dt of the most recent accepted substep
 *   (seconds). Mirrors the ias15_dt_last field for direct access.
 * eps_b_achieved: the dimensionless error norm max|b_6|/max|a_0|
 *   from the most recent PC attempt. Compare against ias15_tol
 *   to diagnose convergence regime.
 * rejected_steps_total: cumulative count of rejected substeps
 *   since k26astro_grav_state_init. Mirrors ias15_rejected_steps. */
int      k26astro_grav_last_pc_iterations  (const K26AstroGravState *state);
double   k26astro_grav_last_dt_taken        (const K26AstroGravState *state);
double   k26astro_grav_last_eps_b_achieved  (const K26AstroGravState *state);
uint32_t k26astro_grav_rejected_steps_total (const K26AstroGravState *state);

#ifdef __cplusplus
}
#endif

#endif
