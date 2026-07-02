/* libk26astro_rt — private world layout.
 *
 * Not installed; lib internals only. */
#ifndef K26ASTRO_RT_WORLD_INTERNAL_H
#define K26ASTRO_RT_WORLD_INTERNAL_H

#include "k26astro_rt/world.h"
#include "k26tick.h"

#define K26ASTRO_WORLD_FRAME_MAX 64

#ifdef __cplusplus
extern "C" {
#endif

/* Per-world saved FPU state. Restored at destroy if this world
 * pinned the FPU. */
typedef struct {
    int  saved_round_mode;
    int  saved_daz;            /* x86 only */
    int  saved_ftz;            /* x86 only */
    int  applied;              /* 1 if pin was applied (so we restore) */
} K26AstroFPUState;

/* Per-pair close-encounter session. Allocated on demand in
 * encounter.c as pairs cross into the MERCURIUS transition region.
 *
 * `k_weight` is the Rein-Tamayo 2019 K(y) value computed at detect
 * time:
 *   K = 1 inside y_inner  → force entirely on the IAS15 "near" side
 *   K = 0 outside y_outer → force entirely on the WH "far" side
 *   0 < K < 1 in the transition window (C² smoothstep)
 * The MERCURIUS orchestrator (orbit_step.c) reads this weight in
 * both directions: the WH kick uses (1-K)·F, the IAS15 sub-step
 * uses K·F. */
typedef struct K26AstroEncounter {
    int       i, j;            /* body indices, i < j */
    double    y_last;          /* last r_ij / hill_radius value */
    double    k_weight;        /* K(y_last) per Rein-Tamayo 2019 */
    uint8_t   active;
} K26AstroEncounter;

struct K26AstroWorld {
    /* Composed integrator state — bodies live here. */
    K26AstroGravState grav;

    /* Optional ephemeris (NULL → no on-rails). */
    K26AstroEphem *ephem;

    /* Per-world frame registry (built-ins are inherited from
     * libk26astro_core; user-registered frames go here). */
    K26AstroFrameInfo frames[K26ASTRO_WORLD_FRAME_MAX];
    int               n_frames;
    int               next_user_id;   /* >= K26A_FRAME_USER_BASE */

    /* Multi-rate scheduler. */
    K26TickWorld   *tick;
    K26TickChannel  chan_orbit;
    K26TickChannel  chan_spin;
    K26TickChannel  chan_render;
    double          spin_hz;
    double          render_hz;

    /* MERCURIUS encounter session state. */
    K26AstroEncounter *encounters;
    int                n_encounters;
    int                cap_encounters;
    double             mercurius_hill_factor;   /* y_inner */
    double             mercurius_outer_factor;  /* y_outer */

    /* Observer-correction mode. */
    K26AstroObserverMode observer_mode;

    /* Atmospheric model (libk26astro_atmos). Caller-owned (non-
     * owning pointer). When NULL, TOPOCENTRIC observer mode
     * degrades to APPARENT (vacuum observer). Set via
     * k26astro_world_set_atmos. */
    struct K26AstroAtmos *atmos;

    /* Determinism + coord modes (immutable). */
    K26AstroWorldMode   mode;
    K26AstroCoordsMode  coord_mode;

    /* FPU state saved at create. */
    K26AstroFPUState    fpu;

    /* Snapshot metadata. 0 if never serialised. */
    uint32_t            snapshot_id;

    /* User context for callbacks. */
    void               *user;

    /* REFERENCED determinism mode op-log context. Opaque
     * K26AstroRefCtx allocated lazily by k26astro_world_set_ref_log_path
     * and freed by k26astro_world_ref_log_close (called from
     * world_destroy). Always NULL when mode != REFERENCED. */
    void               *ref_ctx;

    /* Vehicle registry for substep-close mass-step commit. Non-owning;
     * vehicles outlive worlds typically. Resized doubling from cap 4
     * like the vehicle slot arrays. orbit_step.c walks the registry
     * after k26astro_grav_step and calls k26astro_vehicle_commit_mass_step
     * on each entry to consume the per-substep mass accumulator that
     * propulsion thrust callbacks populate. */
    struct K26AstroVehicle **vehicles;
    int                      n_vehicles;
    int                      cap_vehicles;

    /* Inertial-measurement-unit capture contexts. Each entry is a
     * heap-allocated K26AstroWorldImuCacheCtx_ owned by the world
     * (the grav-state perturb list keeps a borrowed pointer). Freed
     * at world destroy. */
    void                   **imu_cache_ctxs;
    int                      n_imu_cache;
    int                      cap_imu_cache;
};

/* referenced.c — internal emit hooks. No-ops when mode is not
 * REFERENCED or the log path hasn't been set. */
void k26astro_rt_ref_emit_step_begin(struct K26AstroWorld *world, double dt);
void k26astro_rt_ref_emit_step_end  (struct K26AstroWorld *world,
                                      double world_time_s);
void k26astro_rt_ref_emit_encounter (struct K26AstroWorld *world,
                                      int i, int j, double k_weight);

/* fpu_state.c — save the current FPU mode, pin to deterministic mode.
 * Caller-owned struct so the save survives concurrent worlds. */
void k26astro_fpu_save_and_pin(K26AstroFPUState *fs);
void k26astro_fpu_restore     (const K26AstroFPUState *fs);

/* scheduler.c — register the three channels. Called from
 * world_create. */
int  k26astro_rt_scheduler_init(struct K26AstroWorld *world);
void k26astro_rt_scheduler_destroy(struct K26AstroWorld *world);

/* orbit_step.c — orbit-channel callback. Drives the encounter
 * detector + MERCURIUS handoff for `dt_s`. */
void k26astro_rt_orbit_step_cb(double dt_s, void *user);

/* The spin / render channels are stubs for now (the actual work
 * lives in libk26astro_render and in libk26astro_body's attitude
 * propagator). v0.1 callbacks no-op. */
void k26astro_rt_spin_step_cb  (double dt_s, void *user);
void k26astro_rt_render_step_cb(double dt_s, void *user);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_WORLD_INTERNAL_H */
