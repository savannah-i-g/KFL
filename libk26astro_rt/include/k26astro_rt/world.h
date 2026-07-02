/* libk26astro_rt — runtime world lifecycle + accessors.
 *
 * K26AstroWorld composes K26AstroGravState plus the
 * ephem/frame/conics machinery into a single steppable simulation.
 * The world is the KFL `world` opaque type. */
#ifndef K26ASTRO_RT_WORLD_H
#define K26ASTRO_RT_WORLD_H

#include <stdint.h>
#include <stddef.h>

#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"
#include "k26astro_core/frame.h"
#include "k26astro_body/body.h"
#include "k26astro_grav/grav.h"

#ifdef __cplusplus
extern "C" {
#endif

#define K26ASTRO_RT_LIB_VERSION "0.1.0"

/* World-level determinism mode. Set at create; immutable thereafter. */
typedef enum {
    K26ASTRO_MODE_FAST       = 0,   /* SIMD/FMA-tolerant, cosmetic */
    K26ASTRO_MODE_PORTABLE   = 1,   /* bit-identical cross-CPU */
    K26ASTRO_MODE_REFERENCED = 2    /* reserved; op-log variant */
} K26AstroWorldMode;

/* Coordinate-storage mode. See 02-coordinates.md. */
typedef enum {
    K26ASTRO_COORDS_SECTOR_GRID = 0,
    K26ASTRO_COORDS_Q64_64      = 1
} K26AstroCoordsMode;

/* Observer-correction mode. See 09-runtime.md §observer modes. */
typedef enum {
    K26ASTRO_OBS_GEOMETRIC   = 0,   /* no correction */
    K26ASTRO_OBS_ASTROMETRIC = 1,   /* light-time only (default) */
    K26ASTRO_OBS_APPARENT    = 2,   /* light-time + aberration [+ Shapiro] */
    K26ASTRO_OBS_TOPOCENTRIC = 3    /* reserved; requires libk26astro_atmos */
} K26AstroObserverMode;

/* Return codes. Mirrors libk26astro_grav conventions. */
enum {
    K26ASTRO_RT_OK              = 0,
    K26ASTRO_RT_E_NULL          = 1,
    K26ASTRO_RT_E_OOM           = 2,
    K26ASTRO_RT_E_BAD_ARG       = 3,
    K26ASTRO_RT_E_FULL          = 4,
    K26ASTRO_RT_E_NOT_FOUND     = 5,
    K26ASTRO_RT_E_FPU_RACE      = 6,
    K26ASTRO_RT_E_NOT_IMPLEMENTED = 7,
    K26ASTRO_RT_E_INTEGRATOR    = 8
};

/* Opaque. Layout in world_internal.h. */
typedef struct K26AstroWorld K26AstroWorld;

/* Body-iterator callback. Read-only; receives a const pointer to
 * the world's in-place body record. */
typedef void (*K26AstroBodyVisitorFn)(const K26AstroBody *b, int idx, void *user);

/* Forward declarations from libk26astro_ephem (avoid hard include of
 * the ephem header here; the lib pointer is opaque to the world). */
typedef struct K26AstroEphem K26AstroEphem;

/* Lifecycle ------------------------------------------------------- */

/* Create a world. The caller chooses determinism + coord mode at
 * construction; both are pinned for the world's lifetime. On
 * success returns a non-NULL handle and the process FPU state is
 * pinned to round-to-nearest + DAZ/FTZ cleared. Returns NULL on
 * allocation failure or invalid args.
 *
 * Multi-world processes: each `create` saves the caller's FPU
 * state and pins to the world's mode. Destroying restores. If two
 * worlds with different modes coexist live, step calls return
 * K26ASTRO_RT_E_FPU_RACE. Single-world is the primary supported
 * shape; multi-world with the same mode is safe. */
K26AstroWorld *k26astro_world_create(K26AstroWorldMode  mode,
                                      K26AstroCoordsMode coord_mode);

/* Tear down a world. Frees the body array, the scheduler, and
 * encounter state. Restores the caller's FPU state if this was the
 * last live world (or if its mode matched the saved state). */
void k26astro_world_destroy(K26AstroWorld *world);

/* Body collection -------------------------------------------------- */

/* Append a body. Returns its index in the world's array (>= 0) or
 * a negative K26ASTRO_RT_E_* code. The body is copied; the caller
 * does not need to keep `b` live. */
int  k26astro_world_add_body(K26AstroWorld *world, K26AstroBody b);

/* Remove a body by index. Tail-swap; later body indices may shift. */
int  k26astro_world_remove_body(K26AstroWorld *world, int idx);

/* Body count. */
int  k26astro_world_body_count(const K26AstroWorld *world);

/* Lookup a body by name. Returns its index or -1. */
int  k26astro_world_find_body(const K26AstroWorld *world, const char *name);

/* Read-only iterator. Stops if `visit` is NULL. */
void k26astro_world_for_each_body(const K26AstroWorld *world,
                                   K26AstroBodyVisitorFn visit,
                                   void *user);

/* Realloc-safe body handle. Returns a pointer into the world's body
 * array at `idx`, or NULL on out-of-range / NULL world. The returned
 * pointer is stable across step(), observe(), and remove_body() calls
 * (the body array is reallocated only by add_body). KFL builtin
 * `astro_world_body_at` surfaces this for the `body` opaque type. */
K26AstroBody *k26astro_world_body_at(K26AstroWorld *world, int idx);

/* Read the body's name (or NULL). Useful for KFL `for_each` body iter
 * to print / lookup the body by name. */
const char *k26astro_body_name(const K26AstroBody *b);

/* Stepping --------------------------------------------------------- */

/* Advance the world by wallclock_dt_s. Internally drives the
 * libk26tick scheduler which dispatches the orbit, spin, and render
 * channels at their registered rates. Returns 0 on success.
 *
 * Orbit substeps go through the MERCURIUS Rein-Tamayo close-encounter
 * handoff: pairs within `mercurius_outer` Hill-radius units have their
 * pairwise force smoothly partitioned between Wisdom-Holman (far)
 * and IAS15 (near) via the quintic switching function K(y). */
int  k26astro_world_step(K26AstroWorld *world, double wallclock_dt_s);

/* Single-body Kepler advance. Propagates `body_idx`'s state forward
 * by `dt` (seconds) on a Keplerian orbit around its parent_body_idx
 * (SOI parent). Falls back to body 0 (heliocentric convention) when
 * the body has no parent set.
 *
 * Use case: KFL `propagate <body> for <dt>`; a per-body shortcut for
 * the common case of advancing a single spacecraft along its current
 * conic without invoking the full N-body integrator. Cheap (a few
 * Newton-Raphson iterations on the universal anomaly χ); useful for
 * mission-design what-ifs, satellite preview, and ephemeris-style
 * spot checks.
 *
 * Does NOT advance the world clock or any other body; only mutates
 * the named body's pos + vel. The parent body is treated as a fixed
 * point during the step.
 *
 * Returns K26ASTRO_RT_OK on success, or:
 *   E_NULL       - null world
 *   E_BAD_ARG    - out-of-range body_idx, NaN dt, parent==body, or
 *                  parent has non-positive GM
 *   E_INTEGRATOR - Kepler propagator did not converge (rare; usually
 *                  near-parabolic geometry; see kepler_edge.h) */
int  k26astro_world_body_step(K26AstroWorld *world, int body_idx,
                              double dt);

/* Read the world's current epoch (orbit-channel clock). */
int  k26astro_world_now(const K26AstroWorld *world, K26AstroEpoch *out);

/* Configuration --------------------------------------------------- */

/* Attach an ephemeris for on-rails bodies (NULL clears). */
int  k26astro_world_set_ephem(K26AstroWorld *world, K26AstroEphem *ephem);

/* Observer-correction mode. Default ASTROMETRIC. */
int  k26astro_world_set_observer_mode(K26AstroWorld *world,
                                       K26AstroObserverMode mode);

/* Attach a libk26astro_atmos model (forward-declared opaque ptr;
 * the caller owns the lifetime; the world stores a non-owning
 * reference). Setting NULL detaches the atmosphere; TOPOCENTRIC
 * observer mode then degrades to APPARENT. */
struct K26AstroAtmos;
int  k26astro_world_set_atmos(K26AstroWorld *world,
                               struct K26AstroAtmos *atmos);

/* MERCURIUS K(y) transition factors. Default y_inner=3.0,
 * y_outer=5.0 (Rein-Tamayo 2019 §3.1). Both in Hill-radius units. */
int  k26astro_world_set_mercurius_factors(K26AstroWorld *world,
                                           double y_inner, double y_outer);

/* Direct access to the inner grav state (advanced API; integrators
 * + perturbation flags live here). The world owns this state; do
 * not destroy it directly. */
K26AstroGravState *k26astro_world_grav(K26AstroWorld *world);

/* Vehicle registry ------------------------------------------------- */

/* Opaque forward declaration. The concrete struct lives in
 * libk26astro_vehicle; the world holds a non-owning array of
 * pointers and treats vehicles as opaque. */
struct K26AstroVehicle;

/* Register a vehicle for substep-close mass-step commit. The world
 * does not take ownership; the caller retains lifetime. Re-registering
 * the same handle is a no-op (idempotent). Returns 0 on success,
 * negative K26ASTRO_RT_E_* on NULL or allocation failure. */
int  k26astro_world_register_vehicle  (K26AstroWorld *world,
                                       struct K26AstroVehicle *v);

/* Drop a vehicle from the registry. No-op on unknown handles or
 * NULL inputs. A vehicle's destructor should call this before its
 * own teardown to prevent the world dereferencing a freed pointer
 * on subsequent substep walks. */
void k26astro_world_unregister_vehicle(K26AstroWorld *world,
                                       struct K26AstroVehicle *v);

/* Install a non-gravitational acceleration cache for an inertial-
 * measurement-unit sensor. After this call, each evaluation of the
 * gravity state's perturbation registry snapshots accel_out[body_idx]
 * into the vehicle's last_non_grav_accel_inertial field. Inertial-
 * measurement-unit sensors read that field, rotate to body frame, and
 * report the felt acceleration.
 *
 * The capture callback is registered on the grav-state perturb list
 * at call time. Perturbation callbacks registered LATER will
 * contribute to the captured sum; perturbations registered EARLIER
 * will only contribute if their writes precede the capture in
 * registration order — by construction they do, since the capture
 * runs after all previously-registered callbacks.
 *
 * Calling this multiple times for the same (world, vehicle) pair
 * stacks additional capture callbacks; the last one wins. The
 * allocated context is tracked by the world and freed at destroy.
 *
 * Returns 0 on success or a negative K26ASTRO_RT_E_* code. */
int k26astro_world_enable_imu_accel_cache(K26AstroWorld *world,
                                          struct K26AstroVehicle *v,
                                          int body_idx);

/* Frame registry --------------------------------------------------- */

/* Register a user-frame. Returns a frame id >= K26A_FRAME_USER_BASE,
 * or a negative K26ASTRO_RT_E_* code (E_FULL if the per-world
 * registry is at capacity). */
int  k26astro_world_register_frame(K26AstroWorld *world,
                                    const K26AstroFrameInfo *info);

/* Look up a frame by name (case-insensitive). Built-in frames + any
 * registered with this world are searched. Returns the frame id or
 * a negative error code. */
int  k26astro_world_find_frame(const K26AstroWorld *world, const char *name);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_RT_WORLD_H */
