# Scientific library reference

KFL ships a suite of scientific libraries. A KFL program reaches a library
through **builtins** — named functions declared in the library's builtin
manifest (`.kflbi`) — and through **opaque handles**, values of a
library-defined type that the program holds and passes back to the library
(an ephemeris handle, a fitting problem, a simulation world). To use a
library, a program links it via `KFLC_LDLIBS` and the compiler must be
able to find its `.kflbi` manifest (installed by the library's `-dev`
package to `/usr/share/kflc/builtins/`, or pointed at with
`K26_KFL_BUILTIN_PATH`).

The tables below list each library's purpose and a representative slice of
its builtin surface; the authoritative list for any library is its
`.kflbi` manifest and the headers under its `include/` directory.

---

## Core — units, time, frames (`libk26astro_core`)

Foundational value types shared across the suite: positions, epochs,
reference-frame identifiers, and unit constructors.

| Builtin | Purpose |
|---|---|
| `epoch_j2000_tt` | The J2000.0 epoch in Terrestrial Time. |
| `epoch_convert`, `epoch_diff_seconds`, `epoch_to_jd` | Time-scale conversion, differences, Julian date. |
| `frame_by_name`, `frame_transform`, `frame_can_transform` | Reference-frame lookup and coordinate transforms. |
| `au`, `pc`, `ly`, `light_seconds`, `year`, `day`, `deg`, `arcsec` | Unit constructors for physical quantities. |
| `pos_from_m`, `pos_sub`, `pos_dist` | Position construction and geometry. |

## Ephemerides (`libk26astro_ephem`)

Planetary and lunar positions from a JPL SPK ephemeris. The kernel itself
(JPL DE441) is not distributed here; regenerate the inner-system subset
with the tool under `libk26astro_ephem/tools/`.

| Builtin | Purpose |
|---|---|
| `ephem_load`, `ephem_load_default`, `ephem_close` | Open/close an ephemeris (opaque `ephem`). |
| `ephem_lookup_name`, `ephem_id_name` | Body-name/identifier resolution. |
| `ephem_body_pos`, `ephem_body_state` | Position, and full position+velocity state, at an epoch. |
| `ephem_observe` | Apparent observation of a body from an observer. |

## Bodies and orbits (`libk26astro_body`)

Keplerian and equinoctial orbital elements, state-vector conversions,
anomaly transforms, rigid-body attitude propagation, and IAU rotation
models with Earth-orientation (precession/nutation/polar-motion) support.

| Builtin | Purpose |
|---|---|
| `elements_from_state`, `state_from_elements` | Convert between orbital elements and Cartesian state. |
| `equinoctial_from_kep`, `keplerian_from_eq` | Element-set conversions. |
| `anomaly_mean_to_true`, `anomaly_true_to_mean`, `orbital_period` | Anomaly transforms and period. |
| `attitude_step_free`, `attitude_step_torque`, `torque_gravity_gradient` | Rigid-body attitude dynamics. |
| `rotation_lookup`, `rotation_matrix`, `nutation_iau2000a` | IAU body rotation and Earth orientation. |

## Two-body conics (`libk26astro_conics`)

Closed-form two-body mechanics: universal-variable Kepler propagation,
Lambert's problem (single- and multi-revolution), and sphere-of-influence
geometry. This is the core of transfer and launch-window design — see
`kflc/examples/porkchop.kfl`.

| Builtin | Purpose |
|---|---|
| `kepler`, `kepler_from_elements` | Universal-variable Kepler propagation. |
| `lambert`, `lambert_multi` | Lambert boundary-value solver (single / multi-rev). |
| `soi_radius_hill`, `soi_radius_laplace`, `soi_crossing` | Sphere-of-influence radii and crossing detection. |

## Gravitation (`libk26astro_grav`)

N-body integrators (WHFast, IAS15 Gauss–Radau, velocity Verlet) with
optional perturbations: J2 oblateness, solar radiation pressure, GR PPN-1,
and outer-planet effects. Powers the `fn world` simulation surface via the
runtime.

| Builtin | Purpose |
|---|---|
| `step_world`, `grav_advise_step` | Advance a gravitational state; adaptive step advice. |
| `set_integrator`, `set_softening` | Integrator selection and softening. |
| `enable_j2`, `enable_srp`, `enable_gr_ppn1`, `enable_outer_planets` | Toggle perturbations. |
| `grav_register_perturb`, `grav_close_encounter` | Custom perturbations; close-encounter queries. |

## Simulation runtime (`libk26astro_rt`)

The `fn world` runtime manager. It composes gravitation, ephemeris,
frames, and conics into a single opaque `world`, provides a multi-rate
scheduler and the MERCURIUS hybrid integrator, computes light-time /
aberration-corrected observations, and reads and writes versioned binary
snapshots.

| Builtin | Purpose |
|---|---|
| `astro_world_open`, `astro_world_close` | Create/destroy a simulation world. |
| `astro_world_add_body`, `astro_world_find_body`, `astro_world_body_count` | Populate and query the world. |
| `astro_world_step`, `astro_world_observe` | Advance the world; observe a body. |
| `astro_world_snapshot_save`, `astro_world_snapshot_load` | Deterministic binary snapshots. |
| `astro_world_set_mercurius`, `astro_world_set_observer_mode` | Integrator and observer configuration. |

## Spacecraft composition (`libk26astro_vehicle`)

The vehicle composition root used by the runtime to model a spacecraft:
mass properties, engine clusters, propellant tanks, staging events, and
mass/centre-of-mass/inertia queries over time.

| Builtin | Purpose |
|---|---|
| `vehicle_new`, `vehicle_bind_body`, `vehicle_destroy` | Create a vehicle and bind it to a simulated body. |
| `vehicle_set_dry_mass`, `vehicle_set_inertia_diag`, `vehicle_set_com_offset` | Mass properties. |
| `vehicle_add_engine_cluster`, `vehicle_add_tank` | Propulsion composition. |
| `vehicle_schedule_stage_event`, `vehicle_schedule_active_window` | Staging and activity windows. |
| `vehicle_mass_at`, `vehicle_com_at`, `vehicle_inertia_at` | Time-dependent mass-property queries. |

Vehicle exposes further attachment points (thermal networks, sensors,
reaction-control quads, electrical-power buses, attitude controllers) whose
implementing libraries are distributed separately from this scientific
release.

## Curve fitting (`libk26astro_fit`)

Non-linear least-squares via a Levenberg–Marquardt solver (CMINPACK), with
ready-made residual models.

| Builtin | Purpose |
|---|---|
| `fit_problem_linear`, `fit_problem_quadratic`, `fit_problem_power`, `fit_problem_exp`, `fit_problem_gaussian` | Construct a residual problem for a standard model. |
| `fit_lmdif` | Run the Levenberg–Marquardt minimisation. |

## Numerical quadrature (`libk26astro_quad`)

Adaptive one-dimensional integration (QUADPACK) over canned integrand
families.

| Builtin | Purpose |
|---|---|
| `quad_integrand_poly2`, `quad_integrand_power`, `quad_integrand_exp`, `quad_integrand_gaussian`, `quad_integrand_planck` | Construct an integrand (opaque handle). |

## ODE integration (`libk26astro_ode`)

Stiff/non-stiff ODE integration (ODEPACK / LSODA) over standard
right-hand-side families.

| Builtin | Purpose |
|---|---|
| `ode_problem_harmonic`, `ode_problem_van_der_pol`, `ode_problem_lotka_volterra`, `ode_problem_robertson` | Construct an ODE problem. |
| `ode_lsoda_solve` | Integrate with automatic stiffness switching. |

## Geomagnetic field (`libk26astro_geomag`)

The IGRF-14 International Geomagnetic Reference Field.

| Builtin | Purpose |
|---|---|
| `geomag_field_v3` | Field vector at a geodetic position and epoch. |
| `geomag_field_magnitude` | Field magnitude. |

## Atmosphere (`libk26astro_atmos`)

The NRLMSISE-00 empirical atmosphere plus optical refraction and
in-scattering.

| Builtin | Purpose |
|---|---|
| `atmos_earth_standard`, `atmos_density_at` | Standard model; density at altitude. |
| `atmos_refraction_rad`, `atmos_apparent`, `atmos_inscatter` | Astronomical refraction and scattering. |

---

## Plotting and computation

The **plotting** engine (`libk26plot`) and the **computation** primitives
(`libk26compute`, `vector`/`matrix` and their operations) are reached
directly through the language rather than through builtins — `series_*`
statements, the `plot` declaration, and the vector/matrix constructors and
math functions. See [`getting-started.md`](getting-started.md) and
[`kflc/GRAMMAR.md`](../kflc/GRAMMAR.md).
