/* signature.h — emitter / scatterer signature math.
 *
 * IR signatures derive from blackbody / greybody surface emission
 * (Stefan-Boltzmann: P = ε σ A T^4). The Planck-integrated
 * pass-band variant returns the in-band fraction of total emitted
 * power; the bolometric variant returns the total over all
 * wavelengths.
 *
 * Faceted RCS uses geometric-optics summation over facets visible
 * to the look vector. The model is appropriate when target
 * features are large compared to the wavelength; the resonant /
 * sub-wavelength regime is not modelled (stated limitation).
 *
 * References:
 *   - Hudson, R.D. 1969. Infrared System Engineering. Wiley.
 *   - Knott, Shaeffer, Tuley 2004. Radar Cross Section, 2nd ed.
 *   - Born, M. and Wolf, E. 2002. Principles of Optics, 7th ed.
 *
 * All functions are pure (no global state). Determinism: results
 * are bit-identical given identical inputs across reference CPUs
 * (no transcendental functions other than exp / pow which are
 * controlled by the lib's -ffp-contract=off / -fexcess-precision=
 * standard compilation flags).
 */
#ifndef K26ASTRO_DETECT_SIGNATURE_H
#define K26ASTRO_DETECT_SIGNATURE_H

#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- IR signatures -------------------------------------------- */

/* Total bolometric power radiated by a greybody surface:
 *     P = ε σ A T^4
 * Inputs: area_m2 >= 0, epsilon in [0, 1], T_K >= 0. */
double k26astro_signature_ir_bolometric_power(double area_m2,
                                              double epsilon,
                                              double T_K);

/* In-band Planck-integrated power. Integrates the Planck
 * spectral radiance over [lambda_lo_m, lambda_hi_m] using
 * Simpson's rule with n_quad subdivisions (recommended n_quad >= 64
 * for ~5e-5 relative error; the integrand is smooth). Returns
 * the in-band power in watts.
 *
 * The integrand is the Planck spectral radiance multiplied by 4π
 * (Lambertian integration over the full sphere). Caller may
 * choose to apply view-angle / orientation reduction post-hoc
 * for directional emission. */
double k26astro_signature_ir_planck_inband(double area_m2,
                                            double epsilon,
                                            double T_K,
                                            double lambda_lo_m,
                                            double lambda_hi_m,
                                            int n_quad);

/* Lambertian projected intensity (W / sr) from an emitting surface
 * at view-angle θ from the surface normal:
 *     I(θ) = (ε σ A T^4 / π) cos θ.
 * cos_view_angle must be in [0, 1]; negative values (back-facing
 * surface) return 0. */
double k26astro_signature_ir_lambertian(double area_m2,
                                        double epsilon,
                                        double T_K,
                                        double cos_view_angle);

/* ---- Radar cross-section -------------------------------------- */

/* Faceted RCS in m^2 from an array of body-frame facets, each
 * with an outward unit normal and a flat-plate area. The look
 * direction is in the same body frame and must be a unit vector.
 *
 * For each facet:
 *   - Compute cos θ = dot(normal, -look_unit) (look vector points
 *     from observer to target; emission points toward observer).
 *   - If cos θ <= 0, facet does not contribute (faces away).
 *   - Otherwise contribute σ_facet = (4π / λ^2) A^2 cos^4 θ for the
 *     monostatic flat-plate scattering in the geometric-optics
 *     regime (Knott 2004 eq. 5.45). The λ^2 division is folded
 *     into the unitless geometry factor here; the absolute RCS
 *     should be scaled by the caller's wavelength term to match
 *     the radar-equation convention.
 *
 * Returned value is the dimensionless cross-section geometry sum
 * Σ A_i^2 cos^4 θ_i (units m^4). The radar-equation consumer
 * multiplies by 4π / λ^2 to obtain RCS in m^2. */
double k26astro_signature_rcs_facet_geometry(int n_facets,
                                              const K26V3 *normals_body,
                                              const double *areas_m2,
                                              K26V3 look_unit_body);

/* Convenience: compute the conventional monostatic RCS in m^2
 * from the geometric sum plus wavelength. */
double k26astro_signature_rcs_monostatic(int n_facets,
                                          const K26V3 *normals_body,
                                          const double *areas_m2,
                                          K26V3 look_unit_body,
                                          double wavelength_m);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_SIGNATURE_H */
