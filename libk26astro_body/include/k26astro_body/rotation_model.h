/* k26astro_body/rotation_model.h — IAU 2018 rotation model evaluator.
 *
 * The IAU Working Group on Cartographic Coordinates and Rotational
 * Elements publishes pole position (α₀, δ₀) and prime meridian
 * angle (W) as polynomial functions of time, plus body-specific
 * trigonometric corrections. K26 uses these to compute body-fixed
 * attitude (and hence the ICRF↔body-fixed frame transform) for
 * planets, moons, and dwarf planets.
 *
 * Sources:
 *   - Archinal et al. (2018), "Report of the IAU Working Group on
 *     Cartographic Coordinates and Rotational Elements: 2015",
 *     Celestial Mechanics and Dynamical Astronomy 130:22
 *   - Archinal et al. (2019), erratum addressing Phobos prime
 *     meridian + Figure 1/2 sign corrections
 *
 * Coefficient tables are in iau_rotations.c — one entry per named
 * body. Adding new bodies (Saturnian moons not yet covered, future
 * mission targets, asteroids with measured rotation) is a clean
 * append to that table.
 *
 * The evaluator returns a quaternion that takes inertial-frame
 * (ICRF) vectors into body-fixed coordinates:
 *
 *   v_body = q_attitude * v_icrf * q_attitude^-1
 *
 * The composition follows IAU convention:
 *   R = R3(W) · R1(90° - δ₀) · R3(90° + α₀)
 * where R1, R3 are rotations about the x and z axes.
 */
#ifndef K26ASTRO_BODY_ROTATION_MODEL_H
#define K26ASTRO_BODY_ROTATION_MODEL_H

#include <stdint.h>

#include "k26astro_core/epoch.h"
#include "k26m3d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Correction term ----------------------------------------- *
 *
 * Each term contributes A * trig(B + C * T) to one of α₀, δ₀, or W.
 * `arg_index` is a 1-based reference into the body's argument
 * table — many bodies share the same linear functions of T (Mars
 * shares M1..M14 across α, δ, W, etc.), so the table is encoded
 * compactly.
 *
 * Field shapes:
 *   target_idx   0 = α₀, 1 = δ₀, 2 = W
 *   trig_kind    0 = sin, 1 = cos
 *   amplitude    deg (for α / δ) or deg (for W); applied to the
 *                target after the linear part
 *   arg_a        constant part of the argument (deg)
 *   arg_b        linear part (deg/century for α/δ args; varies for
 *                W args — see body comment in iau_rotations.c)
 */
typedef struct {
    uint8_t  target_idx;
    uint8_t  trig_kind;
    double   amplitude;
    double   arg_a;
    double   arg_b;
} K26AstroIAURotCorrection;

/* ---- Rotation model entry ---------------------------------- */

#define K26ASTRO_IAU_MAX_CORRECTIONS 32

typedef struct {
    const char *name;              /* e.g. "iau2018:earth", "iau2018:io" */
    int         naif_id;           /* canonical NAIF id, for cross-ref */

    /* Linear terms — angles in degrees. */
    double      alpha_0;           /* pole RA at J2000 */
    double      delta_0;           /* pole Dec at J2000 */
    double      W_0;               /* prime meridian at J2000 */

    /* Rates — IAU convention: α/δ rates are per Julian century,
     * W rate is per day. */
    double      alpha_rate_per_century;
    double      delta_rate_per_century;
    double      W_rate_per_day;

    int                          n_corrections;
    K26AstroIAURotCorrection     corrections[K26ASTRO_IAU_MAX_CORRECTIONS];
} K26AstroIAURotation;

/* ---- Lookup ------------------------------------------------- */

/* Find a rotation model by name (e.g. "iau2018:earth"). Returns
 * NULL if not registered. */
const K26AstroIAURotation *k26astro_rotation_lookup(const char *name);

/* Find by NAIF id; useful when bridging from ephemeris queries. */
const K26AstroIAURotation *k26astro_rotation_by_naif(int naif_id);

/* Iterate the full table (NULL-terminator at end). */
const K26AstroIAURotation *k26astro_rotation_table(void);

/* ---- Evaluation -------------------------------------------- */

/* Evaluate the IAU rotation model at epoch `t` (any scale —
 * internally converted to TDB). Writes α, δ, W in degrees to the
 * output pointers; any may be NULL. */
void k26astro_rotation_eval_angles(const K26AstroIAURotation *r,
                                   const K26AstroEpoch *t,
                                   double *out_alpha_deg,
                                   double *out_delta_deg,
                                   double *out_W_deg);

/* Evaluate as a body-fixed-from-ICRF quaternion at epoch `t`. */
K26Quat k26astro_rotation_quaternion(const K26AstroIAURotation *r,
                                     const K26AstroEpoch *t);

/* Evaluate the 3x3 rotation matrix (ICRF → body-fixed). Useful for
 * direct vector multiply when the caller already works in matrices.
 * Output is column-major (same convention as libk26m3d K26M4 minus
 * the translation row). */
void k26astro_rotation_matrix(const K26AstroIAURotation *r,
                              const K26AstroEpoch *t,
                              double out_R[9]);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_BODY_ROTATION_MODEL_H */
