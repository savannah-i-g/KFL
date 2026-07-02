/* rotation_iau.c — IAU 2018 rotation model evaluator.
 *
 * Polynomial + trigonometric corrections, then build attitude
 * quaternion from (α₀, δ₀, W) per the IAU convention:
 *
 *   R_ICRF→body = R3(W) · R1(90° - δ₀) · R3(90° + α₀)
 *
 * The 90° offsets reproduce the convention's axis choices:
 *   - +Z of the body frame points along the spin axis (right-hand
 *     rule from rotation direction, "north" by IAU convention).
 *   - +X is in the equatorial plane at the prime-meridian direction
 *     (= W = 0 at J2000, advances at W_rate_per_day).
 *
 * Coefficient tables are in iau_rotations.c; lookups go through
 * the iteration helpers here. */
#include "k26astro_body/rotation_model.h"
#include "k26astro_core/consts.h"

#include <math.h>
#include <string.h>
#include <strings.h>

/* Forward-declare the table (lives in iau_rotations.c). */
extern const K26AstroIAURotation K26ASTRO_IAU_ROTATION_TABLE[];

const K26AstroIAURotation *k26astro_rotation_table(void)
{
    return K26ASTRO_IAU_ROTATION_TABLE;
}

const K26AstroIAURotation *k26astro_rotation_lookup(const char *name)
{
    if (!name) return NULL;
    for (const K26AstroIAURotation *r = K26ASTRO_IAU_ROTATION_TABLE;
         r->name; r++) {
        if (strcasecmp(r->name, name) == 0) return r;
    }
    return NULL;
}

const K26AstroIAURotation *k26astro_rotation_by_naif(int naif_id)
{
    for (const K26AstroIAURotation *r = K26ASTRO_IAU_ROTATION_TABLE;
         r->name; r++) {
        if (r->naif_id == naif_id) return r;
    }
    return NULL;
}

/* ---- Epoch helpers ------------------------------------------- */
/* T_days  = days from J2000 in TDB
 * T_cent  = T_days / 36525 */

static void epoch_to_T_(const K26AstroEpoch *t,
                        double *out_days, double *out_centuries)
{
    K26AstroEpoch tdb = *t;
    if (tdb.scale != K26A_TS_TDB) k26astro_epoch_convert(&tdb, K26A_TS_TDB);
    double days = (double)tdb.days_since_J2000 + tdb.seconds_of_day / 86400.0;
    if (out_days)      *out_days      = days;
    if (out_centuries) *out_centuries = days / 36525.0;
}

/* Trig dispatch — IAU correction terms use either sin or cos. */
static double trig_(uint8_t kind, double arg_deg)
{
    double arg_rad = arg_deg * K26A_RAD_PER_DEG;
    return (kind == 0) ? sin(arg_rad) : cos(arg_rad);
}

void k26astro_rotation_eval_angles(const K26AstroIAURotation *r,
                                   const K26AstroEpoch *t,
                                   double *out_alpha_deg,
                                   double *out_delta_deg,
                                   double *out_W_deg)
{
    if (!r || !t) return;
    double T_d, T_c;
    epoch_to_T_(t, &T_d, &T_c);

    /* Linear (polynomial-of-order-1) parts. */
    double alpha = r->alpha_0 + r->alpha_rate_per_century * T_c;
    double delta = r->delta_0 + r->delta_rate_per_century * T_c;
    double W     = r->W_0     + r->W_rate_per_day        * T_d;

    /* Correction terms. arg = arg_a + arg_b * T_c for α/δ targets;
     * arg = arg_a + arg_b * T_c for W as well — the table stores
     * arg_b in units consistent with the body's argument convention
     * (most bodies use centuries; some use days, encoded in the
     * table's column header in iau_rotations.c). The evaluator
     * treats it as "argument linearly grows with T_c" uniformly;
     * day-rate arguments are pre-scaled by 36525 in the table. */
    for (int i = 0; i < r->n_corrections; i++) {
        const K26AstroIAURotCorrection *c = &r->corrections[i];
        double arg = c->arg_a + c->arg_b * T_c;
        double trig_val = trig_(c->trig_kind, arg);
        double term = c->amplitude * trig_val;
        if      (c->target_idx == 0) alpha += term;
        else if (c->target_idx == 1) delta += term;
        else if (c->target_idx == 2) W     += term;
    }

    if (out_alpha_deg) *out_alpha_deg = alpha;
    if (out_delta_deg) *out_delta_deg = delta;
    if (out_W_deg)     *out_W_deg     = W;
}

/* Build a quaternion from R3(W) · R1(90 - δ) · R3(90 + α). */
static K26Quat quat_from_angles_(double alpha_deg, double delta_deg, double W_deg)
{
    /* Convert to radians. */
    double a = alpha_deg * K26A_RAD_PER_DEG;
    double d = delta_deg * K26A_RAD_PER_DEG;
    double w = W_deg     * K26A_RAD_PER_DEG;

    /* R3 (z-axis rotation by angle φ) as a quaternion:
     *   q = (0, 0, sin(φ/2), cos(φ/2))
     * R1 (x-axis): q = (sin(φ/2), 0, 0, cos(φ/2)). */
    double half_a = 0.5 * (K26A_HALF_PI + a);
    K26Quat qA = { 0.0, 0.0, sin(half_a), cos(half_a) };

    double half_d = 0.5 * (K26A_HALF_PI - d);
    K26Quat qD = { sin(half_d), 0.0, 0.0, cos(half_d) };

    double half_w = 0.5 * w;
    K26Quat qW = { 0.0, 0.0, sin(half_w), cos(half_w) };

    /* Compose: q_total = qW * qD * qA  (applied in reverse to a vector). */
    K26Quat tmp = k26m3d_quat_mul(qD, qA);
    return k26m3d_quat_norm(k26m3d_quat_mul(qW, tmp));
}

K26Quat k26astro_rotation_quaternion(const K26AstroIAURotation *r,
                                     const K26AstroEpoch *t)
{
    double alpha = 0.0, delta = 0.0, W = 0.0;
    k26astro_rotation_eval_angles(r, t, &alpha, &delta, &W);
    return quat_from_angles_(alpha, delta, W);
}

void k26astro_rotation_matrix(const K26AstroIAURotation *r,
                              const K26AstroEpoch *t,
                              double out_R[9])
{
    if (!out_R) return;
    double alpha = 0.0, delta = 0.0, W = 0.0;
    k26astro_rotation_eval_angles(r, t, &alpha, &delta, &W);

    /* Build R = R3(W) · R1(90 - δ) · R3(90 + α) explicitly to avoid
     * any quaternion → matrix conversion noise. Angles in radians. */
    double a = alpha * K26A_RAD_PER_DEG;
    double d = delta * K26A_RAD_PER_DEG;
    double w = W     * K26A_RAD_PER_DEG;

    double cA = cos(K26A_HALF_PI + a), sA = sin(K26A_HALF_PI + a);
    double cD = cos(K26A_HALF_PI - d), sD = sin(K26A_HALF_PI - d);
    double cW = cos(w),                sW = sin(w);

    /* R3(θ) is rotation about +z:
     *   [ cθ -sθ  0 ]
     *   [ sθ  cθ  0 ]
     *   [  0   0  1 ] */
    /* R1(θ) is rotation about +x:
     *   [ 1  0    0  ]
     *   [ 0  cθ  -sθ ]
     *   [ 0  sθ   cθ ] */

    /* Compose R3(W) * R1(90-δ): */
    double M[9];
    M[0] = cW;       M[1] = -sW * cD;   M[2] =  sW * sD;
    M[3] = sW;       M[4] =  cW * cD;   M[5] = -cW * sD;
    M[6] = 0.0;      M[7] =  sD;        M[8] =  cD;

    /* Then * R3(90+α): */
    double a0 = M[0] * cA + M[1] * sA;
    double a1 = -M[0] * sA + M[1] * cA;
    double a2 = M[2];
    double b0 = M[3] * cA + M[4] * sA;
    double b1 = -M[3] * sA + M[4] * cA;
    double b2 = M[5];
    double c0 = M[6] * cA + M[7] * sA;
    double c1 = -M[6] * sA + M[7] * cA;
    double c2 = M[8];

    /* Store column-major (matching libk26m3d K26M4 minus translation). */
    out_R[0] = a0; out_R[1] = b0; out_R[2] = c0;
    out_R[3] = a1; out_R[4] = b1; out_R[5] = c1;
    out_R[6] = a2; out_R[7] = b2; out_R[8] = c2;
}
