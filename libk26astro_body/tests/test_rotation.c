/* test_rotation.c — IAU rotation model evaluator + table coverage. */
#include "k26astro_body/rotation_model.h"
#include "k26astro_core/consts.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* ---- Table coverage ----------------------------------- */
    const K26AstroIAURotation *earth = k26astro_rotation_lookup("iau2018:earth");
    assert(earth != NULL);
    assert(earth->naif_id == 399);

    const K26AstroIAURotation *mars  = k26astro_rotation_lookup("iau2018:mars");
    assert(mars != NULL);
    assert(mars->n_corrections > 5);

    const K26AstroIAURotation *moon  = k26astro_rotation_lookup("iau2018:moon");
    assert(moon != NULL);
    assert(moon->n_corrections > 10);

    const K26AstroIAURotation *titan = k26astro_rotation_lookup("iau2018:titan");
    assert(titan != NULL);

    /* Case-insensitive should work. */
    assert(k26astro_rotation_lookup("IAU2018:EARTH") != NULL);

    /* NAIF-id lookup. */
    const K26AstroIAURotation *mars_by_id = k26astro_rotation_by_naif(499);
    assert(mars_by_id == mars);

    /* Unknown body. */
    assert(k26astro_rotation_lookup("iau2018:fictional") == NULL);

    /* ---- Evaluation at J2000 ----------------------------- */
    K26AstroEpoch j2000 = k26astro_epoch_j2000_tt();
    double alpha = 0.0, delta = 0.0, W = 0.0;
    k26astro_rotation_eval_angles(earth, &j2000, &alpha, &delta, &W);
    /* At J2000, Earth's α₀ = 0 (within the table's polynomial). */
    assert(fabs(alpha) < 1.0e-3);
    /* δ₀ = 90 — Earth's pole points to +z in ICRF. */
    assert(fabs(delta - 90.0) < 1.0e-3);

    /* W advances 360.9856 deg/day; at J2000 itself W = 190.147. */
    assert(fabs(W - 190.147) < 1.0e-3);

    /* ---- W advances correctly over time ----------------- */
    K26AstroEpoch j2001 = j2000;
    k26astro_epoch_add_seconds(&j2001, 86400.0);
    double W1 = 0.0;
    k26astro_rotation_eval_angles(earth, &j2001, NULL, NULL, &W1);
    /* W1 - W0 ≈ 360.9856 deg (mod 360). */
    double dW = fmod(W1 - W, 360.0);
    if (dW > 180.0)  dW -= 360.0;
    if (dW < -180.0) dW += 360.0;
    assert(fabs(dW - 0.9856) < 1e-3);

    /* ---- Quaternion build ------------------------------- */
    K26Quat q = k26astro_rotation_quaternion(earth, &j2000);
    /* Should be unit. */
    double n = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    assert(fabs(n - 1.0) < 1e-10);

    /* ---- Rotation matrix is orthonormal ---------------- */
    double R[9];
    k26astro_rotation_matrix(earth, &j2000, R);
    /* Check orthogonality: R^T * R = I (within floating-point). */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) {
                s += R[k * 3 + i] * R[k * 3 + j];
            }
            double expected = (i == j) ? 1.0 : 0.0;
            assert(fabs(s - expected) < 1e-10);
        }
    }

    /* ---- All registered models build a valid quaternion --- */
    int n_models = 0;
    for (const K26AstroIAURotation *r = k26astro_rotation_table();
         r->name; r++) {
        n_models++;
        K26Quat qm = k26astro_rotation_quaternion(r, &j2000);
        double nm = qm.x * qm.x + qm.y * qm.y + qm.z * qm.z + qm.w * qm.w;
        assert(fabs(nm - 1.0) < 1e-10);
    }
    assert(n_models > 20);   /* At least 20 named bodies. */

    printf("test_rotation: OK (%d models in table)\n", n_models);
    return 0;
}
