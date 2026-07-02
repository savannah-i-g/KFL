/* libk26compute — ODE solvers.
 *
 * RK4: classical fixed-step 4th-order Runge-Kutta. Cheap, robust,
 *      adequate when step size is known to be safe and the user just
 *      wants a result.
 *
 * RK45: Dormand-Prince 5(4) — embedded pair, standard adaptive
 *       step-size controller using local error estimate. The
 *       coefficients are the canonical DP table; do not modify
 *       without re-running ODE test problems (e.g., harmonic
 *       oscillator stays in phase, exponential decay matches
 *       analytical y(t) within tolerance). */

#include "k26compute.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Helper: compute dydt = rhs(t, y) and bounce errors back. */
static int call_rhs(K26CRhsFn rhs, void *user, double t,
                    const K26CVector *y, K26CVector *dydt)
{
    return rhs(t, y, dydt, user);
}

K26CStatus k26c_ode_rk4(K26CRhsFn rhs, void *user,
                        double t0, double t1, size_t n_steps,
                        K26CVector *y)
{
    if (!rhs || !y || !y->data || n_steps == 0) return K26C_ERR_INVAL;
    size_t dim = y->n;
    double h   = (t1 - t0) / (double)n_steps;

    K26CVector k1, k2, k3, k4, ytmp;
    if (k26c_vec_alloc(&k1,   dim) != K26C_OK) return K26C_ERR_OOM;
    if (k26c_vec_alloc(&k2,   dim) != K26C_OK) { k26c_vec_free(&k1); return K26C_ERR_OOM; }
    if (k26c_vec_alloc(&k3,   dim) != K26C_OK) { k26c_vec_free(&k1); k26c_vec_free(&k2); return K26C_ERR_OOM; }
    if (k26c_vec_alloc(&k4,   dim) != K26C_OK) { k26c_vec_free(&k1); k26c_vec_free(&k2); k26c_vec_free(&k3); return K26C_ERR_OOM; }
    if (k26c_vec_alloc(&ytmp, dim) != K26C_OK) { k26c_vec_free(&k1); k26c_vec_free(&k2); k26c_vec_free(&k3); k26c_vec_free(&k4); return K26C_ERR_OOM; }

    double t = t0;
    K26CStatus rc = K26C_OK;
    for (size_t step = 0; step < n_steps; step++) {
        if (call_rhs(rhs, user, t, y, &k1) != 0) { rc = K26C_ERR_CONVERGE; break; }
        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + 0.5 * h * k1.data[i];
        if (call_rhs(rhs, user, t + 0.5 * h, &ytmp, &k2) != 0) { rc = K26C_ERR_CONVERGE; break; }
        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + 0.5 * h * k2.data[i];
        if (call_rhs(rhs, user, t + 0.5 * h, &ytmp, &k3) != 0) { rc = K26C_ERR_CONVERGE; break; }
        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] +       h * k3.data[i];
        if (call_rhs(rhs, user, t +       h, &ytmp, &k4) != 0) { rc = K26C_ERR_CONVERGE; break; }
        for (size_t i = 0; i < dim; i++)
            y->data[i] += (h / 6.0) * (k1.data[i] + 2.0*k2.data[i] + 2.0*k3.data[i] + k4.data[i]);
        t += h;
    }

    k26c_vec_free(&k1); k26c_vec_free(&k2); k26c_vec_free(&k3);
    k26c_vec_free(&k4); k26c_vec_free(&ytmp);
    return rc;
}

/* Dormand-Prince 5(4) coefficients (Butcher tableau). */
static const double DP_C2 = 1.0/5.0;
static const double DP_C3 = 3.0/10.0;
static const double DP_C4 = 4.0/5.0;
static const double DP_C5 = 8.0/9.0;
static const double DP_C6 = 1.0;
static const double DP_C7 = 1.0;

static const double DP_A21 = 1.0/5.0;
static const double DP_A31 = 3.0/40.0,        DP_A32 = 9.0/40.0;
static const double DP_A41 = 44.0/45.0,       DP_A42 = -56.0/15.0,    DP_A43 = 32.0/9.0;
static const double DP_A51 = 19372.0/6561.0,  DP_A52 = -25360.0/2187.0, DP_A53 = 64448.0/6561.0,  DP_A54 = -212.0/729.0;
static const double DP_A61 = 9017.0/3168.0,   DP_A62 = -355.0/33.0,   DP_A63 = 46732.0/5247.0,    DP_A64 = 49.0/176.0,    DP_A65 = -5103.0/18656.0;
/* Row 7 of the A table (DP_A71..DP_A76) is omitted: the FSAL
 * property of Dormand-Prince makes those weights equal to the
 * 5th-order b_i, so y_new IS the row-7 combination of k_1..k_6 by
 * construction. The 7th-stage RHS evaluation reuses y_new directly,
 * which is why this implementation never references a separate
 * row-7 sum. */

/* 5th-order solution weights (b). */
static const double DP_B1 = 35.0/384.0,    DP_B3 = 500.0/1113.0,   DP_B4 = 125.0/192.0,
                    DP_B5 = -2187.0/6784.0, DP_B6 = 11.0/84.0;

/* 4th-order solution weights (b*); b - b* gives the error estimator. */
static const double DP_BS1 = 5179.0/57600.0,  DP_BS3 = 7571.0/16695.0, DP_BS4 = 393.0/640.0,
                    DP_BS5 = -92097.0/339200.0, DP_BS6 = 187.0/2100.0,  DP_BS7 = 1.0/40.0;

K26CStatus k26c_ode_rk45(K26CRhsFn rhs, void *user,
                         double t0, double t1,
                         double rtol, double atol,
                         K26CVector *y)
{
    if (!rhs || !y || !y->data) return K26C_ERR_INVAL;
    if (rtol <= 0.0) rtol = 1e-6;
    if (atol <  0.0) atol = 1e-9;
    if (t1 == t0)    return K26C_OK;

    size_t dim = y->n;
    int direction = (t1 > t0) ? 1 : -1;
    double h = (t1 - t0) / 100.0;            /* initial step guess */
    if (h == 0.0) h = direction * 1e-3;

    /* Per-stage k vectors + ytmp + y_new + err. */
    K26CVector k1, k2, k3, k4, k5, k6, k7, ytmp, y_new, y_err;
    K26CVector *vs[] = { &k1, &k2, &k3, &k4, &k5, &k6, &k7, &ytmp, &y_new, &y_err };
    for (size_t i = 0; i < sizeof(vs)/sizeof(vs[0]); i++) {
        if (k26c_vec_alloc(vs[i], dim) != K26C_OK) {
            for (size_t j = 0; j < i; j++) k26c_vec_free(vs[j]);
            return K26C_ERR_OOM;
        }
    }

    double t = t0;
    size_t max_steps = 100000;     /* belt-and-braces ceiling */
    K26CStatus rc = K26C_OK;

    for (size_t step = 0; step < max_steps; step++) {
        if ((direction > 0 && t >= t1) ||
            (direction < 0 && t <= t1)) break;

        /* Don't overshoot. */
        if ((direction > 0 && t + h > t1) ||
            (direction < 0 && t + h < t1)) {
            h = t1 - t;
        }

        if (call_rhs(rhs, user, t, y, &k1) != 0) { rc = K26C_ERR_CONVERGE; break; }

        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + h*(DP_A21*k1.data[i]);
        if (call_rhs(rhs, user, t + DP_C2*h, &ytmp, &k2) != 0) { rc = K26C_ERR_CONVERGE; break; }

        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + h*(DP_A31*k1.data[i] + DP_A32*k2.data[i]);
        if (call_rhs(rhs, user, t + DP_C3*h, &ytmp, &k3) != 0) { rc = K26C_ERR_CONVERGE; break; }

        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + h*(DP_A41*k1.data[i] + DP_A42*k2.data[i] + DP_A43*k3.data[i]);
        if (call_rhs(rhs, user, t + DP_C4*h, &ytmp, &k4) != 0) { rc = K26C_ERR_CONVERGE; break; }

        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + h*(DP_A51*k1.data[i] + DP_A52*k2.data[i] + DP_A53*k3.data[i] + DP_A54*k4.data[i]);
        if (call_rhs(rhs, user, t + DP_C5*h, &ytmp, &k5) != 0) { rc = K26C_ERR_CONVERGE; break; }

        for (size_t i = 0; i < dim; i++) ytmp.data[i] = y->data[i] + h*(DP_A61*k1.data[i] + DP_A62*k2.data[i] + DP_A63*k3.data[i] + DP_A64*k4.data[i] + DP_A65*k5.data[i]);
        if (call_rhs(rhs, user, t + DP_C6*h, &ytmp, &k6) != 0) { rc = K26C_ERR_CONVERGE; break; }

        /* 5th-order solution. */
        for (size_t i = 0; i < dim; i++)
            y_new.data[i] = y->data[i] + h*(DP_B1*k1.data[i] + DP_B3*k3.data[i] + DP_B4*k4.data[i] + DP_B5*k5.data[i] + DP_B6*k6.data[i]);

        /* k7 used only for the error estimator (FSAL stage). */
        if (call_rhs(rhs, user, t + DP_C7*h, &y_new, &k7) != 0) { rc = K26C_ERR_CONVERGE; break; }

        /* Error estimate = h * Σ (b_i - b*_i) k_i. b2/b7 are zero in
         * the 5th-order weights; b2/b3/b4/b5/b6/b7 in the 4th-order
         * weights.  Element-wise. */
        for (size_t i = 0; i < dim; i++) {
            double err =
                  (DP_B1  - DP_BS1) * k1.data[i]
                + (DP_B3  - DP_BS3) * k3.data[i]
                + (DP_B4  - DP_BS4) * k4.data[i]
                + (DP_B5  - DP_BS5) * k5.data[i]
                + (DP_B6  - DP_BS6) * k6.data[i]
                + (0.0    - DP_BS7) * k7.data[i];
            y_err.data[i] = h * err;
        }

        /* Norm-of-error / scale (per-component tolerance). */
        double err_norm = 0.0;
        for (size_t i = 0; i < dim; i++) {
            double scale = atol + rtol * fmax(fabs(y->data[i]), fabs(y_new.data[i]));
            if (scale == 0.0) scale = 1.0;
            double r = y_err.data[i] / scale;
            err_norm += r * r;
        }
        err_norm = sqrt(err_norm / (double)dim);

        if (err_norm <= 1.0) {
            /* Accept step. */
            t += h;
            memcpy(y->data, y_new.data, dim * sizeof(double));
        }
        /* Adjust step size. A standard PI controller would be better
         * but this elementary I-controller (h_new = 0.9 * h * err^(-1/5))
         * is adequate here and avoids tuning a second history. */
        double factor;
        if (err_norm == 0.0)        factor = 5.0;
        else                         factor = 0.9 * pow(err_norm, -0.2);
        if (factor > 5.0)  factor = 5.0;
        if (factor < 0.1)  factor = 0.1;
        h *= factor;
    }

    for (size_t i = 0; i < sizeof(vs)/sizeof(vs[0]); i++) k26c_vec_free(vs[i]);
    return rc;
}
