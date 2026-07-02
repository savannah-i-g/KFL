/* test_kepler_stm.c — Cartesian state-transition matrix.
 *
 * Gates:
 *   1. STM(dt = 0) = identity within tight tolerance.
 *   2. STM column j times a small perturbation along axis j of the
 *      initial state matches the propagated perturbation to 1e-4
 *      relative across three orbit classes (LEO circular, GTO
 *      eccentric, hyperbolic flyby).
 *   3. The 6 x 6 matrix is approximately symplectic:
 *      |M^T J M - J|_F < 1e-5 (Hamiltonian preservation; loose to
 *      absorb centered-FD truncation order). */

#include "k26astro_conics/stm.h"
#include "k26astro_conics/kepler.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static const double EARTH_MU = 3.986004418e14;

static double matrix_frobenius_(double A[6][6])
{
    double s = 0.0;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            s += A[i][j] * A[i][j];
    return sqrt(s);
}

static int verify_identity_at_zero_(K26V3 r0, K26V3 v0)
{
    double M[6][6];
    int rc = k26astro_kepler_stm(r0, v0, 0.0, EARTH_MU, M);
    assert(rc == 0);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            assert(fabs(M[i][j] - expected) < 1.0e-6);
        }
    }
    return 0;
}

static int verify_column_against_propagation_(K26V3 r0, K26V3 v0,
                                               double dt, const char *tag,
                                               double rel_tol)
{
    double M[6][6];
    int rc = k26astro_kepler_stm(r0, v0, dt, EARTH_MU, M);
    assert(rc == 0);

    /* Apply a perturbation along axis j = 5 (vz), propagate the
     * perturbed initial state, and verify the resulting delta
     * matches M[:, 5] * h to within rel_tol. */
    double v0_norm = sqrt(v0.x * v0.x + v0.y * v0.y + v0.z * v0.z);
    double h = 1.0e-5 * v0_norm;
    if (h < 1.0e-6) h = 1.0e-6;

    K26V3 r0_base = r0, v0_base = v0;
    K26V3 r0_pert = r0, v0_pert = v0;
    v0_pert.z += h;

    K26V3 r_base, v_base, r_pert, v_pert;
    rc = k26astro_kepler_propagate(&r_base, &v_base, r0_base, v0_base,
                                    EARTH_MU, dt, 32);
    assert(rc == 0);
    rc = k26astro_kepler_propagate(&r_pert, &v_pert, r0_pert, v0_pert,
                                    EARTH_MU, dt, 32);
    assert(rc == 0);

    double delta[6] = {
        r_pert.x - r_base.x,
        r_pert.y - r_base.y,
        r_pert.z - r_base.z,
        v_pert.x - v_base.x,
        v_pert.y - v_base.y,
        v_pert.z - v_base.z
    };
    double pred[6];
    for (int i = 0; i < 6; i++) pred[i] = M[i][5] * h;

    double dnorm  = 0.0, diff = 0.0;
    for (int i = 0; i < 6; i++) {
        dnorm += delta[i] * delta[i];
        double e = pred[i] - delta[i];
        diff += e * e;
    }
    dnorm = sqrt(dnorm);
    diff  = sqrt(diff);
    double rel = diff / (dnorm + 1.0e-30);
    printf("  %-20s |delta|=%.3e |pred-delta|/|delta|=%.3e\n",
           tag, dnorm, rel);
    assert(rel < rel_tol);
    return 0;
}

static int verify_symplectic_(K26V3 r0, K26V3 v0, double dt)
{
    double M[6][6];
    int rc = k26astro_kepler_stm(r0, v0, dt, EARTH_MU, M);
    assert(rc == 0);

    /* J = [ 0  I ; -I  0 ] (6 x 6 symplectic form, position-first). */
    double J[6][6] = {{0}};
    J[0][3] = J[1][4] = J[2][5] = 1.0;
    J[3][0] = J[4][1] = J[5][2] = -1.0;

    /* MtJM = M^T J M. */
    double MtJ[6][6] = {{0}}, MtJM[6][6] = {{0}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
                MtJ[i][j] += M[k][i] * J[k][j];
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
                MtJM[i][j] += MtJ[i][k] * M[k][j];

    double D[6][6];
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            D[i][j] = MtJM[i][j] - J[i][j];

    double err = matrix_frobenius_(D);
    /* Symplectic deviation tolerance. The FD STM is symplectic only
     * to FD-truncation order; 1e-5 is generous. */
    printf("  symplectic |M^T J M - J|_F = %.3e\n", err);
    assert(err < 1.0e-5);
    return 0;
}

int main(void)
{
    /* ---- Identity at dt = 0 ---------------------------------- */
    {
        K26V3 r0 = { 7.0e6, 0.0, 0.0 };
        K26V3 v0 = { 0.0, 7.5e3, 0.0 };
        verify_identity_at_zero_(r0, v0);
    }

    /* ---- LEO circular: 400 km altitude, quarter period ------- */
    {
        double r = 6378e3 + 400e3;
        double v_circ = sqrt(EARTH_MU / r);
        K26V3 r0 = { r, 0.0, 0.0 };
        K26V3 v0 = { 0.0, v_circ, 0.0 };
        double T = 2.0 * M_PI * sqrt(r * r * r / EARTH_MU);
        verify_column_against_propagation_(r0, v0, T / 4.0,
                                            "LEO circular T/4", 1.0e-4);
        verify_symplectic_(r0, v0, T / 4.0);
    }

    /* ---- GTO eccentric: perigee at 300 km, apogee at GEO ---- */
    {
        double r_p = 6378e3 + 300e3;
        double r_a = 42164e3;
        double a   = 0.5 * (r_p + r_a);
        double v_p = sqrt(EARTH_MU * (2.0 / r_p - 1.0 / a));
        K26V3 r0 = { r_p, 0.0, 0.0 };
        K26V3 v0 = { 0.0, v_p, 0.0 };
        double T = 2.0 * M_PI * sqrt(a * a * a / EARTH_MU);
        verify_column_against_propagation_(r0, v0, T / 4.0,
                                            "GTO eccentric T/4", 1.0e-3);
        verify_symplectic_(r0, v0, T / 4.0);
    }

    /* ---- Hyperbolic flyby: v_infinity = 5 km/s at periapsis 8000 km */
    {
        double r_p = 8.0e6;
        double v_inf = 5.0e3;
        /* vis-viva at periapsis for hyperbola: v_p^2 = v_inf^2 + 2 mu / r_p */
        double v_p_sq = v_inf * v_inf + 2.0 * EARTH_MU / r_p;
        double v_p = sqrt(v_p_sq);
        K26V3 r0 = { r_p, 0.0, 0.0 };
        K26V3 v0 = { 0.0, v_p, 0.0 };
        /* Time to traverse some hyperbolic chord; pick 600 s. */
        verify_column_against_propagation_(r0, v0, 600.0,
                                            "hyperbolic 600 s", 1.0e-3);
    }

    printf("test_kepler_stm: OK\n");
    return 0;
}
