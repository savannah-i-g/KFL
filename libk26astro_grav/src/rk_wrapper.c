/* rk_wrapper.c — RK4 + RK45 (Dormand-Prince) wrappers over libk26compute.
 *
 * The K26AstroGravState ↔ K26CVector adapter: flattens N-body state
 * into a 6N-double vector (x, y, z, vx, vy, vz × n_bodies), runs
 * libk26compute's RK driver, copies back. Used for non-conservative
 * force integration where the symplectic property of WH / Verlet
 * isn't relevant.
 *
 * The position storage uses K26AstroPos's sector-grid representation
 * during the K26AstroBody phase; the RK state vector flattens it to
 * a doubles offset relative to each body's sector origin. After the
 * RK step completes, k26astro_pos_normalise rebases each body's
 * sector + offset.
 *
 * Determinism: libk26compute's RK4/RK45 are deterministic; we now
 * patched its Makefile to add -ffp-contract=off so cross-platform
 * bit-identity holds in portable mode. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_core/pos.h"
#include "k26astro_core/epoch.h"
#include "k26compute.h"

#include <stdlib.h>
#include <string.h>

/* RHS callback: derivative of y = [x, y, z, vx, vy, vz, ...] is
 * dy/dt = [vx, vy, vz, ax, ay, az, ...] where (ax, ay, az) come from
 * a full accel_total evaluation at the current y. */
static int rhs_(double t, const K26CVector *y, K26CVector *dydt, void *user)
{
    (void)t;
    K26AstroGravState *state = (K26AstroGravState *)user;
    int n = state->n_bodies;
    if (!y || !dydt) return 1;
    if (y->n   != (size_t)(6 * n)) return 1;
    if (dydt->n != (size_t)(6 * n)) return 1;

    /* Unpack y into state->bodies (positions as offsets from sector
     * origin; velocities directly). We use a temporary K26AstroBody
     * shallow-copy so we don't disturb the user's state. */
    K26AstroBody *scratch = (K26AstroBody *)malloc(sizeof(K26AstroBody) * (size_t)n);
    if (!scratch) return 2;
    memcpy(scratch, state->bodies, sizeof(K26AstroBody) * (size_t)n);

    for (int i = 0; i < n; i++) {
        double xi = y->data[6*i + 0];
        double yi = y->data[6*i + 1];
        double zi = y->data[6*i + 2];
        scratch[i].pos = k26astro_pos_from_m(xi, yi, zi);
        scratch[i].vel.x = y->data[6*i + 3];
        scratch[i].vel.y = y->data[6*i + 4];
        scratch[i].vel.z = y->data[6*i + 5];
    }

    K26AstroGravState shadow = *state;
    shadow.bodies = scratch;

    K26V3 *accel = (K26V3 *)calloc((size_t)n, sizeof(K26V3));
    if (!accel) { free(scratch); return 3; }
    k26astro_grav_accel_total(&shadow, accel);

    for (int i = 0; i < n; i++) {
        dydt->data[6*i + 0] = scratch[i].vel.x;
        dydt->data[6*i + 1] = scratch[i].vel.y;
        dydt->data[6*i + 2] = scratch[i].vel.z;
        dydt->data[6*i + 3] = accel[i].x;
        dydt->data[6*i + 4] = accel[i].y;
        dydt->data[6*i + 5] = accel[i].z;
    }

    free(accel);
    free(scratch);
    return 0;
}

/* Pack/unpack helpers between K26AstroBody[] and K26CVector. */
static void pack_(const K26AstroBody *bodies, int n, K26CVector *y)
{
    for (int i = 0; i < n; i++) {
        K26V3 r = k26astro_pos_to_m_approx(&bodies[i].pos);
        y->data[6*i + 0] = r.x;
        y->data[6*i + 1] = r.y;
        y->data[6*i + 2] = r.z;
        y->data[6*i + 3] = bodies[i].vel.x;
        y->data[6*i + 4] = bodies[i].vel.y;
        y->data[6*i + 5] = bodies[i].vel.z;
    }
}

static void unpack_(K26AstroBody *bodies, int n, const K26CVector *y)
{
    for (int i = 0; i < n; i++) {
        bodies[i].pos = k26astro_pos_from_m(
            y->data[6*i + 0],
            y->data[6*i + 1],
            y->data[6*i + 2]);
        bodies[i].vel.x = y->data[6*i + 3];
        bodies[i].vel.y = y->data[6*i + 4];
        bodies[i].vel.z = y->data[6*i + 5];
    }
}

int k26astro_grav_step_rk(K26AstroGravState *state, double dt)
{
    if (!state) return K26ASTRO_E_NULL;
    int n = state->n_bodies;
    size_t dim = (size_t)(6 * n);

    K26CVector y;
    y.n = dim;
    y.data = (double *)calloc(dim, sizeof(double));
    if (!y.data) return K26ASTRO_E_ALLOC;

    pack_(state->bodies, n, &y);

    K26CStatus rc;
    if (state->integrator == K26ASTRO_INTEGRATOR_RK4) {
        rc = k26c_ode_rk4(rhs_, state, 0.0, dt, /*n_steps*/ 1, &y);
    } else {
        /* RK45 (Dormand-Prince) with reasonable default tolerances. */
        rc = k26c_ode_rk45(rhs_, state, 0.0, dt, /*rtol*/ 1e-9,
                            /*atol*/ 1e-12, &y);
    }
    if (rc != K26C_OK) { free(y.data); return K26ASTRO_E_NO_CONVERGE; }

    unpack_(state->bodies, n, &y);
    free(y.data);

    k26astro_epoch_add_seconds(&state->t, dt);
    state->dt_last = dt;
    return K26ASTRO_E_OK;
}
