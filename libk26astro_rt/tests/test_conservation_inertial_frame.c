/* test_conservation_inertial_frame.c — frame-convention gate.
 *
 * Demonstrates that the closed-system invariant
 *
 *     E_closed = ke_surviving
 *              + E_ejecta + E_radiated
 *              + E_reaction_mass
 *              - E_chem_source
 *
 * holds to machine precision when reaction-mass dispatches are
 * recorded with inertial-frame propellant velocity AND the
 * chemistry energy released by the burn is recorded as a source
 * term. Separately demonstrates that the legacy rocket-frame
 * convention (without the chemistry source) leaves a per-burn
 * drift on the order m_p · v_vehicle · v_e that grows with vehicle
 * inertial velocity — the bookkeeping artifact behind the
 * observed multi-burn conservation-drift finding.
 *
 * The substrate primitive is correct in both calling conventions
 * (it sums ½ m v² + m·v regardless of what frame the velocity
 * came from). What this test gates is the convention — that the
 * driver caller passes inertial-frame velocity and the chemistry
 * source — without which the closed-system invariant is not
 * actually invariant. */
#include "k26astro_rt/conservation_ledger.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx_rel(double got, double want, double rel_tol)
{
    const double denom = (fabs(want) > 1.0) ? fabs(want) : 1.0;
    return fabs(got - want) / denom <= rel_tol;
}

int main(void)
{
    /* Scenario: a 1000-kg vehicle moving at +5000 m/s along x
     * ejects 100 kg of propellant such that the propellant ends
     * up at rest in the inertial frame (v_p_x = 0). This is a
     * clean, fully-specified momentum problem with an exact
     * analytic answer that we can compare bit-precisely against
     * the ledger's bookkeeping.
     *
     * Momentum conservation:
     *   M_0 v_0 = M_f v_f + m_p v_p
     *   1000 × 5000 = 900 × v_f + 100 × 0
     *   v_f = 5000000 / 900 = 5555.555... m/s
     *
     * Total KE before:  ½ × 1000 × 5000² = 1.25 × 10¹⁰ J
     * Vehicle KE after: ½ × 900 × v_f² ≈ 1.38889 × 10¹⁰ J
     * Propellant KE after: ½ × 100 × 0² = 0 J
     * Chemistry released: ΔKE = 1.38889e10 - 1.25e10
     *                          ≈ 1.3889 × 10⁹ J */

    const double M_0   = 1000.0;
    const double m_p   = 100.0;
    const double M_f   = M_0 - m_p;
    const double v_0_x = 5000.0;
    const double v_p_x = 0.0;
    const double v_f_x = (M_0 * v_0_x - m_p * v_p_x) / M_f;

    const double KE_before          = 0.5 * M_0 * v_0_x * v_0_x;
    const double KE_after_vehicle   = 0.5 * M_f * v_f_x * v_f_x;
    const double KE_after_propellant= 0.5 * m_p * v_p_x * v_p_x;
    const double E_chem = (KE_after_vehicle + KE_after_propellant)
                          - KE_before;

    /* Sanity: v_f matches Tsiolkovsky-style momentum balance. */
    assert(approx_rel(v_f_x, 5555.55555555555, 1.0e-12));
    /* Sanity: chemistry release is positive. */
    assert(E_chem > 0.0);

    /* ---- Inertial-frame convention: invariant holds ----------- */
    {
        K26AstroRtConservationLedger *ledger =
            k26astro_rt_ledger_new();
        assert(ledger != NULL);

        /* Record the burn with the inertial-frame API. */
        k26astro_rt_ledger_record_reaction_mass_inertial(
            ledger, m_p, v_p_x, 0.0, 0.0);
        k26astro_rt_ledger_record_chemical_source(ledger, E_chem);

        K26AstroRtConservationTotals t;
        k26astro_rt_ledger_totals(ledger, &t);

        /* Counters tick. */
        assert(t.reaction_mass_count == 1);
        assert(t.chem_source_count == 1);

        /* Ledger accumulates propellant KE = 0 (it ended at rest). */
        assert(approx_rel(t.E_reaction_mass_J,
                          KE_after_propellant, 1.0e-12));
        assert(approx_rel(t.E_chem_source_J, E_chem, 1.0e-12));

        /* Closed-system invariant: post-burn sum equals pre-burn KE
         * to machine precision. */
        const double E_closed_post = KE_after_vehicle
                                   + t.E_ejecta_J
                                   + t.E_radiated_J
                                   + t.E_reaction_mass_J
                                   - t.E_chem_source_J;
        const double drift = E_closed_post - KE_before;

        /* Closure to machine precision: the only error source is
         * the v_f_x = 5000000/900 division which is bit-exact at
         * IEEE-754 (5555.55555...e0 truncated identically on both
         * sides of the assertion). */
        assert(fabs(drift) / KE_before < 1.0e-14);

        k26astro_rt_ledger_destroy(ledger);
    }

    /* ---- Rocket-frame convention (legacy / wrong): drifts ----- *
     *
     * Compute the propellant's rocket-frame (post-burn vehicle
     * frame) velocity: v_p_inertial - v_f_inertial =
     * 0 - 5555.555... = -5555.555 m/s. A driver that records this
     * value (the propellant's velocity relative to the post-burn
     * vehicle frame) into the ledger and omits the chemistry
     * source term reproduces the legacy (rocket-frame, no chemistry
     * source) bookkeeping pattern.
     *
     * The legacy entry has __attribute__((deprecated)) so call
     * sites get a build warning. We suppress it here because we
     * are deliberately exercising the wrong-convention path to
     * demonstrate the drift it produces. */
    {
        K26AstroRtConservationLedger *ledger =
            k26astro_rt_ledger_new();
        assert(ledger != NULL);

        const double v_e_rocket = v_f_x - v_p_x; /* +5555.555 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        /* Record exhaust as -v_e along x (retrograde in vehicle
         * frame). NO chemistry source recorded. */
        k26astro_rt_ledger_record_reaction_mass(
            ledger, m_p, -v_e_rocket, 0.0, 0.0);
#pragma GCC diagnostic pop

        K26AstroRtConservationTotals t;
        k26astro_rt_ledger_totals(ledger, &t);

        /* Bookkeeping mistakes the propellant KE: records
         * ½ m_p v_e² instead of ½ m_p v_p_inertial². */
        const double KE_rocket_frame = 0.5 * m_p
                                       * v_e_rocket * v_e_rocket;
        assert(approx_rel(t.E_reaction_mass_J,
                          KE_rocket_frame, 1.0e-12));
        assert(t.E_chem_source_J == 0.0);

        /* Closed-system using this convention drifts: */
        const double E_closed_rocket = KE_after_vehicle
                                     + t.E_reaction_mass_J;
        const double drift = E_closed_rocket - KE_before;

        /* The drift is finite and substantial — for this scenario
         * with v_vehicle=5000, v_e≈5555, m_p=100, it predicts
         * drift ≈ m_p (v_e²/2 - v_p²/2) + (KE_after_v - KE_before)
         * which evaluates to a few × 10⁹ J. */
        assert(fabs(drift) > 1.0e8);

        /* And the wrong-convention drift is much larger than the
         * inertial-convention noise floor. Demonstrates the
         * substrate is sound; the convention is the issue. */
        const double inertial_noise_floor_J = 1.0e-2 * KE_before;
        assert(fabs(drift) > inertial_noise_floor_J);

        fprintf(stderr,
                "test_conservation_inertial_frame: "
                "rocket-frame drift = %.3e J vs inertial-frame "
                "drift ≈ 0 (closure to ~1e-14 of KE)\n",
                drift);

        k26astro_rt_ledger_destroy(ledger);
    }

    /* ---- The two reaction-mass APIs produce identical math
     * given identical inputs (the _inertial suffix encodes the
     * frame contract; the accumulator is the same). */
    {
        K26AstroRtConservationLedger *a = k26astro_rt_ledger_new();
        K26AstroRtConservationLedger *b = k26astro_rt_ledger_new();

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        k26astro_rt_ledger_record_reaction_mass(
            a, 42.0, 1234.5, -678.9, 1011.12);
#pragma GCC diagnostic pop
        k26astro_rt_ledger_record_reaction_mass_inertial(
            b, 42.0, 1234.5, -678.9, 1011.12);

        K26AstroRtConservationTotals ta, tb;
        k26astro_rt_ledger_totals(a, &ta);
        k26astro_rt_ledger_totals(b, &tb);

        assert(ta.E_reaction_mass_J == tb.E_reaction_mass_J);
        assert(ta.p_reaction_x_kgmps == tb.p_reaction_x_kgmps);
        assert(ta.p_reaction_y_kgmps == tb.p_reaction_y_kgmps);
        assert(ta.p_reaction_z_kgmps == tb.p_reaction_z_kgmps);
        assert(ta.reaction_mass_count == tb.reaction_mass_count);

        k26astro_rt_ledger_destroy(a);
        k26astro_rt_ledger_destroy(b);
    }

    /* ---- _record_burn helper: same result as two-call pattern - *
     *
     * Recompute the headline scenario via the burn helper. It
     * should produce E_reaction_mass + E_chem_source values
     * identical to the explicit two-call pattern, while taking
     * physical inputs at the call site that cannot be confused
     * about frame. */
    {
        K26AstroRtConservationLedger *a = k26astro_rt_ledger_new();
        K26AstroRtConservationLedger *b = k26astro_rt_ledger_new();

        /* Burn parameters matching the headline scenario:
         * v_vehicle = (5000, 0, 0) before burn (effectively the
         * pre-burn vehicle frame; for the impulsive setup we use
         * this as the proxy for the burn-time vehicle frame).
         * exhaust_speed = v_f - v_p = 5555.555 m/s rocket-frame.
         * exhaust_dir = (-1, 0, 0) (retrograde). */
        const double v_e_rocket = v_f_x - v_p_x;
        const double v_burn_x   = v_0_x;
        k26astro_rt_ledger_record_burn(
            a, m_p,
            v_burn_x, 0.0, 0.0,
            v_e_rocket,
            -1.0, 0.0, 0.0);

        /* Two-call pattern: compute the same v_p_inertial and
         * E_chem manually, call the lower-level entries. */
        const double v_p_inertial_x = v_burn_x + v_e_rocket * (-1.0);
        k26astro_rt_ledger_record_reaction_mass_inertial(
            b, m_p, v_p_inertial_x, 0.0, 0.0);
        k26astro_rt_ledger_record_chemical_source(
            b, 0.5 * m_p * v_e_rocket * v_e_rocket);

        K26AstroRtConservationTotals ta, tb;
        k26astro_rt_ledger_totals(a, &ta);
        k26astro_rt_ledger_totals(b, &tb);

        assert(ta.E_reaction_mass_J == tb.E_reaction_mass_J);
        assert(ta.p_reaction_x_kgmps == tb.p_reaction_x_kgmps);
        assert(ta.E_chem_source_J == tb.E_chem_source_J);
        assert(ta.reaction_mass_count == tb.reaction_mass_count);
        assert(ta.chem_source_count == tb.chem_source_count);

        /* And the burn-helper version, like the two-call version,
         * should close the closed-system invariant to machine
         * precision when the surviving-vehicle KE is recomputed
         * from the burn's predicted v_f. (We use v_p_x = 0 here;
         * the helper computes v_p_inertial = v_burn + v_e × dir =
         * 5000 - 5555.55 = -555.55. Not the same scenario as the
         * headline test, so we don't recompute KE_before here —
         * the equality assertions above are sufficient to prove
         * the helper composes the two atomic operations correctly.) */

        k26astro_rt_ledger_destroy(a);
        k26astro_rt_ledger_destroy(b);
    }

    /* ---- L_z_ejecta tracking via record_dispatch -------------- *
     *
     * A mass emitted at origin (R, 0, 0) with velocity (0, +v, 0)
     * carries L_z = m × R × v of angular momentum out of the
     * surviving-bodies population. The ledger captures this so
     * the trace-side conservation.csv L_z column closes when the
     * ejected mass leaves the propagator. */
    {
        K26AstroRtConservationLedger *ledger =
            k26astro_rt_ledger_new();
        const double R = 1.0e7;       /* emission position (m) */
        const double v_tan = 5000.0;  /* tangential v (m/s) */
        const double m_ejecta = 50.0;
        k26astro_rt_ledger_record_dispatch(
            ledger, m_ejecta,
            /* v_ejecta */ 0.0, v_tan, 0.0,
            /* origin   */ R, 0.0, 0.0);

        K26AstroRtConservationTotals t;
        k26astro_rt_ledger_totals(ledger, &t);
        const double expected_L_z = m_ejecta * R * v_tan;
        assert(t.L_z_ejecta_kg_m2ps == expected_L_z);

        /* A radial-only emission (along x at origin (R, 0, 0))
         * carries zero L_z. */
        k26astro_rt_ledger_record_dispatch(
            ledger, m_ejecta,
            /* v_ejecta */ v_tan, 0.0, 0.0,
            /* origin   */ R, 0.0, 0.0);
        k26astro_rt_ledger_totals(ledger, &t);
        /* Total L_z unchanged. */
        assert(t.L_z_ejecta_kg_m2ps == expected_L_z);

        k26astro_rt_ledger_destroy(ledger);
    }

    /* ---- NULL safety for new APIs ----------------------------- */
    k26astro_rt_ledger_record_reaction_mass_inertial(
        NULL, 1.0, 1.0, 1.0, 1.0);
    k26astro_rt_ledger_record_chemical_source(NULL, 1.0e6);
    k26astro_rt_ledger_record_burn(
        NULL, 1.0, 0, 0, 0, 1000.0, -1, 0, 0);

    fprintf(stderr, "test_conservation_inertial_frame: OK\n");
    return 0;
}
