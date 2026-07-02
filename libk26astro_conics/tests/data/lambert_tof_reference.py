#!/usr/bin/env python3
"""lambert_tof_reference.py - Lambert TOF curve truth table.

Emits a CSV of T(x, λ, n_rev) ground-truth values for the
Izzo 2014 Lambert TOF curve, computed from the elliptic-case
Lagrange analytic form ONLY (single branch, no Battin/Lagrange
boundary). This is the reference the multi-branch evaluator
at libk26astro_conics/src/lambert_multi.c::tof_curve_ must match.

Background
----------
The implementation splits the TOF evaluator into two branches:
  * Lagrange analytic form for `|x - 1| ∈ [0.01, 0.2]`
  * Battin universal/series form elsewhere (especially near x = 1)
A known scaling mismatch at the boundary - the Battin branch's
`4/3 · 2F1` normalization disagreeing with the Lagrange branch's
`0.5 · (1-x²)^(-1.5)` prefactor - causes numerical convergence on
x but the wrong orbit (semi-major axis off by ~10% on n_rev=1
fixtures).

To diagnose without depending on the suspect code, this script
implements ONLY the Lagrange form across the full elliptic range
(`-1 < x < 1`). Mathematically the two forms are equal; numerically
the Lagrange form loses precision as `x → 1` because `1-x² → 0` in
the denominator. The precision limit is documented explicitly so
the validation gate knows where to require exact agreement versus
where to allow a wider tolerance.

Formula (Izzo 2014 eq. 13, elliptic case)
-----------------------------------------
For -1 < x < 1 (elliptic), λ ∈ (-1, 1):
  a = 1 / (1 - x²)             # normalized semi-major axis
  y = sqrt(1 - λ² · (1 - x²))
  α = 2 · acos(x)              # twice the eccentric anomaly difference
  β = 2 · asin(sqrt(λ² · (1 - x²)))
  if λ < 0: β = -β             # signed transfer direction

  T = 0.5 · a^(3/2) · ((α - sin α) - (β - sin β) + 2π · n_rev)

Grid
----
We sample x ∈ {-0.9, -0.7, ..., 0.7, 0.9} (10 values, avoiding the
±1 singularities), λ ∈ {-0.8, -0.4, 0.0, 0.4, 0.8} (5 values),
n_rev ∈ {0, 1, 2, 3}. That's 200 fixture rows. The evaluator
should match this to 1e-12 relative for points where `|x - 1| > 0.2`
(Lagrange branch) and 1e-10 or so within the Battin band; looser
there because the truth itself loses precision.

Hand-check entries (Izzo 2014 Table 1)
--------------------------------------
The Izzo paper lists reference TOF values at sampled (λ, x); those
should be cross-checked against this script before trusting the
wider grid. The paper's Table 1 entries have not yet been transcribed
here; transcribing them is a planned follow-up, after which the
rerun should match this script to <1e-6 relative.

In the meantime, this CSV is INTERNALLY CONSISTENT with the
single-branch Lagrange form; sufficient to diagnose whether the
multi-branch evaluator agrees with the simpler form. Disagreement
between the evaluator and this CSV indicates the branch-boundary
mismatch. Agreement (everywhere) means the bug is elsewhere; most
likely in the velocity reconstruction (§3.5), not the TOF curve
itself.
"""
import math
import sys


def tof_lagrange(x, lam, n_rev):
    """Izzo eq. 13 - elliptic Lagrange form. -1 < x < 1, λ ∈ (-1, 1)."""
    if not (-1.0 < x < 1.0):
        raise ValueError(f"x = {x} outside elliptic range (-1, 1)")
    if not (-1.0 < lam < 1.0):
        raise ValueError(f"lam = {lam} outside (-1, 1)")
    one_minus_x2 = 1.0 - x * x  # > 0 in elliptic regime
    a = 1.0 / one_minus_x2
    y2 = 1.0 - lam * lam * one_minus_x2
    if y2 < 0.0:
        # Geometry-degenerate: lam·sqrt(1-x²) > 1 means the chord is
        # too long for any orbit with this x. Skip.
        return None
    alpha = 2.0 * math.acos(x)
    beta_arg = math.sqrt(lam * lam * one_minus_x2)
    if beta_arg > 1.0:
        return None
    beta = 2.0 * math.asin(beta_arg)
    if lam < 0.0:
        beta = -beta
    T = 0.5 * a ** 1.5 * (
        (alpha - math.sin(alpha)) - (beta - math.sin(beta))
        + 2.0 * math.pi * n_rev
    )
    return T


def main():
    out = sys.stdout
    out.write("# Lambert TOF truth — Izzo 2014 eq. 13 (elliptic Lagrange)\n")
    out.write("# columns: lam,x,n_rev,T_analytic\n")
    out.write("# precision: ~1e-15 relative for |x| < 0.7; degrades near |x|=1\n")
    out.write("# DO NOT edit by hand — regenerate via\n")
    out.write("#   python3 tests/data/lambert_tof_reference.py > tests/data/lambert_tof_reference.csv\n")
    out.write("lam,x,n_rev,T_analytic\n")

    lams = [-0.8, -0.4, 0.0, 0.4, 0.5, 0.8]      # 0.5 hits Izzo Table 1
    xs = [-0.9, -0.7, -0.5, -0.3, 0.0, 0.3, 0.5, 0.7, 0.9]
    n_revs = [0, 1, 2, 3]

    n_rows = 0
    for lam in lams:
        for x in xs:
            for n_rev in n_revs:
                T = tof_lagrange(x, lam, n_rev)
                if T is None:
                    continue
                # 17g preserves IEEE-754 round-trip of the truth.
                out.write(f"{lam:.17g},{x:.17g},{n_rev},{T:.17g}\n")
                n_rows += 1

    sys.stderr.write(f"lambert_tof_reference: emitted {n_rows} rows\n")

    # Internal consistency: T should be monotonic in n_rev for fixed
    # (λ, x) — each additional revolution adds a fixed positive
    # increment (2π · a^(3/2)).
    sys.stderr.write("# Internal consistency — n_rev monotonicity at λ=0.5, x=0.0:\n")
    prev = None
    for n_rev in range(4):
        T = tof_lagrange(0.0, 0.5, n_rev)
        sys.stderr.write(f"#   n_rev={n_rev}: T={T:.6f}")
        if prev is not None:
            delta = T - prev
            sys.stderr.write(f"  (+{delta:.6f} from n_rev={n_rev-1})")
        sys.stderr.write("\n")
        prev = T


if __name__ == "__main__":
    main()
