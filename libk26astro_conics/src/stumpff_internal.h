/* stumpff_internal.h — universal-variable Stumpff functions C(z), S(z).
 *
 * Private to libk26astro_conics; consumed by kepler.c, kepler_edge.c,
 * lambert.c, lambert_multi.c. The switch on |z| picks between the
 * trigonometric (z > 0, elliptic), hyperbolic (z < 0), and Maclaurin
 * series (|z| ≈ 0) forms; the Maclaurin band avoids cancellation in
 * the (1 - cos√z)/z evaluation at small z.
 *
 * The exact branch thresholds (±1e-6) are inherited from the K26
 * pre-migration impl (libk26astro_body/elements.c through 2026-05-22);
 * tightening to 1e-4 has been considered but the cosine-cancellation
 * loss at z ≈ 1e-5 is observable in test_kepler_edge. */
#ifndef K26ASTRO_CONICS_STUMPFF_INTERNAL_H
#define K26ASTRO_CONICS_STUMPFF_INTERNAL_H

double k26astro_conics_stumpff_C(double z);
double k26astro_conics_stumpff_S(double z);

#endif
