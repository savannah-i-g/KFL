/* de441_naif_ids.h - named NAIF SPICE body codes for the inner-SS
 * gate.
 *
 * Per Acton 1996, "Ancillary Data Services of NASA's Navigation and
 * Ancillary Information Facility", Planet. Space Sci. 44:65: the
 * canonical SPICE id scheme has planet barycentres at 1..9 and the
 * actual planet centres at 100*barycentre+99 (Mercury = 199,
 * Venus = 299, Earth = 399, Mars = 499). Sun = 10.
 *
 * IMPORTANT: DE441 SPK kernel coverage: the JPL-distributed DE441
 * SPK contains entries for ids {10, 1..9, 301, 399}. The planet-
 * centre ids 199, 299, 499 are NOT present because for Mercury and
 * Venus (no significant moons) the planet centre is identical to
 * the barycentre, and for Mars (Phobos+Deimos mass ~1.6e-8 of
 * Mars's mass: sub-metre offset at 1 AU) the barycentre is
 * accurate to far better than our 10/50 km gate tolerance. So we
 * use the BARYCENTRE ids for Mercury, Venus, Mars; Earth uses the
 * planet centre (399) which IS present alongside the Earth-Moon
 * barycentre (3).
 *
 * Earth-specific: id 3 is the Earth-Moon barycentre and id 399 is
 * Earth's centre proper. The integrator wants the planet itself,
 * so we use 399; the difference is ~4670 km on Earth's diameter,
 * well above our 10 km gate tolerance.
 *
 * Reference: NAIF integer codes catalogue at
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
 * (reproduced here for offline reproducibility; we treat these as a
 * cited table rather than magic numbers). */
#ifndef K26_TESTS_ASTRO_DE441_NAIF_IDS_H
#define K26_TESTS_ASTRO_DE441_NAIF_IDS_H

#define K26_NAIF_SUN        10
#define K26_NAIF_MERCURY    1     /* Mercury barycentre (= centre) */
#define K26_NAIF_VENUS      2     /* Venus barycentre (= centre) */
#define K26_NAIF_EMB        3     /* Earth-Moon barycentre */
#define K26_NAIF_EARTH      399   /* Earth centre (NOT EMB = 3) */
#define K26_NAIF_MOON       301   /* Moon centre (for EMB correction) */
#define K26_NAIF_MARS       4     /* Mars barycentre (~ centre, ~1.6e-8) */

#endif /* K26_TESTS_ASTRO_DE441_NAIF_IDS_H */
