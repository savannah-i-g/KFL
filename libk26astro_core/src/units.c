/* units.c — placeholder TU.
 *
 * All conversions live as `static inline` in units.h so callers get
 * constant folding at the call site. This file exists so the linker
 * doesn't have to special-case header-only translation units, and so
 * future non-inlineable helpers (magnitude/bandpass conversions,
 * unit-tagged scalar arithmetic) have a natural home. */
#include "k26astro_core/units.h"
