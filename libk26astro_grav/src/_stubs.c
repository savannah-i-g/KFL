/* _stubs.c — weak-default symbols so the .a links cleanly before the
 * full integrator + perturbation set lands.
 *
 * Each function here is overridden by its real implementation when
 * the corresponding source file is compiled into the .a. As more
 * sources land (wisdom_holman.c, ias15.c, rk_wrapper.c, perturb_*.c),
 * the linker resolves to the real symbol; this file's stubs become
 * dead code in the archive.
 *
 * Currently active stubs are tracked in the comments below. */
#include "k26astro_grav/grav.h"
#include "k26astro_grav/forces.h"
#include "k26astro_grav/perturb.h"

/* step_wh has a real impl in wisdom_holman.c; no stub. */

/* step_ias15 + ias15_reset have real impls in ias15.c. */

/* step_rk has a real impl in rk_wrapper.c; no stub needed. */

/* perturb_j2 / perturb_srp / perturb_gr_ppn1 / srp_shadow_test
 * have real impls in perturb_*.c. */

/* advise_step + close_encounter have real impls in their own files. */
