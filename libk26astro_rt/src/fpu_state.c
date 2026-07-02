/* fpu_state.c — save / pin / restore FPU mode across world lifecycle.
 *
 * The astro suite's PORTABLE determinism contract requires:
 *   - rounding mode = FE_TONEAREST
 *   - DAZ (denormals-as-zero) = OFF on x86
 *   - FTZ (flush-to-zero)     = OFF on x86
 *
 * libk26astro_grav pins per-state at grav_state_init; we wrap with a
 * save/restore pair so the caller's FPU mode survives a world's
 * lifetime regardless of whether the caller cared about it. */
#include "world_internal.h"

#include <fenv.h>

#if defined(__x86_64__) || defined(__i386__)
#include <pmmintrin.h>
#endif

void k26astro_fpu_save_and_pin(K26AstroFPUState *fs)
{
    if (!fs) return;
    fs->saved_round_mode = fegetround();
#if defined(__x86_64__) || defined(__i386__)
    fs->saved_daz = _MM_GET_DENORMALS_ZERO_MODE();
    fs->saved_ftz = _MM_GET_FLUSH_ZERO_MODE();
#else
    fs->saved_daz = 0;
    fs->saved_ftz = 0;
#endif

    fesetround(FE_TONEAREST);
#if defined(__x86_64__) || defined(__i386__)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
#endif

    fs->applied = 1;
}

void k26astro_fpu_restore(const K26AstroFPUState *fs)
{
    if (!fs || !fs->applied) return;
    fesetround(fs->saved_round_mode);
#if defined(__x86_64__) || defined(__i386__)
    _MM_SET_DENORMALS_ZERO_MODE(fs->saved_daz);
    _MM_SET_FLUSH_ZERO_MODE(fs->saved_ftz);
#endif
}
