/* fpu_pin.c — pin FPU state at grav-state init for cross-platform
 * determinism.
 *
 * The astro suite's determinism contract requires:
 *   - rounding mode = FE_TONEAREST (the IEEE-754 default; pinned in
 *     case some upstream library — e.g. some FFT libs — changed it)
 *   - DAZ (denormals-as-zero) = OFF on x86
 *   - FTZ (flush-to-zero) = OFF on x86
 *
 * DAZ/FTZ shortcuts speed denormal arithmetic by zeroing tiny values,
 * which breaks bit-exact reproducibility across CPUs whose MXCSR
 * defaults differ. We force the standard-compliant behaviour.
 *
 * This file is silently no-op on non-x86; aarch64's FZ flag has a
 * similar effect but is not commonly set by libraries, so we don't
 * touch it for now. */
#include <fenv.h>

#if defined(__x86_64__) || defined(__i386__)
#include <pmmintrin.h>
#endif

void k26astro_grav_fpu_pin(void)
{
    fesetround(FE_TONEAREST);
#if defined(__x86_64__) || defined(__i386__)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
#endif
}
