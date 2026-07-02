/* frame.c — frame registry + inertial-only transforms (v0.1).
 *
 * The built-in table is populated lazily on first call. User-
 * registered frames append to the same table; the v0.1 ceiling is
 * 64 entries, generous given the typical astro set + a handful of
 * user-declared frames per scene.
 *
 * Body-fixed transforms (ECEF/MARS_FIXED/LUNAR_PA) are documented as
 * not-implemented and return K26A_FRAME_E_BODY_NOT_REGISTERED. The
 * API surface here is stable; once libk26astro_body registers
 * attitude providers the body-fixed transforms become available
 * without breaking callers. */
#include "k26astro_core/frame.h"

#include <stddef.h>
#include <string.h>
#include <strings.h>

#define K26A_FRAME_REGISTRY_MAX 64

static K26AstroFrameInfo g_table[K26A_FRAME_REGISTRY_MAX];
static int               g_count = 0;
static int               g_builtins_loaded = 0;

static void load_builtins_(void)
{
    if (g_builtins_loaded) return;
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_ICRF, .name = "ICRF",
        .kind = K26A_FRAME_KIND_INERTIAL,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = NULL
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_HCRF, .name = "HCRF",
        .kind = K26A_FRAME_KIND_INERTIAL,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = NULL
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_ECI, .name = "ECI",
        .kind = K26A_FRAME_KIND_INERTIAL,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = NULL
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_GCRS, .name = "GCRS",
        .kind = K26A_FRAME_KIND_INERTIAL,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = NULL
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_ECEF, .name = "ECEF",
        .kind = K26A_FRAME_KIND_BODY_FIXED,
        .parent_inertial = K26A_FRAME_GCRS, .body_name = "earth"
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_MARS_FIXED, .name = "MARS_FIXED",
        .kind = K26A_FRAME_KIND_BODY_FIXED,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = "mars"
    };
    g_table[g_count++] = (K26AstroFrameInfo){
        .id = K26A_FRAME_LUNAR_PA, .name = "LUNAR_PA",
        .kind = K26A_FRAME_KIND_BODY_FIXED,
        .parent_inertial = K26A_FRAME_ICRF, .body_name = "moon"
    };
    g_builtins_loaded = 1;
}

const K26AstroFrameInfo *k26astro_frame_info(K26AstroFrameId id)
{
    load_builtins_();
    if (id == K26A_FRAME_INVALID) return NULL;
    for (int i = 0; i < g_count; i++) {
        if (g_table[i].id == id) return &g_table[i];
    }
    return NULL;
}

K26AstroFrameId k26astro_frame_by_name(const char *name)
{
    if (!name) return K26A_FRAME_INVALID;
    load_builtins_();
    for (int i = 0; i < g_count; i++) {
        if (strcasecmp(g_table[i].name, name) == 0) return g_table[i].id;
    }
    return K26A_FRAME_INVALID;
}

int k26astro_frame_register(K26AstroFrameId id,
                            const char *name,
                            K26AstroFrameKind kind,
                            K26AstroFrameId parent_inertial,
                            const char *body_name)
{
    if (!name) return K26A_FRAME_E_INVALID_INPUT;
    if (id < K26A_FRAME_USER_BASE) return K26A_FRAME_E_INVALID_INPUT;
    load_builtins_();
    for (int i = 0; i < g_count; i++) {
        if (g_table[i].id == id) return K26A_FRAME_E_INVALID_INPUT;
        if (strcasecmp(g_table[i].name, name) == 0)
            return K26A_FRAME_E_INVALID_INPUT;
    }
    if (g_count >= K26A_FRAME_REGISTRY_MAX) return K26A_FRAME_E_REGISTRY_FULL;
    g_table[g_count].id              = id;
    g_table[g_count].name            = name;
    g_table[g_count].kind            = kind;
    g_table[g_count].parent_inertial = parent_inertial;
    g_table[g_count].body_name       = body_name;
    g_count++;
    return K26A_FRAME_OK;
}

void k26astro_frame_clear_user(void)
{
    /* Re-load builtins from scratch: clear, then re-add. The
     * lazy-init flag isn't reset because the next call will short-
     * circuit on count > 0 — instead we directly run the loader. */
    g_count = 0;
    g_builtins_loaded = 0;
    load_builtins_();
}

int k26astro_frame_can_transform(K26AstroFrameId a, K26AstroFrameId b)
{
    const K26AstroFrameInfo *ia = k26astro_frame_info(a);
    const K26AstroFrameInfo *ib = k26astro_frame_info(b);
    if (!ia || !ib) return 0;
    return ia->kind == K26A_FRAME_KIND_INERTIAL
        && ib->kind == K26A_FRAME_KIND_INERTIAL;
}

int k26astro_frame_transform(K26AstroPosTagged       *out,
                             const K26AstroPosTagged *in,
                             K26AstroFrameId          to,
                             const K26AstroEpoch     *at)
{
    (void)at;   /* unused in v0.1 inertial-only path */
    if (!out || !in) return K26A_FRAME_E_INVALID_INPUT;
    const K26AstroFrameInfo *from_info = k26astro_frame_info(in->frame_id);
    const K26AstroFrameInfo *to_info   = k26astro_frame_info(to);
    if (!from_info || !to_info) return K26A_FRAME_E_UNKNOWN_FRAME;

    if (from_info->kind != K26A_FRAME_KIND_INERTIAL
     || to_info->kind   != K26A_FRAME_KIND_INERTIAL) {
        return K26A_FRAME_E_BODY_NOT_REGISTERED;
    }

    /* Inertial→inertial: identity rotation in v0.1. The caller is
     * responsible for body-centric origin shifts (e.g. moving a
     * position from heliocentric to geocentric via the Earth's
     * heliocentric position). */
    out->p        = in->p;
    out->frame_id = to;
    return K26A_FRAME_OK;
}
