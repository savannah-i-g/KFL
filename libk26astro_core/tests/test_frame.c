/* test_frame.c — frame registry + inertial-only transforms (v0.1).
 *
 * Acceptance:
 *   - Built-in frames are queryable by id and by name (case-insens).
 *   - User-frame registration succeeds for ids in the
 *     K26A_FRAME_USER_BASE+ range and rejects collisions / built-in
 *     range / NULL inputs.
 *   - Inertial→inertial transforms succeed and preserve the
 *     position (identity rotation in v0.1).
 *   - Inertial→body-fixed transforms return E_BODY_NOT_REGISTERED. */
#include "k26astro_core/frame.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* ---- Lookups ----------------------------------------------- */
    const K26AstroFrameInfo *icrf = k26astro_frame_info(K26A_FRAME_ICRF);
    assert(icrf != NULL);
    assert(icrf->id == K26A_FRAME_ICRF);
    assert(strcmp(icrf->name, "ICRF") == 0);
    assert(icrf->kind == K26A_FRAME_KIND_INERTIAL);

    const K26AstroFrameInfo *ecef = k26astro_frame_info(K26A_FRAME_ECEF);
    assert(ecef != NULL);
    assert(ecef->kind == K26A_FRAME_KIND_BODY_FIXED);
    assert(ecef->body_name != NULL);
    assert(strcmp(ecef->body_name, "earth") == 0);

    /* Case-insensitive lookup by name. */
    assert(k26astro_frame_by_name("ICRF") == K26A_FRAME_ICRF);
    assert(k26astro_frame_by_name("icrf") == K26A_FRAME_ICRF);
    assert(k26astro_frame_by_name("Ecef") == K26A_FRAME_ECEF);
    assert(k26astro_frame_by_name("not_a_frame") == K26A_FRAME_INVALID);

    /* ---- User registration ------------------------------------ */
    /* Lower id range refused. */
    int rc = k26astro_frame_register(K26A_FRAME_ICRF, "DUP",
                                      K26A_FRAME_KIND_INERTIAL,
                                      K26A_FRAME_ICRF, NULL);
    assert(rc == K26A_FRAME_E_INVALID_INPUT);

    /* User-range registration succeeds. */
    rc = k26astro_frame_register(K26A_FRAME_USER_BASE, "SUN_FIXED",
                                  K26A_FRAME_KIND_BODY_FIXED,
                                  K26A_FRAME_HCRF, "sun");
    assert(rc == K26A_FRAME_OK);
    assert(k26astro_frame_by_name("SUN_FIXED") == K26A_FRAME_USER_BASE);

    /* Duplicate-id rejected. */
    rc = k26astro_frame_register(K26A_FRAME_USER_BASE, "OTHER",
                                  K26A_FRAME_KIND_INERTIAL,
                                  K26A_FRAME_ICRF, NULL);
    assert(rc == K26A_FRAME_E_INVALID_INPUT);

    /* Duplicate-name rejected. */
    rc = k26astro_frame_register(K26A_FRAME_USER_BASE + 1, "sun_fixed",
                                  K26A_FRAME_KIND_INERTIAL,
                                  K26A_FRAME_ICRF, NULL);
    assert(rc == K26A_FRAME_E_INVALID_INPUT);

    /* Clear user frames; built-ins persist. */
    k26astro_frame_clear_user();
    assert(k26astro_frame_by_name("SUN_FIXED") == K26A_FRAME_INVALID);
    assert(k26astro_frame_info(K26A_FRAME_ICRF) != NULL);

    /* ---- Transforms ------------------------------------------- */
    K26AstroPosTagged in_pos;
    in_pos.p = k26astro_pos_from_m(1.5e11, 0.0, 0.0);   /* ~1 AU */
    in_pos.frame_id = K26A_FRAME_ICRF;

    K26AstroPosTagged out_pos;
    K26AstroEpoch at = k26astro_epoch_j2000_tt();

    /* Inertial → inertial succeeds (identity in v0.1). */
    rc = k26astro_frame_transform(&out_pos, &in_pos, K26A_FRAME_HCRF, &at);
    assert(rc == K26A_FRAME_OK);
    assert(out_pos.frame_id == K26A_FRAME_HCRF);
    assert(out_pos.p.sx == in_pos.p.sx);
    assert(out_pos.p.lx == in_pos.p.lx);

    /* Inertial → body-fixed not yet supported. */
    rc = k26astro_frame_transform(&out_pos, &in_pos, K26A_FRAME_ECEF, &at);
    assert(rc == K26A_FRAME_E_BODY_NOT_REGISTERED);

    assert(k26astro_frame_can_transform(K26A_FRAME_ICRF, K26A_FRAME_HCRF) == 1);
    assert(k26astro_frame_can_transform(K26A_FRAME_ICRF, K26A_FRAME_ECEF) == 0);

    /* Round-trip: ICRF → HCRF → ICRF preserves the position. */
    K26AstroPosTagged a = in_pos;
    K26AstroPosTagged b, c;
    rc = k26astro_frame_transform(&b, &a, K26A_FRAME_HCRF, &at); assert(rc == 0);
    rc = k26astro_frame_transform(&c, &b, K26A_FRAME_ICRF, &at); assert(rc == 0);
    assert(c.frame_id == K26A_FRAME_ICRF);
    assert(c.p.sx == a.p.sx);
    assert(c.p.lx == a.p.lx);

    printf("test_frame: OK (registry + inertial transforms + body-fixed stub)\n");
    return 0;
}
