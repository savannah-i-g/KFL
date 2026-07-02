/* builtin_register_astro.c: opaque-type + builtin registration
 * for the astro suite.
 *
 * Acts as the in-process fallback for embedders that do not ship
 * the .kflbi manifest path; the kflc driver calls into this single
 * registration function at startup. It mirrors the manifests at:
 *   - libk26astro_rt/kfl/k26astro_rt.kflbi
 *   - libk26astro_render/kfl/k26astro_render.kflbi
 *
 * Strings passed to the registry must outlive the compiler invocation;
 * we use static string literals, which satisfies that contract. */
#include "kflc.h"

void kflc_register_astro_builtins(void)
{
    /* ---- Opaque types ------------------------------------------- */
    (void)kflc_opaque_register("world",     "K26AstroWorld *");
    (void)kflc_opaque_register("starfield", "K26AstroStarfield *");
    /* Body handle is a pointer into the world's body array.
     * Lifetime semantics: a body pointer is realloc-stable across
     * step()/observe() calls (those don't add bodies) but NOT
     * across subsequent add_body() calls (the array can grow).
     * Inside for_each loops the iterator pointer is valid for one
     * iteration. The parse-time binding for `astro_body NAME` to
     * `propagate NAME` lookup uses an index slot maintained at
     * fn prologue, since pointers are not realloc-stable. */
    (void)kflc_opaque_register("body",      "K26AstroBody *");

    /* ---- libk26astro_rt builtins -------------------------------- */
    (void)kflc_register_builtin("astro_world_open",
                                "k26astro_world_create",              2);
    (void)kflc_register_builtin("astro_world_close",
                                "k26astro_world_destroy",             1);
    (void)kflc_register_builtin("astro_world_add_body",
                                "k26astro_world_add_body",            2);
    (void)kflc_register_builtin("astro_world_body_count",
                                "k26astro_world_body_count",          1);
    (void)kflc_register_builtin("astro_world_find_body",
                                "k26astro_world_find_body",           2);
    (void)kflc_register_builtin("astro_world_step",
                                "k26astro_world_step",                2);
    (void)kflc_register_builtin("astro_world_observe",
                                "k26astro_world_observe",             5);
    (void)kflc_register_builtin("astro_world_snapshot_save",
                                "k26astro_world_snapshot_save",       2);
    (void)kflc_register_builtin("astro_world_snapshot_load",
                                "k26astro_world_snapshot_load",       1);
    (void)kflc_register_builtin("astro_world_set_spin_hz",
                                "k26astro_scheduler_set_spin_hz",     2);
    (void)kflc_register_builtin("astro_world_set_render_hz",
                                "k26astro_scheduler_set_render_hz",   2);
    (void)kflc_register_builtin("astro_world_set_observer_mode",
                                "k26astro_world_set_observer_mode",   2);
    (void)kflc_register_builtin("astro_world_set_mercurius",
                                "k26astro_world_set_mercurius_factors", 3);
    (void)kflc_register_builtin("astro_world_set_ephem",
                                "k26astro_world_set_ephem",           2);
    /* Body-handle accessors. The "body" opaque is a K26AstroBody
     * pointer; these builtins surface the realloc-safe way to
     * obtain one from a world + index. */
    (void)kflc_register_builtin("astro_world_body_at",
                                "k26astro_world_body_at",             2);
    (void)kflc_register_builtin("astro_body_name",
                                "k26astro_body_name",                 1);

    /* ---- libk26astro_render builtins ---------------------------- */
    (void)kflc_register_builtin("astro_starfield",
                                "k26astro_render_starfield_install",  2);
    (void)kflc_register_builtin("astro_render_world",
                                "k26astro_render_install",            2);
    (void)kflc_register_builtin("astro_render_set_log_depth",
                                "k26astro_render_set_log_depth_range", 3);
}
