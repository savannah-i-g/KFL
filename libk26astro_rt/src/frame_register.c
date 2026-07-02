/* frame_register.c — per-world frame registry helpers above
 * libk26astro_core's process-global frame registry.
 *
 * Built-in frames (ICRF/HCRF/ECI/GCRS/ECEF/MARS_FIXED/LUNAR_PA) are
 * always findable. User frames live in the per-world array; we also
 * push the registration into the process-global registry so that
 * existing libk26astro_body frame transforms find them (they query
 * the global registry, not the world). */
#include "k26astro_rt/world.h"
#include "world_internal.h"

#include <string.h>
#include <strings.h>

int k26astro_world_register_frame(K26AstroWorld *world,
                                   const K26AstroFrameInfo *info)
{
    if (!world || !info) return -K26ASTRO_RT_E_NULL;
    if (world->n_frames >= K26ASTRO_WORLD_FRAME_MAX)
        return -K26ASTRO_RT_E_FULL;

    /* Assign next user id (callers may also pass an explicit id in
     * `info->id`; we accept it if it's in the user range). */
    K26AstroFrameId assigned = info->id;
    if (assigned < K26A_FRAME_USER_BASE) {
        assigned = (K26AstroFrameId)(world->next_user_id++);
    } else if ((int)assigned >= world->next_user_id) {
        world->next_user_id = (int)assigned + 1;
    }

    K26AstroFrameInfo *slot = &world->frames[world->n_frames++];
    *slot = *info;
    slot->id = assigned;

    /* Mirror into the process-global registry so libk26astro_body /
     * libk26astro_ephem frame transforms see it. Best-effort; if the
     * global registry rejects (e.g. full) we still keep our local
     * copy. */
    (void)k26astro_frame_register(assigned, info->name, info->kind,
                                   info->parent_inertial, info->body_name);
    return (int)assigned;
}

int k26astro_world_find_frame(const K26AstroWorld *world, const char *name)
{
    if (!world || !name) return -K26ASTRO_RT_E_NULL;
    /* Built-ins first — the global registry handles them. */
    K26AstroFrameId built_in = k26astro_frame_by_name(name);
    if (built_in != K26A_FRAME_INVALID) return (int)built_in;
    /* Then per-world user frames. */
    for (int i = 0; i < world->n_frames; i++) {
        if (world->frames[i].name &&
            strcasecmp(world->frames[i].name, name) == 0)
            return (int)world->frames[i].id;
    }
    return -K26ASTRO_RT_E_NOT_FOUND;
}
