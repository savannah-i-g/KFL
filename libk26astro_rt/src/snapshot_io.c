/* snapshot_io.c — versioned binary save/load.
 *
 * Header (80 bytes):
 *   char     magic[8]       "K26ASNP\0"
 *   uint32_t version        1
 *   uint32_t endian_probe   0x01020304
 *   uint32_t flags          bit0=Q64_64, bit1=FAST_MODE, bit2=PORTABLE
 *   uint32_t body_count
 *   K26AstroEpoch epoch     i64 + dbl + u8 + 7 reserved = 24 bytes
 *   char     integrator[16] e.g. "wisdom_holman\0\0\0"
 *   uint8_t  reserved[8]
 *
 * Body records (field-by-field, no struct memcpy — adds field-reorder
 * resilience):
 *   u64 id; char name[32]; u32 kind;
 *   K26AstroPos (i64 sx,sy,sz; dbl lx,ly,lz);
 *   K26V3 vel (dbl x,y,z);
 *   dbl mass, gm, radius, polar_radius, j2;
 *   K26Quat attitude (dbl w,x,y,z);
 *   K26V3 omega;
 *   u32 attitude_mode;
 *   i32 rotation_model_id;
 *   u8 has_atmos;
 *   i32 ephem_naif_id;
 *   u8 on_rails;
 *   i32 parent_body_idx;
 */
#include "k26astro_rt/snapshot.h"
#include "world_internal.h"
#include "snapshot_internal.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Endian helpers — assume x86_64 host (little-endian native; flip
 * only on cross-platform read). */

int k26astro_snap_write_u16(FILE *f, uint16_t v)
{
    uint8_t b[2] = { (uint8_t)(v & 0xff), (uint8_t)((v >> 8) & 0xff) };
    return fwrite(b, 1, 2, f) == 2 ? 0 : -1;
}

int k26astro_snap_write_u32(FILE *f, uint32_t v)
{
    uint8_t b[4];
    for (int i = 0; i < 4; i++) b[i] = (uint8_t)((v >> (i * 8)) & 0xff);
    return fwrite(b, 1, 4, f) == 4 ? 0 : -1;
}

int k26astro_snap_write_u64(FILE *f, uint64_t v)
{
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (uint8_t)((v >> (i * 8)) & 0xff);
    return fwrite(b, 1, 8, f) == 8 ? 0 : -1;
}

int k26astro_snap_write_i64(FILE *f, int64_t v)
{
    return k26astro_snap_write_u64(f, (uint64_t)v);
}

int k26astro_snap_write_dbl(FILE *f, double v)
{
    uint64_t bits;
    memcpy(&bits, &v, sizeof bits);
    return k26astro_snap_write_u64(f, bits);
}

int k26astro_snap_write_buf(FILE *f, const void *p, size_t n)
{
    return fwrite(p, 1, n, f) == n ? 0 : -1;
}

static uint16_t bswap16_(uint16_t v) { return (uint16_t)((v >> 8) | (v << 8)); }
static uint32_t bswap32_(uint32_t v)
{
    return ((v >> 24) & 0xff)        | ((v >> 8)  & 0xff00) |
           ((v << 8)  & 0xff0000)    | ((v << 24) & 0xff000000);
}
static uint64_t bswap64_(uint64_t v)
{
    return ((uint64_t)bswap32_((uint32_t)(v & 0xffffffffu)) << 32) |
           (uint64_t)bswap32_((uint32_t)(v >> 32));
}

int k26astro_snap_read_u16(FILE *f, uint16_t *out, int swap)
{
    uint8_t b[2];
    if (fread(b, 1, 2, f) != 2) return -1;
    uint16_t v = (uint16_t)b[0] | ((uint16_t)b[1] << 8);
    *out = swap ? bswap16_(v) : v;
    return 0;
}

int k26astro_snap_read_u32(FILE *f, uint32_t *out, int swap)
{
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) return -1;
    uint32_t v = 0;
    for (int i = 0; i < 4; i++) v |= ((uint32_t)b[i]) << (i * 8);
    *out = swap ? bswap32_(v) : v;
    return 0;
}

int k26astro_snap_read_u64(FILE *f, uint64_t *out, int swap)
{
    uint8_t b[8];
    if (fread(b, 1, 8, f) != 8) return -1;
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= ((uint64_t)b[i]) << (i * 8);
    *out = swap ? bswap64_(v) : v;
    return 0;
}

int k26astro_snap_read_i64(FILE *f, int64_t *out, int swap)
{
    uint64_t u;
    int rc = k26astro_snap_read_u64(f, &u, swap);
    if (rc) return rc;
    *out = (int64_t)u;
    return 0;
}

int k26astro_snap_read_dbl(FILE *f, double *out, int swap)
{
    uint64_t bits;
    int rc = k26astro_snap_read_u64(f, &bits, swap);
    if (rc) return rc;
    memcpy(out, &bits, sizeof bits);
    return 0;
}

int k26astro_snap_read_buf(FILE *f, void *p, size_t n)
{
    return fread(p, 1, n, f) == n ? 0 : -1;
}

/* Body serialization ---------------------------------------------- */

static int write_body_(FILE *f, const K26AstroBody *b)
{
    /* Per-body type tag + flags. */
    if (k26astro_snap_write_u16(f, 0)) return -1; /* tag: base layout */
    if (k26astro_snap_write_u16(f, 0)) return -1; /* flags: reserved */

    if (k26astro_snap_write_u64(f, b->id)) return -1;
    if (k26astro_snap_write_buf(f, b->name, K26ASTRO_BODY_NAME_MAX)) return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->kind)) return -1;

    if (k26astro_snap_write_i64(f, b->pos.sx)) return -1;
    if (k26astro_snap_write_i64(f, b->pos.sy)) return -1;
    if (k26astro_snap_write_i64(f, b->pos.sz)) return -1;
    if (k26astro_snap_write_dbl(f, b->pos.lx)) return -1;
    if (k26astro_snap_write_dbl(f, b->pos.ly)) return -1;
    if (k26astro_snap_write_dbl(f, b->pos.lz)) return -1;

    if (k26astro_snap_write_dbl(f, b->vel.x)) return -1;
    if (k26astro_snap_write_dbl(f, b->vel.y)) return -1;
    if (k26astro_snap_write_dbl(f, b->vel.z)) return -1;

    if (k26astro_snap_write_dbl(f, b->mass))          return -1;
    if (k26astro_snap_write_dbl(f, b->gm))            return -1;
    if (k26astro_snap_write_dbl(f, b->radius))        return -1;
    if (k26astro_snap_write_dbl(f, b->polar_radius))  return -1;
    if (k26astro_snap_write_dbl(f, b->j2))            return -1;

    if (k26astro_snap_write_dbl(f, b->attitude.w)) return -1;
    if (k26astro_snap_write_dbl(f, b->attitude.x)) return -1;
    if (k26astro_snap_write_dbl(f, b->attitude.y)) return -1;
    if (k26astro_snap_write_dbl(f, b->attitude.z)) return -1;

    if (k26astro_snap_write_dbl(f, b->omega.x)) return -1;
    if (k26astro_snap_write_dbl(f, b->omega.y)) return -1;
    if (k26astro_snap_write_dbl(f, b->omega.z)) return -1;

    if (k26astro_snap_write_u32(f, (uint32_t)b->attitude_mode))     return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->rotation_model_id)) return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->has_atmos))         return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->ephem_naif_id))     return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->on_rails))          return -1;
    if (k26astro_snap_write_u32(f, (uint32_t)b->parent_body_idx))   return -1;

    return 0;
}

static int read_body_(FILE *f, int swap, K26AstroBody *b)
{
    uint16_t tag = 0, flags = 0;
    if (k26astro_snap_read_u16(f, &tag, swap))   return -1;
    if (k26astro_snap_read_u16(f, &flags, swap)) return -1;
    if (tag != 0) return -1;
    (void)flags;

    if (k26astro_snap_read_u64(f, &b->id, swap)) return -1;
    if (k26astro_snap_read_buf(f, b->name, K26ASTRO_BODY_NAME_MAX)) return -1;
    uint32_t kind = 0;
    if (k26astro_snap_read_u32(f, &kind, swap)) return -1;
    b->kind = (K26AstroBodyKind)kind;

    if (k26astro_snap_read_i64(f, &b->pos.sx, swap)) return -1;
    if (k26astro_snap_read_i64(f, &b->pos.sy, swap)) return -1;
    if (k26astro_snap_read_i64(f, &b->pos.sz, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->pos.lx, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->pos.ly, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->pos.lz, swap)) return -1;

    if (k26astro_snap_read_dbl(f, &b->vel.x, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->vel.y, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->vel.z, swap)) return -1;

    if (k26astro_snap_read_dbl(f, &b->mass, swap))         return -1;
    if (k26astro_snap_read_dbl(f, &b->gm, swap))           return -1;
    if (k26astro_snap_read_dbl(f, &b->radius, swap))       return -1;
    if (k26astro_snap_read_dbl(f, &b->polar_radius, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->j2, swap))           return -1;

    if (k26astro_snap_read_dbl(f, &b->attitude.w, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->attitude.x, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->attitude.y, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->attitude.z, swap)) return -1;

    if (k26astro_snap_read_dbl(f, &b->omega.x, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->omega.y, swap)) return -1;
    if (k26astro_snap_read_dbl(f, &b->omega.z, swap)) return -1;

    uint32_t att_mode = 0, has_atm = 0, on_rails = 0;
    int32_t  rot_id = 0, naif = 0, parent = 0;
    uint32_t tmp;
    if (k26astro_snap_read_u32(f, &att_mode, swap)) return -1;
    if (k26astro_snap_read_u32(f, &tmp, swap))      return -1;
    rot_id   = (int32_t)tmp;
    if (k26astro_snap_read_u32(f, &has_atm, swap))  return -1;
    if (k26astro_snap_read_u32(f, &tmp, swap))      return -1;
    naif     = (int32_t)tmp;
    if (k26astro_snap_read_u32(f, &on_rails, swap)) return -1;
    if (k26astro_snap_read_u32(f, &tmp, swap))      return -1;
    parent   = (int32_t)tmp;
    b->attitude_mode     = (K26AstroAttitudeMode)att_mode;
    b->rotation_model_id = rot_id;
    b->has_atmos         = (uint8_t)has_atm;
    b->ephem_naif_id     = naif;
    b->on_rails          = (uint8_t)on_rails;
    b->parent_body_idx   = parent;
    return 0;
}

/* Public API ------------------------------------------------------- */

int k26astro_world_snapshot_save(const K26AstroWorld *world, const char *path)
{
    if (!world || !path) return -K26ASTRO_RT_E_NULL;
    FILE *f = fopen(path, "wb");
    if (!f) return -K26ASTRO_RT_E_BAD_ARG;

    /* Header. */
    if (k26astro_snap_write_buf(f, K26ASTRO_SNAPSHOT_MAGIC, 8) ||
        k26astro_snap_write_u32(f, K26ASTRO_SNAPSHOT_VERSION_CURRENT) ||
        k26astro_snap_write_u32(f, K26ASTRO_SNAPSHOT_ENDIAN_PROBE))
    {
        fclose(f);
        return -K26ASTRO_RT_E_BAD_ARG;
    }
    uint32_t flags = 0;
    if (world->coord_mode == K26ASTRO_COORDS_Q64_64) flags |= K26ASTRO_SNAPSHOT_FLAG_Q64_64;
    if (world->mode       == K26ASTRO_MODE_FAST)     flags |= K26ASTRO_SNAPSHOT_FLAG_FAST_MODE;
    if (k26astro_snap_write_u32(f, flags)) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }
    if (k26astro_snap_write_u32(f, (uint32_t)world->grav.n_bodies)) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }

    /* Epoch. */
    if (k26astro_snap_write_i64(f, world->grav.t.days_since_J2000)) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }
    if (k26astro_snap_write_dbl(f, world->grav.t.seconds_of_day))   { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }
    uint8_t scale = world->grav.t.scale;
    if (fwrite(&scale, 1, 1, f) != 1) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }
    uint8_t pad[7] = {0};
    if (fwrite(pad, 1, 7, f) != 7) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }

    /* Integrator name (16 bytes, NUL-padded). */
    const char *integrator_name = "wisdom_holman";
    switch (world->grav.integrator) {
    case K26ASTRO_INTEGRATOR_WH:        integrator_name = "wisdom_holman"; break;
    case K26ASTRO_INTEGRATOR_IAS15:     integrator_name = "ias15"; break;
    case K26ASTRO_INTEGRATOR_VERLET:    integrator_name = "verlet"; break;
    case K26ASTRO_INTEGRATOR_RK4:       integrator_name = "rk4"; break;
    case K26ASTRO_INTEGRATOR_RK45:      integrator_name = "rk45"; break;
    case K26ASTRO_INTEGRATOR_MERCURIUS: integrator_name = "mercurius"; break;
    }
    char namebuf[16] = {0};
    size_t nlen = strlen(integrator_name);
    if (nlen > 15) nlen = 15;
    memcpy(namebuf, integrator_name, nlen);
    if (k26astro_snap_write_buf(f, namebuf, 16)) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }

    /* Reserved 8 bytes. */
    uint8_t reserved[8] = {0};
    if (k26astro_snap_write_buf(f, reserved, 8)) { fclose(f); return -K26ASTRO_RT_E_BAD_ARG; }

    /* Body records. */
    for (int i = 0; i < world->grav.n_bodies; i++) {
        if (write_body_(f, &world->grav.bodies[i]) != 0) {
            fclose(f);
            return -K26ASTRO_RT_E_BAD_ARG;
        }
    }

    if (fclose(f) != 0) return -K26ASTRO_RT_E_BAD_ARG;
    return K26ASTRO_RT_OK;
}

K26AstroWorld *k26astro_world_snapshot_load(const char *path)
{
    if (!path) return NULL;
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    char magic[8];
    if (k26astro_snap_read_buf(f, magic, 8) || memcmp(magic, K26ASTRO_SNAPSHOT_MAGIC, 8) != 0) {
        fclose(f); return NULL;
    }
    uint32_t version, probe;
    if (k26astro_snap_read_u32(f, &version, 0)) { fclose(f); return NULL; }
    if (k26astro_snap_read_u32(f, &probe,   0)) { fclose(f); return NULL; }
    int swap = (probe != K26ASTRO_SNAPSHOT_ENDIAN_PROBE);
    if (swap) version = bswap32_(version);
    if (version != K26ASTRO_SNAPSHOT_VERSION_CURRENT) { fclose(f); return NULL; }

    uint32_t flags = 0, n_bodies = 0;
    if (k26astro_snap_read_u32(f, &flags,    swap)) { fclose(f); return NULL; }
    if (k26astro_snap_read_u32(f, &n_bodies, swap)) { fclose(f); return NULL; }

    K26AstroWorldMode  mode = (flags & K26ASTRO_SNAPSHOT_FLAG_FAST_MODE)
                              ? K26ASTRO_MODE_FAST : K26ASTRO_MODE_PORTABLE;
    K26AstroCoordsMode cmod = (flags & K26ASTRO_SNAPSHOT_FLAG_Q64_64)
                              ? K26ASTRO_COORDS_Q64_64 : K26ASTRO_COORDS_SECTOR_GRID;

    K26AstroEpoch t;
    if (k26astro_snap_read_i64(f, &t.days_since_J2000, swap)) { fclose(f); return NULL; }
    if (k26astro_snap_read_dbl(f, &t.seconds_of_day,   swap)) { fclose(f); return NULL; }
    uint8_t scale_byte;
    if (fread(&scale_byte, 1, 1, f) != 1) { fclose(f); return NULL; }
    t.scale = scale_byte;
    uint8_t pad[7];
    if (fread(pad, 1, 7, f) != 7) { fclose(f); return NULL; }

    char namebuf[16];
    if (k26astro_snap_read_buf(f, namebuf, 16)) { fclose(f); return NULL; }
    namebuf[15] = '\0';

    uint8_t reserved[8];
    if (k26astro_snap_read_buf(f, reserved, 8)) { fclose(f); return NULL; }

    K26AstroWorld *world = k26astro_world_create(mode, cmod);
    if (!world) { fclose(f); return NULL; }
    world->grav.t = t;

    /* Pick integrator from name. */
    if      (strncmp(namebuf, "wisdom_holman", 15) == 0) k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_WH);
    else if (strncmp(namebuf, "ias15", 15) == 0)         k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_IAS15);
    else if (strncmp(namebuf, "verlet", 15) == 0)        k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_VERLET);
    else if (strncmp(namebuf, "rk4", 15) == 0)           k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_RK4);
    else if (strncmp(namebuf, "rk45", 15) == 0)          k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_RK45);
    else if (strncmp(namebuf, "mercurius", 15) == 0)     k26astro_grav_set_integrator(&world->grav, K26ASTRO_INTEGRATOR_MERCURIUS);

    /* Bodies. */
    for (uint32_t i = 0; i < n_bodies; i++) {
        K26AstroBody b;
        memset(&b, 0, sizeof b);
        if (read_body_(f, swap, &b) != 0) {
            k26astro_world_destroy(world);
            fclose(f);
            return NULL;
        }
        if (k26astro_world_add_body(world, b) < 0) {
            k26astro_world_destroy(world);
            fclose(f);
            return NULL;
        }
    }

    fclose(f);
    return world;
}
