/* k26astro_core/frame.h — reference frame registry + transforms.
 *
 * Frames are first-class nominal types. A position carries a frame
 * tag; passing an ECEF position where ICRF is expected is a
 * compile-time error in KFL (the type system tracks the subtype) and
 * a runtime check in C (k26astro_frame_transform refuses cross-tagged
 * arithmetic without explicit conversion).
 *
 * Built-in frames (v0.1):
 *   ICRF        — IAU 2009 fundamental; barycentric origin
 *   HCRF        — heliocentric; parallel to ICRF in orientation
 *   ECI         — J2000 mean equator/equinox; Earth-centred
 *   GCRS        — IAU 2000; Earth centre-of-mass origin
 *   ECEF        — Earth-fixed; rotates with Earth (WGS84)
 *   MARS_FIXED  — Mars body-fixed (IAU 2018 rotation model + 2019 correction)
 *   LUNAR_PA    — Moon principal-axis frame
 *
 * Inertial vs body-fixed
 * ----------------------
 * Inertial frames (ICRF/HCRF/ECI/GCRS) differ only in origin —
 * orientation is identical (modulo precession/nutation at the
 * micro-arcsecond level, which v0.1 collapses to identity). A
 * inertial→inertial frame transform is an identity rotation; any
 * origin offset is handled by the application's body bookkeeping.
 *
 * Body-fixed frames (ECEF/MARS_FIXED/LUNAR_PA) rotate with their host
 * body. The transform from any inertial frame depends on the body's
 * attitude at the requested epoch, which lives in `libk26astro_body`.
 * If the body is not yet registered, inertial-to-body-fixed transforms
 * return K26A_FRAME_E_BODY_NOT_REGISTERED so consumers can detect the
 * gap.
 */
#ifndef K26ASTRO_CORE_FRAME_H
#define K26ASTRO_CORE_FRAME_H

#include "k26astro_core/epoch.h"
#include "k26astro_core/pos.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    K26A_FRAME_INVALID    = 0,
    K26A_FRAME_ICRF       = 1,
    K26A_FRAME_HCRF       = 2,
    K26A_FRAME_ECI        = 3,
    K26A_FRAME_GCRS       = 4,
    K26A_FRAME_ECEF       = 5,
    K26A_FRAME_MARS_FIXED = 6,
    K26A_FRAME_LUNAR_PA   = 7,
    /* User-registered frames start at this floor so the enum can be
     * extended without colliding. */
    K26A_FRAME_USER_BASE  = 1024
} K26AstroFrameId;

typedef enum {
    K26A_FRAME_KIND_INERTIAL   = 0,
    K26A_FRAME_KIND_BODY_FIXED = 1
} K26AstroFrameKind;

typedef struct {
    K26AstroFrameId   id;
    const char       *name;        /* e.g. "ICRF" — stable storage */
    K26AstroFrameKind kind;
    K26AstroFrameId   parent_inertial;  /* parent for body_fixed; self for inertial */
    /* For body_fixed: body identifier as a KFL ident (e.g. "earth").
     * Resolution to a K26AstroBody happens at transform time via
     * libk26astro_body's attitude lookup. NULL for inertial. */
    const char       *body_name;
} K26AstroFrameInfo;

/* Frame-tagged position. C-layer guard against passing one frame's
 * value where another is expected. */
typedef struct {
    K26AstroPos     p;
    K26AstroFrameId frame_id;
} K26AstroPosTagged;

/* ---- Error codes ----------------------------------------------- */
#define K26A_FRAME_OK                       0
#define K26A_FRAME_E_INVALID_INPUT          1
#define K26A_FRAME_E_UNKNOWN_FRAME          2
#define K26A_FRAME_E_BODY_NOT_REGISTERED    3
#define K26A_FRAME_E_UNSUPPORTED_PAIR       4
#define K26A_FRAME_E_REGISTRY_FULL          5

/* ---- Frame registry lookup ------------------------------------- */
/* Returns NULL if the id isn't registered. Built-ins are registered
 * automatically on first call to any frame API. */
const K26AstroFrameInfo *k26astro_frame_info(K26AstroFrameId id);

/* Lookup by name (e.g. "ICRF" / "ecef" — case-insensitive). Returns
 * K26A_FRAME_INVALID if not found. */
K26AstroFrameId k26astro_frame_by_name(const char *name);

/* User-frame registration. Pass an id in the K26A_FRAME_USER_BASE+
 * range. `name` and `body_name` must outlive the process (static
 * literals or arena-allocated). Returns 0 on success, non-zero on
 * collision / table full / invalid range. */
int k26astro_frame_register(K26AstroFrameId id,
                            const char *name,
                            K26AstroFrameKind kind,
                            K26AstroFrameId parent_inertial,
                            const char *body_name);

/* Clear user-registered frames (built-ins persist). */
void k26astro_frame_clear_user(void);

/* ---- Transforms ------------------------------------------------ */

/* Transform `in` to the `to` frame at epoch `at`. Writes the result
 * to `out`. Same handle for in/out is supported. Returns one of the
 * K26A_FRAME_E_* codes (0 = success).
 *
 * v0.1 behaviour:
 *   inertial → inertial : identity rotation (origins differ but the
 *                         caller is responsible for translating between
 *                         barycentre / heliocentre / geocentre via
 *                         their body's ephemeris position)
 *   inertial ↔ body-fixed : not supported until libk26astro_body
 *                           registers the body attitude; returns
 *                           K26A_FRAME_E_BODY_NOT_REGISTERED. */
int k26astro_frame_transform(K26AstroPosTagged       *out,
                             const K26AstroPosTagged *in,
                             K26AstroFrameId          to,
                             const K26AstroEpoch     *at);

/* Diagnostic: return 1 if a transform from `a` to `b` is supported by
 * the v0.1 frame layer (inertial↔inertial only). Used by callers
 * that want to gate behaviour rather than try-and-handle the error
 * code. */
int k26astro_frame_can_transform(K26AstroFrameId a, K26AstroFrameId b);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_CORE_FRAME_H */
