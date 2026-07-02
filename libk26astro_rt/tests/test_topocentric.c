/* test_topocentric.c — TOPOCENTRIC observer mode gate.
 *
 * Verifies the libk26astro_atmos integration: setting an atmos
 * object via k26astro_world_set_atmos and switching the world
 * into TOPOCENTRIC mode applies refraction to the observer's
 * apparent direction.
 *
 * Two coverage axes:
 *   1. Sunrise refraction. Observer at Earth surface, looking at
 *      a target near the geometric horizon. With Earth-standard
 *      atmos attached, the apparent elevation should be lifted
 *      by ~34 arcmin (the Bennett horizon plateau).
 *   2. TOPOCENTRIC degrades to APPARENT when no atmos set.
 *      Without calling _set_atmos, TOPOCENTRIC mode should produce
 *      the same answer as APPARENT mode (no refraction).
 *
 * Acceptance:
 *   - apparent elevation lift at horizon-look: 30 to 40 arcmin
 *     (Bennett 1982 closed-form ± atmosphere n0 scale)
 *   - TOPOCENTRIC without atmos == APPARENT */
#include "k26astro_rt/world.h"
#include "k26astro_rt/observer.h"
#include "k26astro_atmos/atmos.h"
#include "k26astro_body/body.h"
#include "k26astro_core/consts.h"
#include "k26astro_core/pos.h"
#include "k26m3d.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define R_EARTH 6.371e6

static double elevation_(K26V3 dir, K26V3 zenith)
{
    double dn = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    double zn = sqrt(zenith.x*zenith.x + zenith.y*zenith.y
                   + zenith.z*zenith.z);
    double cs = (dir.x*zenith.x + dir.y*zenith.y + dir.z*zenith.z)
              / (dn * zn);
    if (cs >  1.0) cs =  1.0;
    if (cs < -1.0) cs = -1.0;
    return asin(cs);
}

static int observe_(K26AstroWorld *w, int tgt, int obs, K26V3 *out_dir)
{
    K26AstroPos out_pos;
    return k26astro_world_observe(w, tgt, obs, &out_pos, out_dir);
}

static int seed_earth_horizon_world_(K26AstroWorld *w,
                                      int *obs_i, int *tgt_i)
{
    /* Explicit Earth body at world origin acts as the observer's
     * parent (zenith reference). The multi-body extension reads
     * observer->parent_body_idx so observers on Moon, Mars, or
     * non-origin Earth-centric frames all work. */
    K26AstroBody earth;
    k26astro_body_init(&earth);
    strncpy(earth.name, "earth", sizeof earth.name - 1);
    earth.mass = 5.972e24; earth.gm = 3.986004418e14;
    earth.radius = R_EARTH;
    earth.pos  = k26astro_pos_from_m(0.0, 0.0, 0.0);
    earth.vel  = k26m3d_v3(0.0, 0.0, 0.0);
    int earth_i = k26astro_world_add_body(w, earth);

    /* Observer at Earth's surface on +x. Zenith direction is +x. */
    K26AstroBody observer;
    k26astro_body_init(&observer);
    strncpy(observer.name, "observer", sizeof observer.name - 1);
    observer.mass = 1.0; observer.gm = 0.0;
    observer.pos  = k26astro_pos_from_m(R_EARTH, 0.0, 0.0);
    observer.vel  = k26m3d_v3(0.0, 0.0, 0.0);  /* No aberration. */
    observer.parent_body_idx = earth_i;
    *obs_i = k26astro_world_add_body(w, observer);

    /* Target body high on +y (perpendicular to zenith → geometric
     * elevation = 0, i.e. exactly on horizon). Far enough away
     * that the geometric horizon angle is essentially 0. */
    K26AstroBody target;
    k26astro_body_init(&target);
    strncpy(target.name, "target", sizeof target.name - 1);
    target.mass = 1.0; target.gm = 0.0;
    target.pos  = k26astro_pos_from_m(R_EARTH, 1.0e9, 0.0);
    target.vel  = k26m3d_v3(0.0, 0.0, 0.0);
    *tgt_i = k26astro_world_add_body(w, target);
    return 0;
}

static int test_sunrise_refraction(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    K26AstroAtmos *atmos = k26astro_atmos_earth_standard();
    assert(atmos);
    int rc = k26astro_world_set_atmos(w, atmos);
    if (rc != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: set_atmos rc=%d\n", rc);
        return 1;
    }
    k26astro_world_set_observer_mode(w, K26ASTRO_OBS_TOPOCENTRIC);

    int obs_i, tgt_i;
    seed_earth_horizon_world_(w, &obs_i, &tgt_i);

    K26V3 dir;
    rc = observe_(w, tgt_i, obs_i, &dir);
    if (rc != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: observe rc=%d\n", rc);
        k26astro_atmos_destroy(atmos);
        k26astro_world_destroy(w);
        return 1;
    }

    K26V3 zenith = { 1.0, 0.0, 0.0 };
    double app_elev = elevation_(dir, zenith);
    double app_arcmin = app_elev * (180.0 * 60.0 / M_PI);
    fprintf(stderr,
        "TOPOCENTRIC horizon-look: apparent elevation = %.3f arcmin\n",
        app_arcmin);
    if (app_arcmin < 30.0 || app_arcmin > 40.0) {
        fprintf(stderr,
            "FAIL: horizon refraction %.3f arcmin out of [30, 40]\n",
            app_arcmin);
        k26astro_atmos_destroy(atmos);
        k26astro_world_destroy(w);
        return 1;
    }

    k26astro_world_destroy(w);
    k26astro_atmos_destroy(atmos);
    return 0;
}

static int test_topocentric_without_atmos(void)
{
    K26AstroWorld *w_top = k26astro_world_create(K26ASTRO_MODE_FAST,
                                                  K26ASTRO_COORDS_SECTOR_GRID);
    K26AstroWorld *w_app = k26astro_world_create(K26ASTRO_MODE_FAST,
                                                  K26ASTRO_COORDS_SECTOR_GRID);
    assert(w_top && w_app);
    k26astro_world_set_observer_mode(w_top, K26ASTRO_OBS_TOPOCENTRIC);
    k26astro_world_set_observer_mode(w_app, K26ASTRO_OBS_APPARENT);

    int obs_i_t, tgt_i_t, obs_i_a, tgt_i_a;
    seed_earth_horizon_world_(w_top, &obs_i_t, &tgt_i_t);
    seed_earth_horizon_world_(w_app, &obs_i_a, &tgt_i_a);

    K26V3 dir_top, dir_app;
    observe_(w_top, tgt_i_t, obs_i_t, &dir_top);
    observe_(w_app, tgt_i_a, obs_i_a, &dir_app);

    double diff = sqrt((dir_top.x - dir_app.x)*(dir_top.x - dir_app.x)
                     + (dir_top.y - dir_app.y)*(dir_top.y - dir_app.y)
                     + (dir_top.z - dir_app.z)*(dir_top.z - dir_app.z));
    fprintf(stderr,
        "TOPOCENTRIC w/o atmos vs APPARENT: |Δdir|=%.3e\n", diff);
    if (diff > 1.0e-12) {
        fprintf(stderr,
            "FAIL: TOPOCENTRIC (no atmos) should equal APPARENT, "
            "got diff %.3e\n", diff);
        k26astro_world_destroy(w_top);
        k26astro_world_destroy(w_app);
        return 1;
    }

    k26astro_world_destroy(w_top);
    k26astro_world_destroy(w_app);
    return 0;
}

/* Multi-body coverage: observer on Mars's surface, world origin
 * at heliocentric. Demonstrates that the topocentric refraction
 * uses position differencing against the named parent body, not
 * the world's coordinate origin. */
static int test_topocentric_multi_body(void)
{
    K26AstroWorld *w = k26astro_world_create(K26ASTRO_MODE_FAST,
                                              K26ASTRO_COORDS_SECTOR_GRID);
    assert(w);
    K26AstroAtmos *atmos = k26astro_atmos_earth_standard();
    assert(atmos);
    k26astro_world_set_atmos(w, atmos);
    k26astro_world_set_observer_mode(w, K26ASTRO_OBS_TOPOCENTRIC);

    /* Mars at a non-origin heliocentric position (~1.5 AU on +x). */
    const double AU = 1.495978707e11;
    const double R_MARS = 3.3895e6;
    K26AstroBody mars;
    k26astro_body_init(&mars);
    strncpy(mars.name, "mars", sizeof mars.name - 1);
    mars.mass = 6.4171e23; mars.gm = 4.282837e13;
    mars.radius = R_MARS;
    mars.pos = k26astro_pos_from_m(1.5 * AU, 0.0, 0.0);
    mars.vel = k26m3d_v3(0.0, 0.0, 0.0);
    int mars_i = k26astro_world_add_body(w, mars);

    /* Observer on Mars's +x surface (heliocentric +x of Mars centre). */
    K26AstroBody obs;
    k26astro_body_init(&obs);
    strncpy(obs.name, "obs", sizeof obs.name - 1);
    obs.mass = 1.0; obs.gm = 0.0;
    obs.pos = k26astro_pos_from_m(1.5 * AU + R_MARS, 0.0, 0.0);
    obs.vel = k26m3d_v3(0.0, 0.0, 0.0);
    obs.parent_body_idx = mars_i;
    int obs_i = k26astro_world_add_body(w, obs);

    /* Target on Mars's local horizon (+y, far). */
    K26AstroBody tgt;
    k26astro_body_init(&tgt);
    strncpy(tgt.name, "tgt", sizeof tgt.name - 1);
    tgt.mass = 1.0; tgt.gm = 0.0;
    tgt.pos = k26astro_pos_from_m(1.5 * AU + R_MARS, 1.0e9, 0.0);
    tgt.vel = k26m3d_v3(0.0, 0.0, 0.0);
    int tgt_i = k26astro_world_add_body(w, tgt);

    K26V3 dir;
    K26AstroPos out_pos;
    int rc = k26astro_world_observe(w, tgt_i, obs_i, &out_pos, &dir);
    if (rc != K26ASTRO_RT_OK) {
        fprintf(stderr, "FAIL: multi-body observe rc=%d\n", rc);
        k26astro_atmos_destroy(atmos);
        k26astro_world_destroy(w);
        return 1;
    }

    /* The local zenith at obs is +x (away from Mars). Apparent
     * elevation should be lifted by the Bennett ~34 arcmin (same
     * physics as the Earth horizon test; the parent-lookup is the
     * only thing being exercised). */
    K26V3 zenith = { 1.0, 0.0, 0.0 };
    double app_elev = elevation_(dir, zenith);
    double app_arcmin = app_elev * (180.0 * 60.0 / M_PI);
    fprintf(stderr,
        "TOPOCENTRIC multi-body (Mars surface): elev = %.3f arcmin\n",
        app_arcmin);
    if (app_arcmin < 30.0 || app_arcmin > 40.0) {
        fprintf(stderr,
            "FAIL: non-origin parent refraction %.3f arcmin out of "
            "[30, 40]\n", app_arcmin);
        k26astro_atmos_destroy(atmos);
        k26astro_world_destroy(w);
        return 1;
    }

    k26astro_atmos_destroy(atmos);
    k26astro_world_destroy(w);
    return 0;
}

int main(void)
{
    if (test_sunrise_refraction()) return 1;
    if (test_topocentric_without_atmos()) return 1;
    if (test_topocentric_multi_body()) return 1;
    fprintf(stderr, "test_topocentric: OK\n");
    return 0;
}
