/* test_detect_smoke.c — sensor opaque construction, slot binding,
 * destructor protocol. Mirrors the vehicle payload-slot gate
 * in libk26astro_vehicle but exercises the detect-sensor slot in
 * place. */

#include "k26astro_detect/detect.h"

#include "k26astro_vehicle/vehicle.h"
#include "k26astro_body/body.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    K26AstroBody body;
    memset(&body, 0, sizeof(body));
    K26AstroVehicle *v = k26astro_vehicle_new();
    assert(v != NULL);
    k26astro_vehicle_bind_body(v, &body);

    /* Construct one of each kind. */
    K26AstroDetectSensor *ir = k26astro_detect_sensor_new_ir(
        /* aperture_m */ 0.3,
        /* integration_s */ 0.1,
        /* passband_lo_um */ 3.0,
        /* passband_hi_um */ 12.0,
        /* throughput */ 0.5,
        /* snr_threshold */ 5.0);
    assert(ir != NULL);
    assert(k26astro_detect_sensor_kind(ir) == K26ASTRO_DETECT_SENSOR_IR);

    K26AstroDetectSensor *radar = k26astro_detect_sensor_new_radar(
        /* p_tx_W */ 1.0e3,
        /* g_tx_db */ 35.0,
        /* g_rx_db */ 35.0,
        /* freq_Hz */ 10.0e9,
        /* loss_sys_db */ 4.0,
        /* bandwidth_Hz */ 1.0e6,
        /* T_sys_K */ 290.0,
        /* noise_figure */ 3.0,
        /* snr_threshold */ 13.0);
    assert(radar != NULL);
    assert(k26astro_detect_sensor_kind(radar) == K26ASTRO_DETECT_SENSOR_RADAR);

    K26AstroDetectSensor *lidar = k26astro_detect_sensor_new_lidar(
        /* pulse_energy_J */ 1.0e-3,
        /* wavelength_nm */ 1064.0,
        /* aperture_rx_m */ 0.2,
        /* atmospheric_tx */ 1.0,
        /* detector_efficiency */ 0.6,
        /* snr_threshold */ 5.0);
    assert(lidar != NULL);
    assert(k26astro_detect_sensor_kind(lidar) == K26ASTRO_DETECT_SENSOR_LIDAR);

    /* Invalid inputs must return NULL. */
    assert(k26astro_detect_sensor_new_ir(-1.0, 0.1, 3.0, 12.0, 0.5, 5.0) == NULL);
    assert(k26astro_detect_sensor_new_ir(0.3, -1.0, 3.0, 12.0, 0.5, 5.0) == NULL);
    assert(k26astro_detect_sensor_new_ir(0.3, 0.1, 12.0, 3.0, 0.5, 5.0) == NULL);
    assert(k26astro_detect_sensor_new_radar(-1.0, 35.0, 35.0, 1e9, 4.0, 1e6, 290.0, 3.0, 13.0) == NULL);
    assert(k26astro_detect_sensor_new_lidar(0.0, 1064.0, 0.2, 1.0, 0.6, 5.0) == NULL);

    /* Vehicle slot binding: attach all three into the PAYLOAD slot. */
    assert(k26astro_detect_sensor_attach(ir, v)    == 0);
    assert(k26astro_detect_sensor_attach(radar, v) == 0);
    assert(k26astro_detect_sensor_attach(lidar, v) == 0);

    /* Sensor-first destruction — slot invalidates in place through
     * the payload back-reference set by attach. */
    k26astro_detect_sensor_destroy(ir);
    k26astro_detect_sensor_destroy(radar);
    k26astro_detect_sensor_destroy(lidar);

    /* Vehicle teardown sees three dead slots; should not
     * dereference. */
    k26astro_vehicle_destroy(v);

    /* ---- Reverse order: vehicle destroyed first --------------- */
    K26AstroBody body2;
    memset(&body2, 0, sizeof(body2));
    K26AstroVehicle *v2 = k26astro_vehicle_new();
    k26astro_vehicle_bind_body(v2, &body2);

    /* Vehicle-first: the payload base's on_owner_destroy path blanks
     * the sensor's owner back-reference during vehicle teardown, so
     * the sensor's own destructor does not invalidate a freed
     * vehicle's slot. */
    K26AstroDetectSensor *ir2 = k26astro_detect_sensor_new_ir(
        0.3, 0.1, 3.0, 12.0, 0.5, 5.0);
    assert(k26astro_detect_sensor_attach(ir2, v2) == 0);
    k26astro_vehicle_destroy(v2);
    k26astro_detect_sensor_destroy(ir2);

    /* ---- NULL safety ----------------------------------------- */
    k26astro_detect_sensor_destroy(NULL);
    assert(k26astro_detect_sensor_kind(NULL) == 0);

    printf("test_detect_smoke: OK\n");
    return 0;
}
