/* detect_sensor.c — K26AstroDetectSensor opaque + vehicle
 * composition wiring.
 *
 * The opaque carries the per-sensor specification (modality +
 * kind-specific spec). One vehicle may carry many detection
 * sensors; each binds into the vehicle's generic PAYLOAD array
 * slot through the shared K26AstroPayload base (payload.h). Vehicle
 * teardown reaches back through the base's on_owner_destroy
 * callback, so no per-kind vehicle hook is needed. */

#include "k26astro_detect/detect.h"
#include "k26astro_vehicle/vehicle.h"
#include "k26astro_vehicle/payload.h"
#include "k26astro_defense/defense_kinds.h"

#include <stdlib.h>
#include <string.h>

/* Per-kind specification — tagged union via the public kind enum.
 * The K26AstroPayload base is the first member so a sensor pointer
 * downcasts to the base for vehicle composition. */
struct K26AstroDetectSensor {
    K26AstroPayload          base;
    K26AstroDetectSensorKind kind;
    union {
        struct {
            double aperture_m;
            double integration_s;
            double passband_lo_um;
            double passband_hi_um;
            double throughput;
            double snr_threshold;
        } ir;
        struct {
            double p_tx_W;
            double g_tx_db;
            double g_rx_db;
            double freq_Hz;
            double loss_sys_db;
            double bandwidth_Hz;
            double T_sys_K;
            double noise_figure;
            double snr_threshold;
        } radar;
        struct {
            double pulse_energy_J;
            double wavelength_nm;
            double aperture_rx_m;
            double atmospheric_tx;
            double detector_efficiency;
            double snr_threshold;
        } lidar;
    } spec;
};

/* ---- Constructors --------------------------------------------- */

K26AstroDetectSensor *k26astro_detect_sensor_new_ir(
    double aperture_m,
    double integration_s,
    double passband_lo_um,
    double passband_hi_um,
    double throughput,
    double snr_threshold)
{
    if (aperture_m <= 0.0 || integration_s <= 0.0) return NULL;
    if (passband_hi_um <= passband_lo_um) return NULL;
    if (throughput <= 0.0 || snr_threshold <= 0.0) return NULL;

    K26AstroDetectSensor *s = (K26AstroDetectSensor *)calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->base.kind       = K26ASTRO_DEFENSE_KIND_DETECT_SENSOR;
    s->base.owner_slot = -1;
    s->base.generation = 1;
    s->kind = K26ASTRO_DETECT_SENSOR_IR;
    s->spec.ir.aperture_m      = aperture_m;
    s->spec.ir.integration_s   = integration_s;
    s->spec.ir.passband_lo_um  = passband_lo_um;
    s->spec.ir.passband_hi_um  = passband_hi_um;
    s->spec.ir.throughput      = throughput;
    s->spec.ir.snr_threshold   = snr_threshold;
    return s;
}

K26AstroDetectSensor *k26astro_detect_sensor_new_radar(
    double p_tx_W,
    double g_tx_db,
    double g_rx_db,
    double freq_Hz,
    double loss_sys_db,
    double bandwidth_Hz,
    double T_sys_K,
    double noise_figure,
    double snr_threshold)
{
    if (p_tx_W <= 0.0 || freq_Hz <= 0.0 || bandwidth_Hz <= 0.0) return NULL;
    if (T_sys_K <= 0.0 || snr_threshold <= 0.0) return NULL;

    K26AstroDetectSensor *s = (K26AstroDetectSensor *)calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->base.kind       = K26ASTRO_DEFENSE_KIND_DETECT_SENSOR;
    s->base.owner_slot = -1;
    s->base.generation = 1;
    s->kind = K26ASTRO_DETECT_SENSOR_RADAR;
    s->spec.radar.p_tx_W         = p_tx_W;
    s->spec.radar.g_tx_db        = g_tx_db;
    s->spec.radar.g_rx_db        = g_rx_db;
    s->spec.radar.freq_Hz        = freq_Hz;
    s->spec.radar.loss_sys_db    = loss_sys_db;
    s->spec.radar.bandwidth_Hz   = bandwidth_Hz;
    s->spec.radar.T_sys_K        = T_sys_K;
    s->spec.radar.noise_figure   = (noise_figure < 1.0) ? 1.0 : noise_figure;
    s->spec.radar.snr_threshold  = snr_threshold;
    return s;
}

K26AstroDetectSensor *k26astro_detect_sensor_new_lidar(
    double pulse_energy_J,
    double wavelength_nm,
    double aperture_rx_m,
    double atmospheric_tx,
    double detector_efficiency,
    double snr_threshold)
{
    if (pulse_energy_J <= 0.0 || wavelength_nm <= 0.0) return NULL;
    if (aperture_rx_m <= 0.0 || snr_threshold <= 0.0) return NULL;

    K26AstroDetectSensor *s = (K26AstroDetectSensor *)calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->base.kind       = K26ASTRO_DEFENSE_KIND_DETECT_SENSOR;
    s->base.owner_slot = -1;
    s->base.generation = 1;
    s->kind = K26ASTRO_DETECT_SENSOR_LIDAR;
    s->spec.lidar.pulse_energy_J      = pulse_energy_J;
    s->spec.lidar.wavelength_nm       = wavelength_nm;
    s->spec.lidar.aperture_rx_m       = aperture_rx_m;
    s->spec.lidar.atmospheric_tx      =
        (atmospheric_tx <= 0.0)   ? 1.0 :
        (atmospheric_tx >  1.0)   ? 1.0 : atmospheric_tx;
    s->spec.lidar.detector_efficiency =
        (detector_efficiency <= 0.0) ? 1.0 :
        (detector_efficiency >  1.0) ? 1.0 : detector_efficiency;
    s->spec.lidar.snr_threshold       = snr_threshold;
    return s;
}

void k26astro_detect_sensor_destroy(K26AstroDetectSensor *s)
{
    if (!s) return;
    k26astro_payload_unlink_(&s->base);
    free(s);
}

/* Bind the sensor into a vehicle's PAYLOAD array slot. Records the
 * back-reference so a subsequent sensor-first destroy nulls the
 * vehicle slot in place. Returns 0 on success, -1 on NULL inputs
 * or allocation failure. */
int k26astro_detect_sensor_attach(K26AstroDetectSensor *s,
                                  struct K26AstroVehicle *v)
{
    if (!s) return -1;
    return k26astro_payload_bind(&s->base, v);
}

K26AstroDetectSensorKind
k26astro_detect_sensor_kind(const K26AstroDetectSensor *s)
{
    return s ? s->kind : (K26AstroDetectSensorKind)0;
}
