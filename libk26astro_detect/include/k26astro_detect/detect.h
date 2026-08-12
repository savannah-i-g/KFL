/* k26astro_detect/detect.h — tactical detection.
 *
 * libk26astro_detect provides per-vehicle tactical detection
 * geometry — passive infrared, monostatic radar, photon-counting
 * lidar, faceted RCS — plus the counter-detection helper that
 * scores how loud a vehicle's own active emissions make it.
 *
 * The lib composes against:
 *   - libk26astro_vehicle       (generic PAYLOAD slot)
 *   - libk26astro_body          (body-frame geometry)
 *   - libk26astro_ephem         (frame transformations)
 *   - libk26astro_rt            (world tick + RNG)
 *   - libk26compute             (K26CRng, K26V3)
 *
 * Surface temperature is a caller-supplied input (set by the
 * platform's thermal model).
 *
 * The umbrella header pulls in every sub-surface so consumers
 * include a single header. Sub-surfaces remain individually
 * includable for callers that want only one capability.
 *
 * Determinism: noise generators thread a K26CRng* explicitly.
 * Passing NULL yields mean-only deterministic output; passing a
 * world-seeded RNG yields REFERENCED-mode bit-identical replay.
 *
 * Stated limitations:
 *   - IR signature assumes greybody / Lambertian emission. Sharply
 *     spectrally-structured plumes are first-order modelled.
 *   - Radar RCS is geometric-optics faceted; the resonant /
 *     subwavelength regime is not modelled.
 *   - Lidar atmospheric attenuation is a scalar caller-supplied
 *     factor; spectrally-structured atmospheres are not modelled.
 *   - Aperture diffraction sets angular resolution but does not
 *     bound the integrated-flux measurement; consumers needing
 *     to discriminate point sources from extended targets should
 *     track angular size separately.
 *
 * References:
 *   Hudson 1969          Infrared System Engineering
 *   Skolnik 2008         Radar Handbook, 3rd ed.
 *   Knott et al. 2004    Radar Cross Section, 2nd ed.
 *   Measures 1992        Laser Remote Sensing
 *   "There Ain't No Stealth In Space" — synthesis on W. Chung's
 *     Atomic Rockets site (projectrho.com)
 */
#ifndef K26ASTRO_DETECT_DETECT_H
#define K26ASTRO_DETECT_DETECT_H

#include "k26astro_detect/detect_consts.h"
#include "k26astro_detect/signature.h"
#include "k26astro_detect/detect_ir_passive.h"
#include "k26astro_detect/detect_active_rf.h"
#include "k26astro_detect/detect_active_lidar.h"
#include "k26astro_detect/counter_detect.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations. */
struct K26AstroVehicle;

/* ---- Detect-sensor opaque ------------------------------------- *
 *
 * A K26AstroDetectSensor is the per-vehicle handle bound into
 * libk26astro_vehicle's generic PAYLOAD slot. Each handle carries
 * one sensor's specifications (modality + aperture + integration
 * + pass-band) — multiple sensors per vehicle add to the slot
 * array.
 *
 * The sensor opaque is distinct from navigation-class sensors
 * (star tracker, IMU, magnetometer, ...), which measure the
 * platform's own state; detection sensors observe other
 * platforms. */
typedef struct K26AstroDetectSensor K26AstroDetectSensor;

typedef enum {
    K26ASTRO_DETECT_SENSOR_IR     = 1,
    K26ASTRO_DETECT_SENSOR_RADAR  = 2,
    K26ASTRO_DETECT_SENSOR_LIDAR  = 3
} K26AstroDetectSensorKind;

/* Construct a passive-IR detection sensor.
 *
 *   aperture_m         optical aperture diameter
 *   integration_s      dwell / integration time
 *   passband_lo_um     lower pass-band edge
 *   passband_hi_um     upper pass-band edge
 *   throughput         optical throughput × QE in [0, 1]
 *   snr_threshold      detection trigger
 *
 * Returns NULL on invalid inputs or allocation failure. */
K26AstroDetectSensor *k26astro_detect_sensor_new_ir(
    double aperture_m,
    double integration_s,
    double passband_lo_um,
    double passband_hi_um,
    double throughput,
    double snr_threshold);

/* Construct an active-radar detection sensor. */
K26AstroDetectSensor *k26astro_detect_sensor_new_radar(
    double p_tx_W,
    double g_tx_db,
    double g_rx_db,
    double freq_Hz,
    double loss_sys_db,
    double bandwidth_Hz,
    double T_sys_K,
    double noise_figure,
    double snr_threshold);

/* Construct an active-lidar detection sensor. atmospheric_tx and
 * detector_efficiency values <= 0 (or > 1) are stored as 1.0
 * (vacuum path / ideal detector defaults). */
K26AstroDetectSensor *k26astro_detect_sensor_new_lidar(
    double pulse_energy_J,
    double wavelength_nm,
    double aperture_rx_m,
    double atmospheric_tx,
    double detector_efficiency,
    double snr_threshold);

/* Destroy a detection sensor. Unlinks from the bound vehicle's
 * payload slot before freeing. Safe on NULL. */
void k26astro_detect_sensor_destroy(K26AstroDetectSensor *s);

/* Bind the sensor into a vehicle's PAYLOAD array slot. Returns 0 on
 * success, -1 on NULL inputs or allocation failure. */
int k26astro_detect_sensor_attach(K26AstroDetectSensor *s,
                                  struct K26AstroVehicle *v);

/* Sensor accessors. */
K26AstroDetectSensorKind
k26astro_detect_sensor_kind(const K26AstroDetectSensor *s);

#ifdef __cplusplus
}
#endif

#endif /* K26ASTRO_DETECT_DETECT_H */
