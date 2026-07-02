/* internal.h — private definitions shared across libk26astro_body
 * translation units. Not installed.
 */
#ifndef K26ASTRO_BODY_INTERNAL_H
#define K26ASTRO_BODY_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "k26astro_body/star.h"

/* Shared catalogue struct shape. Defined here so the file-backed
 * loader (star_catalogue.c) and the baked-in subset
 * (star_catalogue_hip.c) agree on the layout exactly. */
struct K26AstroCatalogue {
    void                *map;
    size_t               map_size;
    const K26AstroStar  *records;
    size_t               n_records;
    double               epoch_jd_tt;
    uint64_t             catalogue_id;
    uint8_t              file_backed;
};

#endif /* K26ASTRO_BODY_INTERNAL_H */
