/* libk26plot — Cairo-backed scientific plotting for K26.
 *
 * This library renders 2-D scientific figures from caller-supplied
 * series data. The supported series kinds are line, scatter, errorbar,
 * pre-binned histogram, box plot (precomputed quartiles, whiskers, and
 * outliers), and heatmap (2-D matrix with linear color scale). Per-axis
 * log-scale is engaged via config flags. Output formats are PNG (image
 * surface) and SVG (vector); both render from a single Cairo draw pass
 * for byte-stable cross-format agreement. Defaults are publication
 * grade: sans-serif typeface, visible tick marks, "nice" tick steps
 * (1 / 2 / 5 x 10^n for linear, integer powers of 10 for log), 5%
 * data-range padding, optional grid, palette borrowed from matplotlib's
 * tab10.
 *
 * Determinism: same K26PlotConfig + same series produces the same PNG
 * bytes and the same SVG bytes (modulo Cairo's own version metadata).
 *
 * Not provided: contour plots, 3-D rendering, multiple subplots in a
 * single figure, polar / radial coordinates, animated frames. */
#ifndef K26PLOT_H
#define K26PLOT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define K26PLOT_LIB_VERSION "0.2.0"

typedef enum {
    K26P_OK          = 0,
    K26P_ERR_INVAL,
    K26P_ERR_OOM,
    K26P_ERR_IO,        /* Cairo couldn't write surface to path */
    K26P_ERR_INTERNAL,
} K26PStatus;

const char *k26plot_status_str(K26PStatus s);

/* Series kinds — extend with care; downstream switches on this. */
typedef enum {
    K26P_LINE      = 0,
    K26P_SCATTER   = 1,
    K26P_ERRORBAR  = 2,
    /* Pre-binned histogram: xs holds N+1 edges, ys holds N counts. */
    K26P_HISTOGRAM = 3,
    /* Box plot: caller pre-computes quartile statistics per box;
     * see the box_* fields below. */
    K26P_BOX       = 4,
    /* Heatmap: 2-D row-major matrix with linear color scale; see
     * the heat_* fields below. */
    K26P_HEATMAP   = 5,
} K26PSeriesKind;

typedef struct {
    K26PSeriesKind kind;
    const char    *label;            /* legend text; NULL → no legend entry */

    /* Common point-data — caller-owned, must outlive k26plot_render.
     * Used by line / scatter / errorbar / histogram (and by box plot
     * for the per-box x positions stored in `xs`, with `n` = box count). */
    const double  *xs;
    const double  *ys;
    const double  *yerr;             /* errorbar only; NULL otherwise */
    size_t         n;                /* count of (x,y) pairs / bins / boxes */

    /* Style. RGB in [0, 1]; if all three are 0 the renderer picks the
     * next default-palette color for this series index. */
    double r, g, b;
    double linewidth;                /* line/errorbar/box outline; 0 → default 1.5 */
    double point_size;               /* scatter/errorbar/outliers; 0 → default 4.0 */

    /* ---- Box plot extras (kind == K26P_BOX) ----
     * `n` is the box count. `xs` is the per-box x position. The five
     * box_* arrays each have length n. Outliers are flat — parallel
     * arrays of (x, y) coords across all boxes; n_outliers can be 0.
     * Caller pre-computes Q1/median/Q3 + whisker extents + outliers
     * (see libk26compute's k26c_stats_quantile for the standard path). */
    const double  *box_q1;
    const double  *box_median;
    const double  *box_q3;
    const double  *box_whisker_lo;
    const double  *box_whisker_hi;
    const double  *outliers_x;
    const double  *outliers_y;
    size_t         n_outliers;

    /* ---- Heatmap extras (kind == K26P_HEATMAP) ----
     * Row-major: `heat_data[row * heat_cols + col]`. The grid is
     * placed in axis-space spanning [heat_x0, heat_x1] × [heat_y0,
     * heat_y1] — caller controls position and aspect. Color scale is
     * linear from heat_vmin → heat_vmax; pass NaN for either to
     * auto-fit from the data's min/max. */
    const double  *heat_data;
    size_t         heat_rows;
    size_t         heat_cols;
    double         heat_x0, heat_x1;
    double         heat_y0, heat_y1;
    double         heat_vmin, heat_vmax;
} K26PSeries;

typedef struct {
    int width_px, height_px;

    const char *title;
    const char *xlabel, *ylabel;

    /* Set both NaN to auto-fit a 5% padded range from the data; set a
     * finite value to pin that bound. (NaN/finite mix is allowed — pin
     * one end, auto-fit the other.) */
    double xmin, xmax;
    double ymin, ymax;

    int show_grid;
    int show_legend;

    /* Log-scale axes — non-zero engages log10 mapping for that axis.
     * Data points must be strictly positive on a log axis (the
     * renderer drops non-positive values rather than NaNing). Tick
     * generation switches to integer powers of 10. */
    int log_x;
    int log_y;
} K26PlotConfig;

/* Initialise to publication-grade defaults: 800×600, no labels, auto
 * range, grid + legend on. The caller fills in title / labels / data
 * after this call. */
void k26plot_config_defaults(K26PlotConfig *cfg);

/* Render every series into one figure. `png_path` and `svg_path` are
 * each optional (NULL skips that format). Returns K26P_OK on success;
 * any I/O failure on either path returns K26P_ERR_IO without
 * attempting the other format. */
K26PStatus k26plot_render(const K26PlotConfig *cfg,
                          const K26PSeries *series, size_t n_series,
                          const char *png_path,
                          const char *svg_path);

/* Render every series onto a caller-provided cairo context. Used by
 * callers that render a plot live into an existing surface, in-place,
 * without a temporary PNG/SVG file. The function clips to the
 * destination rectangle, translates, draws, and restores the cairo
 * state. `cr` is owned by the caller; this routine does not destroy
 * it. Forward-declared via the cairo_t typedef so callers do not
 * need to include cairo.h transitively. */
typedef struct _cairo cairo_t;
K26PStatus k26plot_draw_to_cairo(cairo_t *cr,
                                  int dst_x, int dst_y,
                                  int dst_w, int dst_h,
                                  const K26PlotConfig *cfg,
                                  const K26PSeries *series, size_t n_series);

#ifdef __cplusplus
}
#endif

#endif /* K26PLOT_H */
