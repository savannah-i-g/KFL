/* libk26plot — Cairo-backed scientific plotting.
 *
 * Single-file implementation: status table, config defaults, tick
 * generation (linear + log), auto-extent computation, axis frame +
 * grid + title + labels, legend, per-kind series drawing (line /
 * scatter / errorbar / pre-binned histogram / box plot / heatmap),
 * and PNG + SVG export wrappers.
 *
 * The renderer is structured so the SAME draw routine is re-run for
 * each output format — PNG and SVG come from one logical figure
 * description, no format-specific drift. */

#include "k26plot.h"

#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Layout constants ---- */
#define MARGIN_LEFT     60.0
#define MARGIN_RIGHT    24.0
#define MARGIN_TOP      40.0
#define MARGIN_BOTTOM   54.0
#define LEGEND_PAD      8.0
#define LEGEND_LINE_LEN 18.0
#define DEFAULT_WIDTH   800
#define DEFAULT_HEIGHT  600

/* matplotlib's tab10 palette — the de-facto modern scientific
 * default. Series with all-zero (r,g,b) cycle through these by index. */
static const double DEFAULT_PALETTE[10][3] = {
    { 0.121, 0.466, 0.705 },        /* blue   */
    { 1.000, 0.498, 0.054 },        /* orange */
    { 0.172, 0.627, 0.172 },        /* green  */
    { 0.839, 0.152, 0.156 },        /* red    */
    { 0.580, 0.403, 0.741 },        /* purple */
    { 0.549, 0.337, 0.294 },        /* brown  */
    { 0.890, 0.466, 0.760 },        /* pink   */
    { 0.498, 0.498, 0.498 },        /* grey   */
    { 0.737, 0.741, 0.133 },        /* olive  */
    { 0.090, 0.745, 0.811 },        /* cyan   */
};

/* Heatmap color scale — viridis-approximating 4-stop gradient.
 * Perceptually monotonic, colourblind-safe, the standard scientific
 * default for sequential data. */
static const double HEAT_STOPS[4][4] = {
    /* t,    r,     g,     b */
    { 0.00, 0.267, 0.005, 0.329 },        /* dark purple */
    { 0.33, 0.230, 0.318, 0.546 },        /* blue-violet */
    { 0.66, 0.119, 0.601, 0.541 },        /* teal-green  */
    { 1.00, 0.993, 0.906, 0.144 },        /* yellow      */
};

/* ---- Status ---- */
const char *k26plot_status_str(K26PStatus s)
{
    switch (s) {
    case K26P_OK:           return "ok";
    case K26P_ERR_INVAL:    return "invalid argument";
    case K26P_ERR_OOM:      return "out of memory";
    case K26P_ERR_IO:       return "i/o failure (cairo write)";
    case K26P_ERR_INTERNAL: return "internal error";
    }
    return "unknown status";
}

/* ---- Defaults ---- */
void k26plot_config_defaults(K26PlotConfig *cfg)
{
    if (!cfg) return;
    memset(cfg, 0, sizeof *cfg);
    cfg->width_px    = DEFAULT_WIDTH;
    cfg->height_px   = DEFAULT_HEIGHT;
    cfg->title       = NULL;
    cfg->xlabel      = NULL;
    cfg->ylabel      = NULL;
    cfg->xmin = cfg->xmax = NAN;
    cfg->ymin = cfg->ymax = NAN;
    cfg->show_grid   = 1;
    cfg->show_legend = 1;
    cfg->log_x       = 0;
    cfg->log_y       = 0;
}

/* ---- Axis-scale transform ----
 * Linear axis: identity. Log axis: log10 of strictly-positive value;
 * non-positive values map to NaN (caller-side: data should already be
 * filtered, but defensive code below treats NaN cells as "skip"). */
static inline double scale_axis(double v, int log_axis)
{
    if (!log_axis) return v;
    return v > 0.0 ? log10(v) : NAN;
}

/* ---- "Nice" tick step: 1 / 2 / 5 × 10ⁿ (linear) ---- */
static double nice_step(double range, int target_ticks)
{
    if (target_ticks < 2) target_ticks = 2;
    if (range <= 0.0)     return 1.0;
    double raw  = range / (double)target_ticks;
    double mag  = pow(10.0, floor(log10(raw)));
    double norm = raw / mag;
    if      (norm < 1.5) return 1.0  * mag;
    else if (norm < 3.0) return 2.0  * mag;
    else if (norm < 7.0) return 5.0  * mag;
    else                 return 10.0 * mag;
}

static size_t compute_ticks_lin(double lo, double hi, int target,
                                double *out, size_t out_cap)
{
    if (out_cap == 0 || hi <= lo) return 0;
    double step  = nice_step(hi - lo, target);
    double first = ceil(lo / step) * step;
    if (fabs(first) < 1e-12 * step) first = 0.0;
    size_t n = 0;
    for (double t = first; t <= hi + step * 1e-9 && n < out_cap; t += step) {
        out[n++] = t;
    }
    return n;
}

/* For log axes: ticks at integer powers of 10 within the log-space
 * range. The values returned are in LOG-TRANSFORMED space (matching
 * how lo/hi are stored), but represent decades — render the label
 * as 10^t in data space. */
static size_t compute_ticks_log(double lo, double hi,
                                double *out, size_t out_cap)
{
    if (out_cap == 0 || hi <= lo) return 0;
    int e0 = (int)ceil(lo);
    int e1 = (int)floor(hi);
    if (e1 < e0) {
        /* Range spans <1 decade — fall back to linear ticks in log space. */
        return compute_ticks_lin(lo, hi, 5, out, out_cap);
    }
    size_t n = 0;
    for (int e = e0; e <= e1 && n < out_cap; e++) {
        out[n++] = (double)e;
    }
    return n;
}

static size_t compute_ticks(double lo, double hi, int log_axis,
                            int target_lin,
                            double *out, size_t out_cap)
{
    return log_axis ? compute_ticks_log(lo, hi, out, out_cap)
                    : compute_ticks_lin(lo, hi, target_lin, out, out_cap);
}

static void format_tick(double t, int log_axis, char *buf, size_t cap)
{
    if (log_axis) {
        double v = pow(10.0, t);
        snprintf(buf, cap, "%.6g", v);
    } else {
        snprintf(buf, cap, "%.6g", t);
    }
}

/* ---- Auto-extent computation ----
 * Walk every series, update lo/hi for x and y in DATA space (linear).
 * Errorbar series extend the y-range by ±yerr; histograms include the
 * 0 baseline; box plots use whisker extents; heatmaps use the
 * heat_{x0,x1,y0,y1} spatial bounds. Empty data falls back to [0, 1]. */
static void compute_extents(const K26PSeries *series, size_t n,
                            double *xmin_out, double *xmax_out,
                            double *ymin_out, double *ymax_out)
{
    int seen_x = 0, seen_y = 0;
    double xlo = 0, xhi = 1, ylo = 0, yhi = 1;

    for (size_t i = 0; i < n; i++) {
        const K26PSeries *s = &series[i];

        if (s->kind == K26P_HEATMAP) {
            double xs[2] = { s->heat_x0, s->heat_x1 };
            double ys[2] = { s->heat_y0, s->heat_y1 };
            for (int k = 0; k < 2; k++) {
                if (!seen_x) { xlo = xhi = xs[k]; seen_x = 1; }
                else { if (xs[k] < xlo) xlo = xs[k]; if (xs[k] > xhi) xhi = xs[k]; }
                if (!seen_y) { ylo = yhi = ys[k]; seen_y = 1; }
                else { if (ys[k] < ylo) ylo = ys[k]; if (ys[k] > yhi) yhi = ys[k]; }
            }
            continue;
        }

        if (!s->xs || s->n == 0) continue;

        size_t x_count = (s->kind == K26P_HISTOGRAM) ? s->n + 1 : s->n;
        for (size_t j = 0; j < x_count; j++) {
            double x = s->xs[j];
            if (!seen_x) { xlo = xhi = x; seen_x = 1; }
            else { if (x < xlo) xlo = x; if (x > xhi) xhi = x; }
        }

        if (s->kind == K26P_BOX) {
            /* Y range covers the whisker extents per box. */
            if (s->box_whisker_lo && s->box_whisker_hi) {
                for (size_t j = 0; j < s->n; j++) {
                    double a = s->box_whisker_lo[j];
                    double b = s->box_whisker_hi[j];
                    if (!seen_y) { ylo = a; yhi = b; seen_y = 1; }
                    if (a < ylo) ylo = a;
                    if (b > yhi) yhi = b;
                }
            }
            /* Outliers can extend further. */
            if (s->outliers_y) {
                for (size_t j = 0; j < s->n_outliers; j++) {
                    double y = s->outliers_y[j];
                    if (!seen_y) { ylo = yhi = y; seen_y = 1; }
                    else { if (y < ylo) ylo = y; if (y > yhi) yhi = y; }
                }
            }
            continue;
        }

        if (!s->ys) continue;
        for (size_t j = 0; j < s->n; j++) {
            double y = s->ys[j];
            if (!seen_y) { ylo = yhi = y; seen_y = 1; }
            else { if (y < ylo) ylo = y; if (y > yhi) yhi = y; }
            if (s->kind == K26P_ERRORBAR && s->yerr) {
                double e = s->yerr[j];
                if (y - e < ylo) ylo = y - e;
                if (y + e > yhi) yhi = y + e;
            }
            if (s->kind == K26P_HISTOGRAM) {
                if (0.0 < ylo) ylo = 0.0;
                if (0.0 > yhi) yhi = 0.0;
            }
        }
    }
    if (!seen_x) { xlo = 0; xhi = 1; }
    if (!seen_y) { ylo = 0; yhi = 1; }
    if (xlo == xhi) { xlo -= 0.5; xhi += 0.5; }
    if (ylo == yhi) { ylo -= 0.5; yhi += 0.5; }

    /* 5% padding both directions. */
    double xpad = (xhi - xlo) * 0.05;
    double ypad = (yhi - ylo) * 0.05;
    *xmin_out = xlo - xpad;
    *xmax_out = xhi + xpad;
    *ymin_out = ylo - ypad;
    *ymax_out = yhi + ypad;
}

/* Resolve final extents in TRANSFORMED axis space.
 * For log axes the returned bounds are log10-of-data values; coord
 * mapping below applies scale_axis() to data points before linear
 * placement, keeping the rest of the pipeline scale-agnostic. */
static void resolve_extents(const K26PlotConfig *cfg,
                            const K26PSeries *series, size_t n,
                            double *xmin, double *xmax,
                            double *ymin, double *ymax)
{
    double axmin, axmax, aymin, aymax;
    compute_extents(series, n, &axmin, &axmax, &aymin, &aymax);
    *xmin = isnan(cfg->xmin) ? axmin : cfg->xmin;
    *xmax = isnan(cfg->xmax) ? axmax : cfg->xmax;
    *ymin = isnan(cfg->ymin) ? aymin : cfg->ymin;
    *ymax = isnan(cfg->ymax) ? aymax : cfg->ymax;

    /* Log-transform bounds. Clamp non-positive to a small positive
     * before log to avoid -inf; data points still get filtered at
     * draw time via scale_axis returning NaN. */
    if (cfg->log_x) {
        if (*xmin <= 0) *xmin = (axmin > 0) ? axmin : 1e-9;
        if (*xmax <= *xmin) *xmax = *xmin * 10.0;
        *xmin = log10(*xmin);
        *xmax = log10(*xmax);
    }
    if (cfg->log_y) {
        if (*ymin <= 0) *ymin = (aymin > 0) ? aymin : 1e-9;
        if (*ymax <= *ymin) *ymax = *ymin * 10.0;
        *ymin = log10(*ymin);
        *ymax = log10(*ymax);
    }
    if (*xmax <= *xmin) *xmax = *xmin + 1.0;
    if (*ymax <= *ymin) *ymax = *ymin + 1.0;
}

/* ---- Coord mapping ----
 * Cairo's y axis grows downward, so y_to_pixel inverts the data y
 * before mapping into [top, bottom]. Both axes apply scale_axis()
 * for the log case so callers pass raw data values throughout. */
static inline double x_to_px(double x, int log_x,
                             double xmin, double xmax,
                             double left, double right)
{
    double t = scale_axis(x, log_x);
    if (isnan(t)) return left - 1e6;       /* off-canvas; clipped */
    return left + (t - xmin) / (xmax - xmin) * (right - left);
}
static inline double y_to_px(double y, int log_y,
                             double ymin, double ymax,
                             double top, double bottom)
{
    double t = scale_axis(y, log_y);
    if (isnan(t)) return bottom + 1e6;
    return bottom - (t - ymin) / (ymax - ymin) * (bottom - top);
}

/* ---- Resolve series color ----
 * If r==g==b==0, fall back to palette[idx % 10]. */
static void series_color(const K26PSeries *s, size_t idx,
                         double *cr, double *cg, double *cb)
{
    if (s->r == 0.0 && s->g == 0.0 && s->b == 0.0) {
        const double *p = DEFAULT_PALETTE[idx % 10];
        *cr = p[0]; *cg = p[1]; *cb = p[2];
    } else {
        *cr = s->r; *cg = s->g; *cb = s->b;
    }
}

/* ---- Draw axes (frame + grid + ticks + labels + title) ---- */
static void draw_axes(cairo_t *cr, double w, double h,
                      const K26PlotConfig *cfg,
                      double xmin, double xmax,
                      double ymin, double ymax)
{
    double left   = MARGIN_LEFT;
    double right  = w - MARGIN_RIGHT;
    double top    = MARGIN_TOP;
    double bottom = h - MARGIN_BOTTOM;

    double xticks[64];
    size_t nx = compute_ticks(xmin, xmax, cfg->log_x, 7, xticks, 64);
    double yticks[64];
    size_t ny = compute_ticks(ymin, ymax, cfg->log_y, 6, yticks, 64);

    cairo_select_font_face(cr, "sans-serif",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);

    /* Grid. Tick values are in transformed space; map directly via
     * the linear formula since we already log-transformed bounds. */
    if (cfg->show_grid) {
        cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.25);
        cairo_set_line_width(cr, 0.5);
        for (size_t i = 0; i < nx; i++) {
            double t = xticks[i];
            double px = left + (t - xmin) / (xmax - xmin) * (right - left);
            cairo_move_to(cr, px, top);
            cairo_line_to(cr, px, bottom);
        }
        for (size_t i = 0; i < ny; i++) {
            double t = yticks[i];
            double py = bottom - (t - ymin) / (ymax - ymin) * (bottom - top);
            cairo_move_to(cr, left,  py);
            cairo_line_to(cr, right, py);
        }
        cairo_stroke(cr);
    }

    /* Frame. */
    cairo_set_source_rgb(cr, 0.15, 0.15, 0.15);
    cairo_set_line_width(cr, 1.0);
    cairo_rectangle(cr, left, top, right - left, bottom - top);
    cairo_stroke(cr);

    /* Tick marks + numeric labels. */
    cairo_set_font_size(cr, 11.0);
    char buf[32];
    cairo_text_extents_t te;

    for (size_t i = 0; i < nx; i++) {
        double t = xticks[i];
        double px = left + (t - xmin) / (xmax - xmin) * (right - left);
        cairo_move_to(cr, px, bottom);
        cairo_line_to(cr, px, bottom + 5);
        cairo_stroke(cr);
        format_tick(t, cfg->log_x, buf, sizeof buf);
        cairo_text_extents(cr, buf, &te);
        cairo_move_to(cr, px - te.width / 2 - te.x_bearing, bottom + 18);
        cairo_show_text(cr, buf);
    }
    for (size_t i = 0; i < ny; i++) {
        double t = yticks[i];
        double py = bottom - (t - ymin) / (ymax - ymin) * (bottom - top);
        cairo_move_to(cr, left - 5, py);
        cairo_line_to(cr, left,     py);
        cairo_stroke(cr);
        format_tick(t, cfg->log_y, buf, sizeof buf);
        cairo_text_extents(cr, buf, &te);
        cairo_move_to(cr, left - 8 - te.width - te.x_bearing,
                      py + te.height / 2);
        cairo_show_text(cr, buf);
    }

    /* Title. */
    if (cfg->title && cfg->title[0]) {
        cairo_set_font_size(cr, 14.0);
        cairo_text_extents(cr, cfg->title, &te);
        cairo_move_to(cr,
                      (left + right) / 2 - te.width / 2 - te.x_bearing,
                      MARGIN_TOP - 14);
        cairo_show_text(cr, cfg->title);
    }
    /* X-axis label. */
    if (cfg->xlabel && cfg->xlabel[0]) {
        cairo_set_font_size(cr, 12.0);
        cairo_text_extents(cr, cfg->xlabel, &te);
        cairo_move_to(cr,
                      (left + right) / 2 - te.width / 2 - te.x_bearing,
                      h - 14);
        cairo_show_text(cr, cfg->xlabel);
    }
    /* Y-axis label (rotated). */
    if (cfg->ylabel && cfg->ylabel[0]) {
        cairo_set_font_size(cr, 12.0);
        cairo_text_extents(cr, cfg->ylabel, &te);
        cairo_save(cr);
        cairo_translate(cr, 18, (top + bottom) / 2);
        cairo_rotate(cr, -M_PI / 2);
        cairo_move_to(cr, -te.width / 2 - te.x_bearing, 0);
        cairo_show_text(cr, cfg->ylabel);
        cairo_restore(cr);
    }
}

/* ---- Series draw routines ----
 * All take cfg so log_x / log_y propagate into coord mapping. */

static void draw_line(cairo_t *cr, const K26PSeries *s,
                      const K26PlotConfig *cfg,
                      double xmin, double xmax, double ymin, double ymax,
                      double left, double right, double top, double bottom)
{
    if (s->n < 1 || !s->xs || !s->ys) return;
    cairo_set_line_width(cr, s->linewidth > 0 ? s->linewidth : 1.5);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
    cairo_set_line_cap (cr, CAIRO_LINE_CAP_ROUND);
    int started = 0;
    for (size_t i = 0; i < s->n; i++) {
        double tx = scale_axis(s->xs[i], cfg->log_x);
        double ty = scale_axis(s->ys[i], cfg->log_y);
        if (isnan(tx) || isnan(ty)) { started = 0; continue; }
        double px = left + (tx - xmin) / (xmax - xmin) * (right - left);
        double py = bottom - (ty - ymin) / (ymax - ymin) * (bottom - top);
        if (!started) { cairo_move_to(cr, px, py); started = 1; }
        else          cairo_line_to(cr, px, py);
    }
    cairo_stroke(cr);
}

static void draw_scatter(cairo_t *cr, const K26PSeries *s,
                         const K26PlotConfig *cfg,
                         double xmin, double xmax, double ymin, double ymax,
                         double left, double right, double top, double bottom)
{
    if (!s->xs || !s->ys) return;
    double r = (s->point_size > 0 ? s->point_size : 4.0) * 0.5;
    for (size_t i = 0; i < s->n; i++) {
        double px = x_to_px(s->xs[i], cfg->log_x, xmin, xmax, left, right);
        double py = y_to_px(s->ys[i], cfg->log_y, ymin, ymax, top, bottom);
        cairo_arc(cr, px, py, r, 0, 2 * M_PI);
        cairo_fill(cr);
    }
}

static void draw_errorbar(cairo_t *cr, const K26PSeries *s,
                          const K26PlotConfig *cfg,
                          double xmin, double xmax, double ymin, double ymax,
                          double left, double right, double top, double bottom)
{
    if (!s->xs || !s->ys) return;
    double cap = 4.0;
    double dot = (s->point_size > 0 ? s->point_size : 4.0) * 0.5;
    cairo_set_line_width(cr, s->linewidth > 0 ? s->linewidth : 1.0);
    for (size_t i = 0; i < s->n; i++) {
        double px = x_to_px(s->xs[i], cfg->log_x, xmin, xmax, left, right);
        double py = y_to_px(s->ys[i], cfg->log_y, ymin, ymax, top, bottom);
        double e  = s->yerr ? s->yerr[i] : 0.0;
        double py_lo = y_to_px(s->ys[i] - e, cfg->log_y, ymin, ymax, top, bottom);
        double py_hi = y_to_px(s->ys[i] + e, cfg->log_y, ymin, ymax, top, bottom);
        cairo_move_to(cr, px,        py_hi);
        cairo_line_to(cr, px,        py_lo);
        cairo_move_to(cr, px - cap,  py_hi);
        cairo_line_to(cr, px + cap,  py_hi);
        cairo_move_to(cr, px - cap,  py_lo);
        cairo_line_to(cr, px + cap,  py_lo);
        cairo_stroke(cr);
        cairo_arc(cr, px, py, dot, 0, 2 * M_PI);
        cairo_fill(cr);
    }
}

/* Pre-binned: xs has n+1 edges, ys has n bin counts. */
static void draw_histogram(cairo_t *cr, const K26PSeries *s,
                           const K26PlotConfig *cfg,
                           double xmin, double xmax, double ymin, double ymax,
                           double left, double right, double top, double bottom)
{
    if (!s->xs || !s->ys) return;
    double y0_px = y_to_px(0.0, cfg->log_y, ymin, ymax, top, bottom);
    for (size_t i = 0; i < s->n; i++) {
        double pxa = x_to_px(s->xs[i],     cfg->log_x, xmin, xmax, left, right);
        double pxb = x_to_px(s->xs[i + 1], cfg->log_x, xmin, xmax, left, right);
        double pyb = y_to_px(s->ys[i],     cfg->log_y, ymin, ymax, top, bottom);
        double rect_h = y0_px - pyb;
        cairo_rectangle(cr, pxa, pyb, pxb - pxa, rect_h);
        cairo_fill_preserve(cr);
        cairo_save(cr);
        cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 0.5);
        cairo_set_line_width(cr, 0.5);
        cairo_stroke(cr);
        cairo_restore(cr);
    }
}

/* Box plot. Caller pre-computes Q1 / median / Q3 / whisker extents.
 * Body fills with translucent series color; outline + median line +
 * whiskers + outlier dots use the solid series color. */
static void draw_box(cairo_t *cr, const K26PSeries *s,
                     const K26PlotConfig *cfg,
                     double pr, double pg, double pb,
                     double xmin, double xmax, double ymin, double ymax,
                     double left, double right, double top, double bottom)
{
    if (s->n == 0 || !s->xs || !s->box_q1 || !s->box_median ||
        !s->box_q3 || !s->box_whisker_lo || !s->box_whisker_hi) return;

    const double box_half = 12.0;          /* fixed pixel half-width */
    const double cap_half = 6.0;
    cairo_set_line_width(cr, s->linewidth > 0 ? s->linewidth : 1.5);

    for (size_t i = 0; i < s->n; i++) {
        double px    = x_to_px(s->xs[i],            cfg->log_x, xmin, xmax, left, right);
        double py_q1 = y_to_px(s->box_q1[i],        cfg->log_y, ymin, ymax, top, bottom);
        double py_md = y_to_px(s->box_median[i],    cfg->log_y, ymin, ymax, top, bottom);
        double py_q3 = y_to_px(s->box_q3[i],        cfg->log_y, ymin, ymax, top, bottom);
        double py_wl = y_to_px(s->box_whisker_lo[i],cfg->log_y, ymin, ymax, top, bottom);
        double py_wh = y_to_px(s->box_whisker_hi[i],cfg->log_y, ymin, ymax, top, bottom);
        double bx0 = px - box_half;
        double bx1 = px + box_half;
        /* Cairo y grows down; q3 is higher data → smaller py. */
        double rect_y = (py_q3 < py_q1) ? py_q3 : py_q1;
        double rect_h = fabs(py_q1 - py_q3);

        /* Body fill. */
        cairo_save(cr);
        cairo_set_source_rgba(cr, pr, pg, pb, 0.28);
        cairo_rectangle(cr, bx0, rect_y, bx1 - bx0, rect_h);
        cairo_fill(cr);
        cairo_restore(cr);

        /* Body outline. */
        cairo_set_source_rgb(cr, pr, pg, pb);
        cairo_rectangle(cr, bx0, rect_y, bx1 - bx0, rect_h);
        cairo_stroke(cr);

        /* Median line (heavier). */
        cairo_save(cr);
        cairo_set_line_width(cr, (s->linewidth > 0 ? s->linewidth : 1.5) + 1.0);
        cairo_move_to(cr, bx0, py_md);
        cairo_line_to(cr, bx1, py_md);
        cairo_stroke(cr);
        cairo_restore(cr);

        /* Whisker stems. */
        cairo_move_to(cr, px, py_q1); cairo_line_to(cr, px, py_wl);
        cairo_move_to(cr, px, py_q3); cairo_line_to(cr, px, py_wh);
        /* Whisker caps. */
        cairo_move_to(cr, px - cap_half, py_wl);
        cairo_line_to(cr, px + cap_half, py_wl);
        cairo_move_to(cr, px - cap_half, py_wh);
        cairo_line_to(cr, px + cap_half, py_wh);
        cairo_stroke(cr);
    }

    /* Outlier dots — flat (x, y) arrays across all boxes. */
    if (s->outliers_x && s->outliers_y && s->n_outliers > 0) {
        double dot = (s->point_size > 0 ? s->point_size : 4.0) * 0.5;
        cairo_set_source_rgb(cr, pr, pg, pb);
        for (size_t i = 0; i < s->n_outliers; i++) {
            double px = x_to_px(s->outliers_x[i], cfg->log_x, xmin, xmax, left, right);
            double py = y_to_px(s->outliers_y[i], cfg->log_y, ymin, ymax, top, bottom);
            cairo_arc(cr, px, py, dot, 0, 2 * M_PI);
            cairo_fill(cr);
        }
    }
}

/* Map normalised t in [0,1] through the heatmap colour gradient. */
static void heat_color(double t, double *r, double *g, double *b)
{
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    int i = 0;
    while (i < 3 && t > HEAT_STOPS[i + 1][0]) i++;
    double t0 = HEAT_STOPS[i][0];
    double t1 = HEAT_STOPS[i + 1][0];
    double f = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
    *r = HEAT_STOPS[i][1] + f * (HEAT_STOPS[i + 1][1] - HEAT_STOPS[i][1]);
    *g = HEAT_STOPS[i][2] + f * (HEAT_STOPS[i + 1][2] - HEAT_STOPS[i][2]);
    *b = HEAT_STOPS[i][3] + f * (HEAT_STOPS[i + 1][3] - HEAT_STOPS[i][3]);
}

/* Heatmap. Renders heat_data as a grid of coloured cells in
 * [heat_x0, heat_x1] × [heat_y0, heat_y1] axis space. NaN cells are
 * skipped (transparent). */
static void draw_heatmap(cairo_t *cr, const K26PSeries *s,
                         const K26PlotConfig *cfg,
                         double xmin, double xmax, double ymin, double ymax,
                         double left, double right, double top, double bottom)
{
    if (!s->heat_data || s->heat_rows == 0 || s->heat_cols == 0) return;

    /* Resolve color-scale bounds (NaN-fill from data). */
    double vmin = s->heat_vmin, vmax = s->heat_vmax;
    if (isnan(vmin) || isnan(vmax)) {
        double mn =  INFINITY, mx = -INFINITY;
        size_t total = s->heat_rows * s->heat_cols;
        for (size_t i = 0; i < total; i++) {
            double v = s->heat_data[i];
            if (isnan(v)) continue;
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        if (isnan(vmin)) vmin = isfinite(mn) ? mn : 0.0;
        if (isnan(vmax)) vmax = isfinite(mx) ? mx : 1.0;
    }
    if (vmax <= vmin) vmax = vmin + 1.0;

    double cell_w = (s->heat_x1 - s->heat_x0) / (double)s->heat_cols;
    double cell_h = (s->heat_y1 - s->heat_y0) / (double)s->heat_rows;

    for (size_t r = 0; r < s->heat_rows; r++) {
        for (size_t c = 0; c < s->heat_cols; c++) {
            double v = s->heat_data[r * s->heat_cols + c];
            if (isnan(v)) continue;
            double t = (v - vmin) / (vmax - vmin);
            double cr0, cg0, cb0;
            heat_color(t, &cr0, &cg0, &cb0);

            double dx0 = s->heat_x0 + (double)c       * cell_w;
            double dx1 = s->heat_x0 + (double)(c + 1) * cell_w;
            double dy0 = s->heat_y0 + (double)r       * cell_h;
            double dy1 = s->heat_y0 + (double)(r + 1) * cell_h;

            double pxa = x_to_px(dx0, cfg->log_x, xmin, xmax, left, right);
            double pxb = x_to_px(dx1, cfg->log_x, xmin, xmax, left, right);
            double pya = y_to_px(dy0, cfg->log_y, ymin, ymax, top, bottom);
            double pyb = y_to_px(dy1, cfg->log_y, ymin, ymax, top, bottom);
            double y_top = (pya < pyb) ? pya : pyb;
            double cell_h_px = fabs(pyb - pya);

            cairo_set_source_rgb(cr, cr0, cg0, cb0);
            /* +0.5px overdraw on each edge so adjacent cells don't show
             * 1-pixel seams from rasterisation rounding. */
            cairo_rectangle(cr, pxa - 0.5, y_top - 0.5,
                                (pxb - pxa) + 1.0, cell_h_px + 1.0);
            cairo_fill(cr);
        }
    }
}

/* ---- Legend ---- */
static void draw_legend(cairo_t *cr, double w,
                        const K26PSeries *series, size_t n_series,
                        double top)
{
    size_t rows = 0;
    for (size_t i = 0; i < n_series; i++)
        if (series[i].label && series[i].label[0]) rows++;
    if (rows == 0) return;

    cairo_select_font_face(cr, "sans-serif",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 11.0);

    double max_label_w = 0;
    cairo_text_extents_t te;
    for (size_t i = 0; i < n_series; i++) {
        const K26PSeries *s = &series[i];
        if (!s->label || !s->label[0]) continue;
        cairo_text_extents(cr, s->label, &te);
        if (te.width > max_label_w) max_label_w = te.width;
    }
    double row_h = 16.0;
    double box_w = LEGEND_PAD * 2 + LEGEND_LINE_LEN + 6 + max_label_w;
    double box_h = LEGEND_PAD * 2 + row_h * (double)rows;
    double box_x = w - MARGIN_RIGHT - box_w - 6;
    double box_y = top + 6;

    cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 0.85);
    cairo_rectangle(cr, box_x, box_y, box_w, box_h);
    cairo_fill_preserve(cr);
    cairo_set_source_rgba(cr, 0.4, 0.4, 0.4, 0.7);
    cairo_set_line_width(cr, 0.5);
    cairo_stroke(cr);

    double row_y = box_y + LEGEND_PAD + row_h * 0.5;
    size_t rendered = 0;
    for (size_t i = 0; i < n_series; i++) {
        const K26PSeries *s = &series[i];
        if (!s->label || !s->label[0]) continue;
        double pr, pg, pb;
        series_color(s, i, &pr, &pg, &pb);
        cairo_set_source_rgb(cr, pr, pg, pb);
        cairo_set_line_width(cr, 2.0);

        double x0 = box_x + LEGEND_PAD;
        double x1 = x0 + LEGEND_LINE_LEN;
        if (s->kind == K26P_HISTOGRAM || s->kind == K26P_BOX) {
            cairo_rectangle(cr, x0, row_y - 5, LEGEND_LINE_LEN, 10);
            cairo_fill(cr);
        } else if (s->kind == K26P_HEATMAP) {
            /* Mini-gradient as the swatch — three cells, low/mid/high. */
            for (int j = 0; j < 3; j++) {
                double t = (double)j / 2.0;
                double cr0, cg0, cb0;
                heat_color(t, &cr0, &cg0, &cb0);
                cairo_set_source_rgb(cr, cr0, cg0, cb0);
                cairo_rectangle(cr, x0 + j * (LEGEND_LINE_LEN / 3.0),
                                    row_y - 5,
                                    LEGEND_LINE_LEN / 3.0, 10);
                cairo_fill(cr);
            }
        } else {
            cairo_move_to(cr, x0, row_y);
            cairo_line_to(cr, x1, row_y);
            cairo_stroke(cr);
            if (s->kind == K26P_SCATTER || s->kind == K26P_ERRORBAR) {
                double dot = (s->point_size > 0 ? s->point_size : 4.0) * 0.5;
                cairo_arc(cr, (x0 + x1) / 2, row_y, dot, 0, 2 * M_PI);
                cairo_fill(cr);
            }
        }
        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_move_to(cr, x1 + 6, row_y + 4);
        cairo_show_text(cr, s->label);

        row_y += row_h;
        rendered++;
        if (rendered >= rows) break;
    }
}

/* ---- Master draw routine ---- */
static void render_to_cairo(cairo_t *cr, double w, double h,
                            const K26PlotConfig *cfg,
                            const K26PSeries *series, size_t n_series)
{
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_rectangle(cr, 0, 0, w, h);
    cairo_fill(cr);

    double xmin, xmax, ymin, ymax;
    resolve_extents(cfg, series, n_series, &xmin, &xmax, &ymin, &ymax);

    draw_axes(cr, w, h, cfg, xmin, xmax, ymin, ymax);

    double left   = MARGIN_LEFT;
    double right  = w - MARGIN_RIGHT;
    double top    = MARGIN_TOP;
    double bottom = h - MARGIN_BOTTOM;

    cairo_save(cr);
    cairo_rectangle(cr, left, top, right - left, bottom - top);
    cairo_clip(cr);

    for (size_t i = 0; i < n_series; i++) {
        const K26PSeries *s = &series[i];
        double pr, pg, pb;
        series_color(s, i, &pr, &pg, &pb);
        cairo_set_source_rgb(cr, pr, pg, pb);
        switch (s->kind) {
        case K26P_LINE:
            draw_line(cr, s, cfg, xmin, xmax, ymin, ymax, left, right, top, bottom);
            break;
        case K26P_SCATTER:
            draw_scatter(cr, s, cfg, xmin, xmax, ymin, ymax, left, right, top, bottom);
            break;
        case K26P_ERRORBAR:
            draw_errorbar(cr, s, cfg, xmin, xmax, ymin, ymax, left, right, top, bottom);
            break;
        case K26P_HISTOGRAM:
            draw_histogram(cr, s, cfg, xmin, xmax, ymin, ymax, left, right, top, bottom);
            break;
        case K26P_BOX:
            draw_box(cr, s, cfg, pr, pg, pb,
                     xmin, xmax, ymin, ymax, left, right, top, bottom);
            break;
        case K26P_HEATMAP:
            draw_heatmap(cr, s, cfg, xmin, xmax, ymin, ymax,
                         left, right, top, bottom);
            break;
        }
    }
    cairo_restore(cr);

    if (cfg->show_legend) {
        draw_legend(cr, w, series, n_series, top);
    }
}

/* ---- Public live-render entry: caller-provided Cairo context ----
 *
 * Renders the plot at (dst_x, dst_y) with dimensions (dst_w, dst_h)
 * onto the caller's cairo_t. The caller is responsible for the
 * surface; this routine clips, translates, draws, and restores.
 *
 * Designed for live in-place rendering into a caller surface: open a
 * cairo_xlib_surface against fl_window, call this, destroy. The
 * widget appears with the plot rendered in-place. */
K26PStatus k26plot_draw_to_cairo(cairo_t *cr,
                                  int dst_x, int dst_y,
                                  int dst_w, int dst_h,
                                  const K26PlotConfig *cfg,
                                  const K26PSeries *series, size_t n_series)
{
    if (!cr || !cfg) return K26P_ERR_INVAL;
    if (n_series > 0 && !series) return K26P_ERR_INVAL;
    if (dst_w <= 0 || dst_h <= 0) return K26P_ERR_INVAL;

    cairo_save(cr);
    cairo_translate(cr, dst_x, dst_y);
    cairo_rectangle(cr, 0, 0, dst_w, dst_h);
    cairo_clip(cr);
    render_to_cairo(cr, (double)dst_w, (double)dst_h,
                    cfg, series, n_series);
    cairo_restore(cr);
    return K26P_OK;
}

/* ---- Public render entry: PNG and/or SVG ---- */
K26PStatus k26plot_render(const K26PlotConfig *cfg,
                          const K26PSeries *series, size_t n_series,
                          const char *png_path,
                          const char *svg_path)
{
    if (!cfg || cfg->width_px <= 0 || cfg->height_px <= 0) return K26P_ERR_INVAL;
    if (n_series > 0 && !series) return K26P_ERR_INVAL;
    if (!png_path && !svg_path)  return K26P_ERR_INVAL;

    double w = (double)cfg->width_px;
    double h = (double)cfg->height_px;

    if (png_path) {
        cairo_surface_t *surf = cairo_image_surface_create(
            CAIRO_FORMAT_ARGB32, cfg->width_px, cfg->height_px);
        if (!surf || cairo_surface_status(surf) != CAIRO_STATUS_SUCCESS) {
            if (surf) cairo_surface_destroy(surf);
            return K26P_ERR_OOM;
        }
        cairo_t *cr = cairo_create(surf);
        if (!cr || cairo_status(cr) != CAIRO_STATUS_SUCCESS) {
            if (cr) cairo_destroy(cr);
            cairo_surface_destroy(surf);
            return K26P_ERR_OOM;
        }
        render_to_cairo(cr, w, h, cfg, series, n_series);
        cairo_status_t wrc = cairo_surface_write_to_png(surf, png_path);
        cairo_destroy(cr);
        cairo_surface_destroy(surf);
        if (wrc != CAIRO_STATUS_SUCCESS) return K26P_ERR_IO;
    }

    if (svg_path) {
        cairo_surface_t *surf = cairo_svg_surface_create(svg_path, w, h);
        if (!surf || cairo_surface_status(surf) != CAIRO_STATUS_SUCCESS) {
            if (surf) cairo_surface_destroy(surf);
            return K26P_ERR_IO;
        }
        cairo_t *cr = cairo_create(surf);
        if (!cr || cairo_status(cr) != CAIRO_STATUS_SUCCESS) {
            if (cr) cairo_destroy(cr);
            cairo_surface_destroy(surf);
            return K26P_ERR_OOM;
        }
        render_to_cairo(cr, w, h, cfg, series, n_series);
        cairo_destroy(cr);
        cairo_surface_finish(surf);
        cairo_surface_destroy(surf);
    }

    return K26P_OK;
}
