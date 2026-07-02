/* k26plot-test — standalone smoke test for libk26plot.
 *
 * Produces three reference figures into a caller-supplied output dir:
 *   1. fit.png/svg       — synthetic Hooke's-law data (scatter) + OLS
 *                          fit line + error bars; uses libk26compute
 *                          for the regression.
 *   2. histogram.png/svg — histogram of normal-distribution samples.
 *   3. multi.png/svg     — line plot of three trig functions.
 *
 * Each figure exercises a different draw path. Eyeball test: open the
 * PNGs and verify axes, grid, ticks, legend, and series rendering all
 * read as publication-grade. Exit 0 if every render returns OK. */

#include "k26plot.h"
#include "k26compute.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int run_fit(const char *outdir);
static int run_hist(const char *outdir);
static int run_multi(const char *outdir);

/* IMPORTANT: do NOT reuse a static buffer across the two argument
 * evaluations of a single render call — both args would resolve to
 * the same address with whatever content the second call wrote,
 * causing the PNG path to silently become the SVG path. Both render
 * outputs would then collide on one file. Use per-call stack
 * buffers instead. */
#define MAKE_PATH(buf, outdir, name) \
    snprintf((buf), sizeof(buf), "%s/%s", (outdir), (name))

int main(int argc, char **argv)
{
    const char *outdir = (argc >= 2) ? argv[1] : "/tmp";
    int fail = 0;
    fail += run_fit(outdir);
    fail += run_hist(outdir);
    fail += run_multi(outdir);
    if (fail == 0) {
        printf("\nlibk26plot v%s smoke OK — outputs in %s/\n",
               K26PLOT_LIB_VERSION, outdir);
    } else {
        printf("\n%d figure(s) failed\n", fail);
    }
    return fail ? 1 : 0;
}

/* ---- 1. Fit + scatter + errorbars ---- */
static int run_fit(const char *outdir)
{
    /* Synthetic Hooke's-law data with deterministic noise. */
    const size_t N = 12;
    double xs[12], ys[12], yerr[12];
    K26CRng rng; k26c_rng_init(&rng, 0xC0FFEE);
    for (size_t i = 0; i < N; i++) {
        xs[i] = (double)i;
        double noise = k26c_rng_normal(&rng, 0.0, 0.3);
        ys[i]   = 2.5 * xs[i] + noise;
        yerr[i] = 0.3;
    }

    /* OLS fit through libk26compute. */
    K26CLinReg fit;
    if (k26c_stats_linreg(&fit, xs, ys, N) != K26C_OK) {
        fprintf(stderr, "linreg failed\n");
        return 1;
    }
    /* Fit-line samples for the line series. */
    double fxs[2] = { xs[0] - 0.5, xs[N-1] + 0.5 };
    double fys[2] = { fit.a + fit.b * fxs[0], fit.a + fit.b * fxs[1] };

    char title[160];
    snprintf(title, sizeof title,
             "Hooke's law fit:  y = %.3f + %.3f x   (r² = %.4f)",
             fit.a, fit.b, fit.r2);

    K26PlotConfig cfg;
    k26plot_config_defaults(&cfg);
    cfg.width_px  = 900;
    cfg.height_px = 600;
    cfg.title  = title;
    cfg.xlabel = "displacement x";
    cfg.ylabel = "force F";

    K26PSeries series[2];
    memset(series, 0, sizeof series);
    series[0].kind = K26P_ERRORBAR;
    series[0].label = "measured";
    series[0].xs = xs; series[0].ys = ys; series[0].yerr = yerr;
    series[0].n = N;
    series[1].kind = K26P_LINE;
    series[1].label = "OLS fit";
    series[1].xs = fxs; series[1].ys = fys; series[1].n = 2;
    series[1].linewidth = 2.0;

    char png_path[1024], svg_path[1024];
    MAKE_PATH(png_path, outdir, "fit.png");
    MAKE_PATH(svg_path, outdir, "fit.svg");
    K26PStatus rc = k26plot_render(&cfg, series, 2, png_path, svg_path);
    if (rc != K26P_OK) {
        fprintf(stderr, "fit render: %s\n", k26plot_status_str(rc));
        return 1;
    }
    printf("OK fit         → %s/fit.{png,svg}\n", outdir);
    return 0;
}

/* ---- 2. Histogram of normal samples ---- */
static int run_hist(const char *outdir)
{
    /* Sample, bin, plot. Deterministic. */
    const size_t N = 5000;
    K26CRng rng; k26c_rng_init(&rng, 0xBADCAFE);
    double samples[5000];
    for (size_t i = 0; i < N; i++) samples[i] = k26c_rng_normal(&rng, 0.0, 1.0);

    /* 30 equal-width bins from -4 to +4. */
    enum { NBINS = 30 };
    double edges[NBINS + 1];
    double counts[NBINS] = { 0 };
    double lo = -4.0, hi = 4.0;
    double step = (hi - lo) / NBINS;
    for (int i = 0; i <= NBINS; i++) edges[i] = lo + i * step;
    for (size_t i = 0; i < N; i++) {
        int b = (int)floor((samples[i] - lo) / step);
        if (b < 0 || b >= NBINS) continue;
        counts[b] += 1.0;
    }

    K26PlotConfig cfg;
    k26plot_config_defaults(&cfg);
    cfg.width_px  = 900;
    cfg.height_px = 600;
    cfg.title  = "Normal(0, 1)  — 5000 samples";
    cfg.xlabel = "value";
    cfg.ylabel = "count";

    K26PSeries series[1];
    memset(series, 0, sizeof series);
    series[0].kind = K26P_HISTOGRAM;
    series[0].label = "samples";
    series[0].xs = edges;        /* NBINS+1 entries */
    series[0].ys = counts;       /* NBINS entries */
    series[0].n  = NBINS;

    char png_path[1024], svg_path[1024];
    MAKE_PATH(png_path, outdir, "histogram.png");
    MAKE_PATH(svg_path, outdir, "histogram.svg");
    K26PStatus rc = k26plot_render(&cfg, series, 1, png_path, svg_path);
    if (rc != K26P_OK) {
        fprintf(stderr, "histogram render: %s\n", k26plot_status_str(rc));
        return 1;
    }
    printf("OK histogram   → %s/histogram.{png,svg}\n", outdir);
    return 0;
}

/* ---- 3. Multi-series line plot ---- */
static int run_multi(const char *outdir)
{
    enum { N = 200 };
    static double xs[N], y_sin[N], y_cos[N], y_decay[N];
    for (int i = 0; i < N; i++) {
        double t = (double)i * (4.0 * M_PI) / (N - 1);
        xs[i]     = t;
        y_sin[i]  = sin(t);
        y_cos[i]  = cos(t);
        y_decay[i] = exp(-t * 0.3) * sin(t);
    }

    K26PlotConfig cfg;
    k26plot_config_defaults(&cfg);
    cfg.width_px  = 900;
    cfg.height_px = 600;
    cfg.title  = "Three references";
    cfg.xlabel = "t";
    cfg.ylabel = "amplitude";

    K26PSeries series[3];
    memset(series, 0, sizeof series);
    series[0].kind = K26P_LINE; series[0].label = "sin(t)";
    series[0].xs = xs; series[0].ys = y_sin;   series[0].n = N;
    series[1].kind = K26P_LINE; series[1].label = "cos(t)";
    series[1].xs = xs; series[1].ys = y_cos;   series[1].n = N;
    series[2].kind = K26P_LINE; series[2].label = "exp(-0.3 t) sin(t)";
    series[2].xs = xs; series[2].ys = y_decay; series[2].n = N;

    char png_path[1024], svg_path[1024];
    MAKE_PATH(png_path, outdir, "multi.png");
    MAKE_PATH(svg_path, outdir, "multi.svg");
    K26PStatus rc = k26plot_render(&cfg, series, 3, png_path, svg_path);
    if (rc != K26P_OK) {
        fprintf(stderr, "multi render: %s\n", k26plot_status_str(rc));
        return 1;
    }
    printf("OK multi-line  → %s/multi.{png,svg}\n", outdir);
    return 0;
}
