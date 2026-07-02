/* star_catalogue_hip.c — embedded Hipparcos-2 (van Leeuwen 2007)
 * baked-in naked-eye subset.
 *
 * This file ships a fixed catalogue of ~60 of the brightest /
 * most-recognisable stars so that K26 has a usable default
 * starfield without requiring the external apk add-ons. When
 * `/usr/share/k26astro/stars/hip_naked.k26stars` is installed
 * (the ISO default), the file-backed catalogue supersedes this
 * subset.
 *
 * Source: Hipparcos-2 (van Leeuwen 2007 New Reduction), epoch
 * J1991.25 TT. RA / Dec in degrees, converted to radians at
 * compile time. Parallax in mas, proper motion in mas/yr.
 *
 * Why this subset: the magnitude V ≤ 2.0 stars (the brightest
 * ~30) define the dominant naked-eye sky from any solar-system
 * observer; the additional ~30 below extend to the V ≤ 3 set
 * including notable navigation stars (Polaris, Mizar, etc.).
 *
 * Updates: when a new catalogue release lands (Gaia DR4 = 2026-12-02),
 * the file-backed `.k26stars` loader picks it up via the apk;
 * this baked-in table doesn't need to change. To replace this
 * embedded set wholesale, regenerate the .c file from a fresh
 * Hipparcos query (or, in a Gaia DR4-flavoured K26 release, port
 * the entries to DR4 source ids + epoch J2017.5).
 */
#include "k26astro_body/star.h"
#include "k26astro_core/consts.h"
#include "internal.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Convenience: degrees → radians at file scope. */
#define DEG (K26A_RAD_PER_DEG)
/* HMS → degrees converter for RA, used inline below. */
#define HMS(h, m, s) (((h) + (m) / 60.0 + (s) / 3600.0) * 15.0)
/* DMS → degrees converter for Dec (positive or negative). */
#define DMS(d, m, s) ((d) + ((d) < 0 ? -1.0 : 1.0) * ((m) / 60.0 + (s) / 3600.0))

/* Compile-time RA in degrees → radians. */
#define RA_DEG(d)  ((d) * DEG)
#define DEC_DEG(d) ((d) * DEG)

/* J1991.25 = JD 2448349.0625 (TT). */
#define HIP_EPOCH_JD_TT 2448349.0625

/* The baked-in star records. Each entry's RA/Dec are at J1991.25;
 * apparent-direction queries at simulation epoch propagate via
 * proper motion in star_catalogue.c:k26astro_star_apparent_dir. */
static const K26AstroStar HIP_STARS[] = {
    /* ---- magnitude ≤ 0 — the standout stars ------------------ */
    /* Sirius (α CMa) — HIP 32349 */
    { 32349, "HIP 32349",
      RA_DEG(101.28715533), DEC_DEG(-16.71611586),
      379.21, -546.05, -1223.14, -7.6,
      -1.46, 0.01 },
    /* Canopus (α Car) — HIP 30438 */
    { 30438, "HIP 30438",
      RA_DEG( 95.98787790), DEC_DEG(-52.69566138),
       10.43,   19.99,   23.67, 20.5,
      -0.72, 0.16 },
    /* Arcturus (α Boo) — HIP 69673 */
    { 69673, "HIP 69673",
      RA_DEG(213.91528842), DEC_DEG( 19.18240916),
       88.83, -1093.45, -1999.40, -5.2,
      -0.05, 1.23 },
    /* Vega (α Lyr) — HIP 91262 */
    { 91262, "HIP 91262",
      RA_DEG(279.23473479), DEC_DEG( 38.78368896),
      130.23,  200.94,  286.23, -13.9,
       0.03, 0.00 },
    /* Capella (α Aur) — HIP 24608 */
    { 24608, "HIP 24608",
      RA_DEG( 79.17232862), DEC_DEG( 45.99799100),
       77.29,   75.52, -427.13, 29.4,
       0.08, 0.80 },
    /* Rigel (β Ori) — HIP 24436 */
    { 24436, "HIP 24436",
      RA_DEG( 78.63446707), DEC_DEG( -8.20163836),
        3.78,    1.31,   0.50, 17.8,
       0.13, -0.03 },
    /* Procyon (α CMi) — HIP 37279 */
    { 37279, "HIP 37279",
      RA_DEG(114.82549322), DEC_DEG(  5.22498756),
      284.56, -716.57, -1034.58, -3.2,
       0.34, 0.42 },
    /* Achernar (α Eri) — HIP 7588 */
    {  7588, "HIP 7588",
      RA_DEG( 24.42852898), DEC_DEG(-57.23675792),
       23.39,   88.02, -40.08, 16.0,
       0.46, -0.16 },

    /* ---- magnitude 0.0–1.0 ------------------------------- */
    /* Betelgeuse (α Ori) — HIP 27989 */
    { 27989, "HIP 27989",
      RA_DEG( 88.79293875), DEC_DEG(  7.40706418),
        6.55,   24.95,    9.56, 21.0,
       0.50, 1.85 },
    /* Hadar (β Cen) — HIP 68702 */
    { 68702, "HIP 68702",
      RA_DEG(210.95583970), DEC_DEG(-60.37303959),
        8.32,   -33.96, -25.06, 5.9,
       0.61, -0.23 },
    /* Altair (α Aql) — HIP 97649 */
    { 97649, "HIP 97649",
      RA_DEG(297.69582730), DEC_DEG(  8.86832120),
      194.95,  536.23,  385.29, -26.1,
       0.77, 0.22 },
    /* Aldebaran (α Tau) — HIP 21421 */
    { 21421, "HIP 21421",
      RA_DEG( 68.98016279), DEC_DEG( 16.50930208),
       48.94,   62.78,  -189.36, 54.3,
       0.86, 1.54 },
    /* Acrux (α Cru) — HIP 60718 */
    { 60718, "HIP 60718",
      RA_DEG(186.64965749), DEC_DEG(-63.09909340),
       10.13,  -35.83,  -14.86,  -11.2,
       0.77, -0.24 },
    /* Antares (α Sco) — HIP 80763 */
    { 80763, "HIP 80763",
      RA_DEG(247.35196792), DEC_DEG(-26.43200447),
        5.40,   -10.16, -23.21, -3.2,
       0.91, 1.83 },
    /* Spica (α Vir) — HIP 65474 */
    { 65474, "HIP 65474",
      RA_DEG(201.29824733), DEC_DEG(-11.16131803),
       13.06,   -42.50, -31.73,  1.0,
       0.98, -0.24 },
    /* Pollux (β Gem) — HIP 37826 */
    { 37826, "HIP 37826",
      RA_DEG(116.32896464), DEC_DEG( 28.02609710),
       96.54,  -625.69, -45.95, 3.2,
       1.14, 0.99 },
    /* Fomalhaut (α PsA) — HIP 113368 */
    {113368, "HIP 113368",
      RA_DEG(344.41269282), DEC_DEG(-29.62223812),
      129.81,  329.22, -164.22, 6.5,
       1.16, 0.09 },
    /* Deneb (α Cyg) — HIP 102098 */
    {102098, "HIP 102098",
      RA_DEG(310.35797800), DEC_DEG( 45.28033827),
        2.31,    1.99,    1.85, -4.9,
       1.25, 0.09 },
    /* Mimosa (β Cru) — HIP 62434 */
    { 62434, "HIP 62434",
      RA_DEG(191.93028841), DEC_DEG(-59.68873242),
        9.25,   -42.97, -16.18, 15.6,
       1.25, -0.24 },
    /* Regulus (α Leo) — HIP 49669 */
    { 49669, "HIP 49669",
      RA_DEG(152.09296233), DEC_DEG( 11.96720889),
       41.13,   -249.40,    4.91, 5.9,
       1.36, -0.11 },
    /* Adhara (ε CMa) — HIP 33579 */
    { 33579, "HIP 33579",
      RA_DEG(104.65645141), DEC_DEG(-28.97208866),
        7.57,    2.63,    2.29, 27.3,
       1.50, -0.21 },
    /* Castor (α Gem) — HIP 36850 */
    { 36850, "HIP 36850",
      RA_DEG(113.64950419), DEC_DEG( 31.88827661),
       64.12,  -191.45,  -145.20, 5.4,
       1.58, 0.04 },
    /* Gacrux (γ Cru) — HIP 61084 */
    { 61084, "HIP 61084",
      RA_DEG(187.79150467), DEC_DEG(-57.11317290),
       36.83,    27.94, -264.33, 21.4,
       1.59, 1.59 },
    /* Shaula (λ Sco) — HIP 85927 */
    { 85927, "HIP 85927",
      RA_DEG(263.40216853), DEC_DEG(-37.10374756),
        5.71,   -8.90,  -29.95, -3.0,
       1.62, -0.22 },

    /* ---- navigation + recognition stars ------------------------ */
    /* Polaris (α UMi) — HIP 11767 */
    { 11767, "HIP 11767",
      RA_DEG( 37.95456067), DEC_DEG( 89.26410897),
        7.54,   44.22,  -11.74, -16.4,
       1.97, 0.60 },
    /* Algol (β Per) — HIP 14576 */
    { 14576, "HIP 14576",
      RA_DEG( 47.04221820), DEC_DEG( 40.95564838),
       35.14,    2.39,  -1.44, 4.0,
       2.09, -0.05 },
    /* Mizar (ζ UMa) — HIP 65378 */
    { 65378, "HIP 65378",
      RA_DEG(200.98142250), DEC_DEG( 54.92541625),
       41.73,   119.01, -25.91, -6.3,
       2.27, 0.15 },
    /* Albireo (β Cyg) — HIP 95947 */
    { 95947, "HIP 95947",
      RA_DEG(292.68034250), DEC_DEG( 27.96168833),
        8.46,    -7.17, -5.79, -24.0,
       3.05, 1.13 },
    /* Bellatrix (γ Ori) — HIP 25336 */
    { 25336, "HIP 25336",
      RA_DEG( 81.28276408), DEC_DEG(  6.34970223),
       13.42,    -8.75, -13.28, 18.2,
       1.64, -0.22 },
    /* Elnath (β Tau) — HIP 25428 */
    { 25428, "HIP 25428",
      RA_DEG( 81.57291278), DEC_DEG( 28.60745169),
       24.89,    23.28, -174.22, 9.2,
       1.65, -0.13 },
    /* Mintaka (δ Ori) — HIP 25930 */
    { 25930, "HIP 25930",
      RA_DEG( 83.00166759), DEC_DEG( -0.29909161),
        4.71,    1.67,    0.56, 16.0,
       2.23, -0.18 },
    /* Alnilam (ε Ori) — HIP 26311 */
    { 26311, "HIP 26311",
      RA_DEG( 84.05338889), DEC_DEG( -1.20191947),
        1.65,    1.49,   -1.06, 25.9,
       1.69, -0.18 },
    /* Alnitak (ζ Ori) — HIP 26727 */
    { 26727, "HIP 26727",
      RA_DEG( 85.18968762), DEC_DEG( -1.94257119),
        4.43,    3.99,    2.54, 18.5,
       1.79, -0.21 },
    /* Saiph (κ Ori) — HIP 27366 */
    { 27366, "HIP 27366",
      RA_DEG( 86.93911679), DEC_DEG( -9.66961011),
        4.52,    1.55,   -1.20, 20.5,
       2.07, -0.17 },
    /* Rigil Kentaurus (α Cen A) — HIP 71683 */
    { 71683, "HIP 71683",
      RA_DEG(219.90205833), DEC_DEG(-60.83397222),
      747.23, -3678.19,  481.84, -22.4,
       0.00, 0.67 },
    /* Proxima Centauri — HIP 70890 */
    { 70890, "HIP 70890",
      RA_DEG(217.42895408), DEC_DEG(-62.67949933),
      771.64, -3775.40,  768.10, -22.4,
      11.13, 1.84 },
    /* Beta Centauri (Agena) — HIP 68702 already above */
    /* Hamal (α Ari) — HIP 9884 */
    {  9884, "HIP 9884",
      RA_DEG( 31.79336764), DEC_DEG( 23.46241608),
       49.56,    190.73, -145.77, -14.0,
       2.00, 1.15 },
    /* Diphda (β Cet) — HIP 3419 */
    {  3419, "HIP 3419",
      RA_DEG( 10.89738378), DEC_DEG(-17.98660442),
       33.86,    232.79,   31.99, 13.0,
       2.04, 1.02 },
    /* Schedar (α Cas) — HIP 3179 */
    {  3179, "HIP 3179",
      RA_DEG( 10.12683870), DEC_DEG( 56.53733114),
       14.27,    50.34, -32.40, -3.8,
       2.24, 1.17 },
    /* Caph (β Cas) — HIP 746 */
    {   746, "HIP 746",
      RA_DEG(  2.29452167), DEC_DEG( 59.14977889),
       59.89,    523.50, -179.81, 11.0,
       2.27, 0.34 },
    /* Mirfak (α Per) — HIP 15863 */
    { 15863, "HIP 15863",
      RA_DEG( 51.08070756), DEC_DEG( 49.86124194),
        6.44,    24.11, -26.21, -2.0,
       1.79, 0.48 },
    /* Cor Caroli (α CVn) — HIP 63125 */
    { 63125, "HIP 63125",
      RA_DEG(193.93443672), DEC_DEG( 38.31821194),
       29.60,    -239.40,  56.65, -3.0,
       2.89, -0.12 },
    /* Cor Hydrae (α Hya, Alphard) — HIP 46390 */
    { 46390, "HIP 46390",
      RA_DEG(141.89697053), DEC_DEG( -8.65866756),
       18.40,    -14.49,  33.49, -4.3,
       1.99, 1.44 },
    /* Etamin (γ Dra) — HIP 87833 */
    { 87833, "HIP 87833",
      RA_DEG(269.15146533), DEC_DEG( 51.48894989),
       21.14,   -8.52, -23.05, -28.0,
       2.23, 1.52 },
    /* Vindemiatrix (ε Vir) — HIP 63608 */
    { 63608, "HIP 63608",
      RA_DEG(195.54414125), DEC_DEG( 10.95914519),
       29.76,   -275.05,   19.96, -14.3,
       2.83, 0.94 },
    /* Sadr (γ Cyg) — HIP 100453 */
    {100453, "HIP 100453",
      RA_DEG(305.55709694), DEC_DEG( 40.25667808),
        1.78,    2.39,   -0.91, -7.5,
       2.20, 0.67 },
    /* Sadalmelik (α Aqr) — HIP 109074 */
    {109074, "HIP 109074",
      RA_DEG(331.44570958), DEC_DEG( -0.31983044),
        6.23,    18.25, -9.40, 8.0,
       2.95, 0.97 },
    /* Sadalsuud (β Aqr) — HIP 106278 */
    {106278, "HIP 106278",
      RA_DEG(322.88971008), DEC_DEG( -5.57117272),
        5.33,    18.77,   -8.21, 6.5,
       2.87, 0.83 },
    /* Markab (α Peg) — HIP 113963 */
    {113963, "HIP 113963",
      RA_DEG(346.19022547), DEC_DEG( 15.20531697),
       23.36,    61.10,  -42.56, -4.0,
       2.49, -0.04 },
    /* Scheat (β Peg) — HIP 113881 */
    {113881, "HIP 113881",
      RA_DEG(345.94358950), DEC_DEG( 28.08246917),
       16.37,    187.65,  137.61, 9.0,
       2.42, 1.67 },
    /* Algenib (γ Peg) — HIP 1067 */
    {  1067, "HIP 1067",
      RA_DEG(  3.30889758), DEC_DEG( 15.18361139),
        9.79,    4.20, -8.79, 4.1,
       2.83, -0.23 },
    /* Enif (ε Peg) — HIP 107315 */
    {107315, "HIP 107315",
      RA_DEG(326.04632175), DEC_DEG(  9.87501292),
        4.85,    27.31,   1.50, 5.1,
       2.40, 1.53 },
    /* Wezen (δ CMa) — HIP 34444 */
    { 34444, "HIP 34444",
      RA_DEG(107.09785275), DEC_DEG(-26.39320042),
        1.82,    -2.75,  3.33, 34.3,
       1.83, 0.67 },
    /* Mirzam (β CMa) — HIP 30324 */
    { 30324, "HIP 30324",
      RA_DEG( 95.67492933), DEC_DEG(-17.95591725),
        6.53,    -3.45,  -0.47, 33.7,
       1.98, -0.24 },
    /* Aludra (η CMa) — HIP 35904 */
    { 35904, "HIP 35904",
      RA_DEG(111.02373722), DEC_DEG(-29.30311111),
        1.02,    -2.39,   6.05, 41.1,
       2.45, -0.08 },
    /* Mirach (β And) — HIP 5447 */
    {  5447, "HIP 5447",
      RA_DEG( 17.43298047), DEC_DEG( 35.62057989),
       16.36,    175.59, -112.23, 3.0,
       2.07, 1.58 },
    /* Almach (γ And) — HIP 9640 */
    {  9640, "HIP 9640",
      RA_DEG( 30.97473878), DEC_DEG( 42.32970694),
        9.19,    43.08, -49.85, -12.0,
       2.18, 1.37 },
    /* Alkaid (η UMa) — HIP 67301 */
    { 67301, "HIP 67301",
      RA_DEG(206.88571378), DEC_DEG( 49.31334089),
       31.38,    -121.17, -15.61, -10.9,
       1.85, -0.10 },
    /* Dubhe (α UMa) — HIP 54061 */
    { 54061, "HIP 54061",
      RA_DEG(165.93196614), DEC_DEG( 61.75103339),
       26.38,    -134.11, -34.70, -9.0,
       1.81, 1.06 },
    /* Merak (β UMa) — HIP 53910 */
    { 53910, "HIP 53910",
      RA_DEG(165.46031783), DEC_DEG( 56.38234389),
       40.90,    81.66,    33.74, -12.0,
       2.34, 0.03 }
};
/* N_HIP_STARS computed inline in k26astro_catalogue_hip_default below. */

K26AstroCatalogue *k26astro_catalogue_hip_default(void)
{
    static struct K26AstroCatalogue impl = {
        .map          = NULL,
        .map_size     = 0,
        .records      = HIP_STARS,
        .n_records    = (size_t)(sizeof(HIP_STARS) / sizeof(HIP_STARS[0])),
        .epoch_jd_tt  = HIP_EPOCH_JD_TT,
        .catalogue_id = K26_CAT_ID_HIPPARCOS_2007,
        .file_backed  = 0
    };
    return &impl;
}
