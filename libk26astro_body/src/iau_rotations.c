/* iau_rotations.c — IAU 2018 rotation-element coefficient tables.
 *
 * Source: Archinal et al. (2018), "Report of the IAU Working Group
 * on Cartographic Coordinates and Rotational Elements: 2015",
 * Celestial Mechanics and Dynamical Astronomy 130:22.
 *
 * 2019 erratum (Archinal et al. 2019, CeMec 131:61) applied where
 * relevant — Phobos prime meridian + Figure 1/2 sign corrections.
 *
 * Updates: when the IAU Working Group publishes a new convention,
 * edit individual entries in place. Each body's correction terms
 * are commented with their term identifier (M1..M14 for Mars,
 * E1..E13 for Moon's E_n series, Ja..Je for Jupiter, etc.) keyed
 * to the original publication so the cross-reference is explicit.
 *
 * Argument-rate convention: all `arg_b` fields are in deg / Julian
 * century, regardless of whether the target variable's rate is
 * per-day (W) or per-century (α / δ). The evaluator (rotation_iau.c)
 * scales by T_c uniformly.
 *
 * Amplitudes are in degrees (deg), consistent with the linear
 * terms. (IAU papers often tabulate the corrections in degrees
 * even when they're sub-arcsecond — we follow the same convention
 * here so the table reads against the publication directly.) */
#include "k26astro_body/rotation_model.h"

/* ---- Argument shorthand ------------------------------------- *
 *
 * Many bodies share their argument tables (the linear functions of
 * T_c that drive sin/cos terms). The IAU paper defines:
 *
 *   E1..E13 — for the Moon, Mars satellites
 *   J1..J8  — for Jupiter satellites
 *   Sa..Sh  — for Saturn satellites (the lower-case Greek S series)
 *   U1..U16 — for Uranus satellites
 *   N1..N7  — for Neptune satellites
 *
 * In the table below, each correction is encoded as (target, trig,
 * amplitude, arg_a, arg_b) where the arg variable's identity is
 * implicit in the comment. */

const K26AstroIAURotation K26ASTRO_IAU_ROTATION_TABLE[] = {
    /* ====================================================== */
    /* Sun                                                    */
    /* ====================================================== */
    {
        .name    = "iau2018:sun",
        .naif_id = 10,
        .alpha_0 = 286.13,
        .delta_0 =  63.87,
        .W_0     =  84.176,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 14.18440,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Mercury (with M1..M5 corrections)                      */
    /* M_n = a + b T_c, used in the W-target sin terms.       */
    /* ====================================================== */
    {
        .name    = "iau2018:mercury",
        .naif_id = 199,
        .alpha_0 = 281.0103,
        .delta_0 =  61.4155,
        .W_0     = 329.5988,
        .alpha_rate_per_century = -0.0328,
        .delta_rate_per_century = -0.0049,
        .W_rate_per_day         =   6.1385108,
        .n_corrections          = 5,
        .corrections = {
            /* M1 = 174.7910857 + 4.092335 T_c */
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.00993822,
              .arg_a = 174.7910857, .arg_b =  4.092335 },
            /* M2 = 349.5821714 + 8.184670 T_c */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.00104581,
              .arg_a = 349.5821714, .arg_b =  8.184670 },
            /* M3 = 164.3732571 + 12.277005 T_c */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.00010280,
              .arg_a = 164.3732571, .arg_b = 12.277005 },
            /* M4 = 339.1643429 + 16.369340 T_c */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.00002364,
              .arg_a = 339.1643429, .arg_b = 16.369340 },
            /* M5 = 153.9554286 + 20.461675 T_c */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.00000532,
              .arg_a = 153.9554286, .arg_b = 20.461675 }
        }
    },

    /* ====================================================== */
    /* Venus (retrograde rotation)                            */
    /* ====================================================== */
    {
        .name    = "iau2018:venus",
        .naif_id = 299,
        .alpha_0 = 272.76,
        .delta_0 =  67.16,
        .W_0     = 160.20,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -1.4813688,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Earth (IAU simple model; for high-precision GCRS↔ECEF  */
    /* use the full IAU 2006/2000A pipeline in                 */
    /* earth_orientation.c)                                   */
    /* ====================================================== */
    {
        .name    = "iau2018:earth",
        .naif_id = 399,
        .alpha_0 = 0.00,
        .delta_0 = 90.00,
        .W_0     = 190.147,
        .alpha_rate_per_century = -0.641,
        .delta_rate_per_century = -0.557,
        .W_rate_per_day         = 360.9856235,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Mars (with M1..M14 corrections — Archinal 2018         */
    /* §11, Table 5)                                          */
    /* ====================================================== */
    {
        .name    = "iau2018:mars",
        .naif_id = 499,
        .alpha_0 = 317.269202,
        .delta_0 =  54.432516,
        .W_0     = 176.049863,
        .alpha_rate_per_century = -0.10927547,
        .delta_rate_per_century = -0.05827105,
        .W_rate_per_day         = 350.891982443297,
        .n_corrections          = 17,
        .corrections = {
            /* α target. */
            /* M1 = 198.991226 + 19139.4819985 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.000068,
              .arg_a = 198.991226, .arg_b = 19139.4819985 },
            /* M2 = 226.292679 + 38280.8511281 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.000238,
              .arg_a = 226.292679, .arg_b = 38280.8511281 },
            /* M3 = 249.663391 + 57420.7251593 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.000052,
              .arg_a = 249.663391, .arg_b = 57420.7251593 },
            /* M4 = 266.183510 + 76560.6367950 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.000009,
              .arg_a = 266.183510, .arg_b = 76560.6367950 },
            /* M5 = 79.398797 + 0.5042615 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.419057,
              .arg_a =  79.398797, .arg_b =      0.5042615 },

            /* δ target — uses cos M_n. */
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000051,
              .arg_a = 122.433576, .arg_b = 19139.9407476 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000141,
              .arg_a =  43.058401, .arg_b = 38280.8753272 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000031,
              .arg_a = 57.663379, .arg_b = 57420.7517205 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000005,
              .arg_a = 79.476401, .arg_b = 76560.6495004 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  1.591274,
              .arg_a = 166.325722, .arg_b =     0.5042615 },

            /* W target. */
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.000145,
              .arg_a = 129.071773, .arg_b = 19140.0328244 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.000157,
              .arg_a =  36.352167, .arg_b = 38281.0473591 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.000040,
              .arg_a = 56.668646, .arg_b = 57420.9295360 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.000001,
              .arg_a = 67.364003, .arg_b = 76560.2552215 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.000001,
              .arg_a = 104.792680, .arg_b = 95700.4387578 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.584542,
              .arg_a =  95.391654, .arg_b =     0.5042615 },
            /* Second order W correction. */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.000035,
              .arg_a = 190.500000, .arg_b =     1.0085230 }
        }
    },

    /* ====================================================== */
    /* Jupiter (Ja..Je corrections)                           */
    /* ====================================================== */
    {
        .name    = "iau2018:jupiter",
        .naif_id = 599,
        .alpha_0 = 268.056595,
        .delta_0 =  64.495303,
        .W_0     = 284.95,
        .alpha_rate_per_century = -0.006499,
        .delta_rate_per_century =  0.002413,
        .W_rate_per_day         = 870.5360000,
        .n_corrections          = 10,
        .corrections = {
            /* α: -0.000913 sin Ja - 0.000325 sin Jb -... */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.000913,
              .arg_a =  99.360714, .arg_b = 4850.4046 },     /* Ja */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.000325,
              .arg_a = 175.895369, .arg_b = 1191.9605 },     /* Jb */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.000284,
              .arg_a = 300.323162, .arg_b =  262.5475 },     /* Jc */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.000099,
              .arg_a = 114.012305, .arg_b = 6070.2476 },     /* Jd */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.000266,
              .arg_a =  49.511251, .arg_b =   64.3000 },     /* Je */
            /* δ: cos of same arguments. */
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000389,
              .arg_a =  99.360714, .arg_b = 4850.4046 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000139,
              .arg_a = 175.895369, .arg_b = 1191.9605 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000121,
              .arg_a = 300.323162, .arg_b =  262.5475 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.000042,
              .arg_a = 114.012305, .arg_b = 6070.2476 },
            { .target_idx = 1, .trig_kind = 1, .amplitude = -0.000114,
              .arg_a =  49.511251, .arg_b =   64.3000 }
        }
    },

    /* ====================================================== */
    /* Saturn                                                 */
    /* ====================================================== */
    {
        .name    = "iau2018:saturn",
        .naif_id = 699,
        .alpha_0 =  40.589,
        .delta_0 =  83.537,
        .W_0     =  38.90,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 810.7939024,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Uranus (retrograde rotation)                           */
    /* ====================================================== */
    {
        .name    = "iau2018:uranus",
        .naif_id = 799,
        .alpha_0 = 257.311,
        .delta_0 = -15.175,
        .W_0     = 203.81,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -501.1600928,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Neptune                                                */
    /* ====================================================== */
    {
        .name    = "iau2018:neptune",
        .naif_id = 899,
        .alpha_0 = 299.36,
        .delta_0 =  43.46,
        .W_0     = 249.978,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 541.1397757,
        .n_corrections          = 3,
        .corrections = {
            /* Neptune's α/δ have a single sin N correction:
             *   N = 357.85 + 52.316 T_c */
            { .target_idx = 0, .trig_kind = 0, .amplitude = 0.70,
              .arg_a = 357.85, .arg_b = 52.316 },
            { .target_idx = 1, .trig_kind = 1, .amplitude = -0.51,
              .arg_a = 357.85, .arg_b = 52.316 },
            /* W gets the same N argument. */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.48,
              .arg_a = 357.85, .arg_b = 52.316 }
        }
    },

    /* ====================================================== */
    /* Pluto                                                  */
    /* ====================================================== */
    {
        .name    = "iau2018:pluto",
        .naif_id = 999,
        .alpha_0 = 132.993,
        .delta_0 =  -6.163,
        .W_0     = 302.695,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 56.3625225,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Moon (E1..E13 corrections — Archinal 2018 §11.2)       */
    /* E_n = a + b T_c, used across α, δ, W targets.          */
    /* ====================================================== */
    {
        .name    = "iau2018:moon",
        .naif_id = 301,
        .alpha_0 = 269.9949,
        .delta_0 =  66.5392,
        .W_0     =  38.3213,
        .alpha_rate_per_century =  0.0031,
        .delta_rate_per_century =  0.0130,
        .W_rate_per_day         = 13.17635815,
        .n_corrections          = 24,
        .corrections = {
            /* E1 = 125.045 - 0.0529921 T_c — actually -1934.1361849 deg/century */
            /* IAU table publishes E_n with arg_b in deg/century. */
            /* α: -3.8787 sin E1, -0.1204 sin E2, +0.0700 sin E3, -0.0172 sin E4,
             *    +0.0072 sin E6, -0.0052 sin E10, +0.0043 sin E13 */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -3.8787,
              .arg_a = 125.045, .arg_b = -1934.136185 },     /* E1 */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.1204,
              .arg_a = 250.089, .arg_b = -3868.272370 },     /* E2 */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.0700,
              .arg_a = 260.008, .arg_b = 475263.3328 },     /* E3 */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.0172,
              .arg_a = 176.625, .arg_b = 487269.6797 },     /* E4 */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.0072,
              .arg_a = 357.529, .arg_b =  35999.0500 },     /* E6 */
            { .target_idx = 0, .trig_kind = 0, .amplitude = -0.0052,
              .arg_a = 311.589, .arg_b = 964468.4990 },     /* E10 */
            { .target_idx = 0, .trig_kind = 0, .amplitude =  0.0043,
              .arg_a = 134.963, .arg_b = 477198.8700 },     /* E13 */
            /* δ: cos same arguments. */
            { .target_idx = 1, .trig_kind = 1, .amplitude =  1.5419,
              .arg_a = 125.045, .arg_b = -1934.136185 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.0239,
              .arg_a = 250.089, .arg_b = -3868.272370 },
            { .target_idx = 1, .trig_kind = 1, .amplitude = -0.0278,
              .arg_a = 260.008, .arg_b = 475263.3328 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.0068,
              .arg_a = 176.625, .arg_b = 487269.6797 },
            { .target_idx = 1, .trig_kind = 1, .amplitude = -0.0029,
              .arg_a = 357.529, .arg_b =  35999.0500 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.0009,
              .arg_a = 311.589, .arg_b = 964468.4990 },
            { .target_idx = 1, .trig_kind = 1, .amplitude =  0.0008,
              .arg_a = 134.963, .arg_b = 477198.8700 },
            /* W: sin same. */
            { .target_idx = 2, .trig_kind = 0, .amplitude =  3.5610,
              .arg_a = 125.045, .arg_b = -1934.136185 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.1208,
              .arg_a = 250.089, .arg_b = -3868.272370 },
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.0642,
              .arg_a = 260.008, .arg_b = 475263.3328 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.0158,
              .arg_a = 176.625, .arg_b = 487269.6797 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.0252,
              .arg_a = 357.529, .arg_b =  35999.0500 },
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.0066,
              .arg_a = 311.589, .arg_b = 964468.4990 },
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.0047,
              .arg_a = 134.963, .arg_b = 477198.8700 },
            /* W gets three extra terms — E11, E12, E5. */
            { .target_idx = 2, .trig_kind = 0, .amplitude = -0.0046,
              .arg_a = 357.529, .arg_b = 1934.1361849 + 35999.0500 },     /* E11 */
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.0028,
              .arg_a = 125.045, .arg_b = -3868.2723700 },
            { .target_idx = 2, .trig_kind = 0, .amplitude =  0.0052,
              .arg_a = 297.850, .arg_b = 445267.1115 }      /* E5 */
        }
    },

    /* ====================================================== */
    /* Mars satellites — Phobos + Deimos                       */
    /* ====================================================== */
    {
        .name    = "iau2018:phobos",
        .naif_id = 401,
        .alpha_0 = 317.67071657,
        .delta_0 =  52.88627266,
        .W_0     =  35.18774440,    /* 2019 erratum value */
        .alpha_rate_per_century = -0.10844326,
        .delta_rate_per_century = -0.06134706,
        .W_rate_per_day         = 1128.84475928,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:deimos",
        .naif_id = 402,
        .alpha_0 = 316.65705808,
        .delta_0 =  53.50992033,
        .W_0     =  79.39932954,
        .alpha_rate_per_century = -0.10518014,
        .delta_rate_per_century = -0.05979094,
        .W_rate_per_day         = 285.16188899,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Galilean moons                                          */
    /* ====================================================== */
    {
        .name    = "iau2018:io",
        .naif_id = 501,
        .alpha_0 = 268.05,
        .delta_0 =  64.50,
        .W_0     = 200.39,
        .alpha_rate_per_century = -0.009,
        .delta_rate_per_century =  0.003,
        .W_rate_per_day         = 203.4889538,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:europa",
        .naif_id = 502,
        .alpha_0 = 268.08,
        .delta_0 =  64.51,
        .W_0     =  36.022,
        .alpha_rate_per_century = -0.009,
        .delta_rate_per_century =  0.003,
        .W_rate_per_day         = 101.3747235,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:ganymede",
        .naif_id = 503,
        .alpha_0 = 268.20,
        .delta_0 =  64.57,
        .W_0     =  44.064,
        .alpha_rate_per_century = -0.009,
        .delta_rate_per_century =  0.003,
        .W_rate_per_day         = 50.3176081,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:callisto",
        .naif_id = 504,
        .alpha_0 = 268.72,
        .delta_0 =  64.83,
        .W_0     = 259.51,
        .alpha_rate_per_century = -0.009,
        .delta_rate_per_century =  0.003,
        .W_rate_per_day         = 21.5710715,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Major Saturnian moons                                   */
    /* ====================================================== */
    {
        .name    = "iau2018:mimas",
        .naif_id = 601,
        .alpha_0 =  40.66,
        .delta_0 =  83.52,
        .W_0     = 333.46,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 381.9945550,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:enceladus",
        .naif_id = 602,
        .alpha_0 =  40.66,
        .delta_0 =  83.52,
        .W_0     =   6.32,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 262.7318996,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:tethys",
        .naif_id = 603,
        .alpha_0 =  40.66,
        .delta_0 =  83.52,
        .W_0     =   8.95,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 190.6979085,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:dione",
        .naif_id = 604,
        .alpha_0 =  40.66,
        .delta_0 =  83.52,
        .W_0     = 357.6,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 131.5349316,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:rhea",
        .naif_id = 605,
        .alpha_0 =  40.38,
        .delta_0 =  83.55,
        .W_0     = 235.16,
        .alpha_rate_per_century = -0.036,
        .delta_rate_per_century = -0.004,
        .W_rate_per_day         = 79.6900478,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:titan",
        .naif_id = 606,
        .alpha_0 =  39.4827,
        .delta_0 =  83.4279,
        .W_0     = 186.5855,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 22.5769768,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:iapetus",
        .naif_id = 608,
        .alpha_0 = 318.16,
        .delta_0 =  75.03,
        .W_0     = 350.20,
        .alpha_rate_per_century = -3.949,
        .delta_rate_per_century = -1.143,
        .W_rate_per_day         = 4.5379572,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Uranian moons                                           */
    /* ====================================================== */
    {
        .name    = "iau2018:miranda",
        .naif_id = 705,
        .alpha_0 = 257.43,
        .delta_0 = -15.08,
        .W_0     =  30.70,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -254.6906892,    /* retrograde */
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:ariel",
        .naif_id = 701,
        .alpha_0 = 257.43,
        .delta_0 = -15.10,
        .W_0     = 156.22,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -142.8356681,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:umbriel",
        .naif_id = 702,
        .alpha_0 = 257.43,
        .delta_0 = -15.10,
        .W_0     = 108.05,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -86.8688923,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:titania",
        .naif_id = 703,
        .alpha_0 = 257.43,
        .delta_0 = -15.10,
        .W_0     =  77.74,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -41.3514316,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },
    {
        .name    = "iau2018:oberon",
        .naif_id = 704,
        .alpha_0 = 257.43,
        .delta_0 = -15.10,
        .W_0     =   6.77,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -26.7394932,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Neptune's Triton (retrograde)                           */
    /* ====================================================== */
    {
        .name    = "iau2018:triton",
        .naif_id = 801,
        .alpha_0 = 299.36,
        .delta_0 =  41.17,
        .W_0     = 296.53,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = -61.2572637,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* Pluto's Charon                                          */
    /* ====================================================== */
    {
        .name    = "iau2018:charon",
        .naif_id = 901,
        .alpha_0 = 132.993,
        .delta_0 =  -6.163,
        .W_0     = 122.695,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 56.3625225,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    },

    /* ====================================================== */
    /* NULL terminator                                         */
    /* ====================================================== */
    {
        .name    = NULL,
        .naif_id = 0,
        .alpha_0 = 0.0, .delta_0 = 0.0, .W_0 = 0.0,
        .alpha_rate_per_century = 0.0,
        .delta_rate_per_century = 0.0,
        .W_rate_per_day         = 0.0,
        .n_corrections          = 0,
        .corrections            = {{ 0 }}
    }
};
