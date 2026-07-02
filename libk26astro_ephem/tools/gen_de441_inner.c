/* gen_de441_inner.c — verify-and-copy tool for the DE441 inner-SS
 * kernel. Host-only.
 *
 * Usage:
 *   gen_de441_inner <input.bsp> <output.spk>
 *
 * Reads `input.bsp` (a JPL DE441 SPK file from
 *   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
 *   de441.bsp
 * or a user-extracted inner-planet subset), verifies that each inner
 * planet (Mercury / Venus / Earth / Mars / Sun) can be queried at
 * J2000 TDB, and copies the file byte-for-byte to `output.spk`.
 *
 * v0.1 deliberately doesn't extract a subset — the upstream
 * de441_inner files commonly shipped are already inner-planet-only.
 * Extraction logic for "trim a full DE441 to inner planets" would
 * need to rewrite SPK records; tracked in BACKLOG.md.
 *
 * Output is suitable for installing at
 *   /usr/share/k26astro/ephem/de441_inner.spk
 * (the path libk26astro_ephem's `_load_default` reads). */
#include "k26astro_ephem/ephem.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* NAIF SPICE body ids per Acton 1996. The JPL DE441 SPK kernel
 * contains entries for ids {10, 1..9, 301, 399} but NOT 199/299/499
 * — for moonless planets (Mercury, Venus) the planet centre is
 * identical to the barycentre, and for Mars the offset from the
 * barycentre is ~1.6e-8 of the mass-equivalent (sub-metre at 1 AU),
 * far below our gate tolerance. So we query the barycentre ids:
 *   Sun     = 10
 *   Mercury = 1   (barycentre = planet centre, no moons)
 *   Venus   = 2   (barycentre = planet centre, no moons)
 *   Earth   = 399 (planet centre; NOT 3 = Earth-Moon barycentre)
 *   Mars    = 4   (barycentre ≈ planet centre, error << gate tol) */
#define K26_NAIF_SUN      10
#define K26_NAIF_MERCURY  1
#define K26_NAIF_VENUS    2
#define K26_NAIF_EARTH    399
#define K26_NAIF_MARS     4

static int copy_file_(const char *in_path, const char *out_path)
{
    FILE *fi = fopen(in_path, "rb");
    if (!fi) { fprintf(stderr, "open %s: %s\n", in_path, strerror(errno));
        return -1; }
    FILE *fo = fopen(out_path, "wb");
    if (!fo) { fprintf(stderr, "open %s: %s\n", out_path, strerror(errno));
        fclose(fi); return -1; }
    char buf[64 * 1024];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, fi)) > 0) {
        if (fwrite(buf, 1, n, fo) != n) {
            fprintf(stderr, "write %s: %s\n", out_path, strerror(errno));
            fclose(fi); fclose(fo); return -1;
        }
    }
    fclose(fi);
    fclose(fo);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr,
            "usage: %s <input.bsp> <output.spk>\n"
            "verifies the inner-planet states at J2000 TDB, then copies\n"
            "the input file to the output path. Output is suitable for\n"
            "installation at /usr/share/k26astro/ephem/de441_inner.spk.\n",
            argv[0]);
        return 1;
    }
    const char *in_path  = argv[1];
    const char *out_path = argv[2];

    /* Verify the input. */
    K26AstroEphem *e = k26astro_ephem_load(in_path);
    if (!e) {
        fprintf(stderr, "gen_de441_inner: cannot open %s as an SPK file\n",
                in_path);
        return 2;
    }
    K26AstroEpoch t = k26astro_epoch_j2000_tt();
    const int targets[5] = {
        K26_NAIF_SUN, K26_NAIF_MERCURY, K26_NAIF_VENUS,
        K26_NAIF_EARTH, K26_NAIF_MARS
    };
    const char *names[5] = { "sun", "mercury", "venus", "earth", "mars" };
    int fails = 0;
    for (int i = 0; i < 5; i++) {
        K26AstroStateXV s = k26astro_ephem_body_state(e, targets[i], &t);
        K26V3 r = k26astro_pos_to_m_approx(&s.pos);
        if (r.x == 0.0 && r.y == 0.0 && r.z == 0.0
            && s.vel.x == 0.0) {
            fprintf(stderr,
                "gen_de441_inner: input doesn't cover %s (naif=%d) at J2000\n",
                names[i], targets[i]);
            fails++;
        } else {
            double r_au = sqrt(r.x*r.x + r.y*r.y + r.z*r.z) / 1.495978707e11;
            fprintf(stderr,
                "  %-10s r=%.4f AU  (verified at J2000)\n",
                names[i], r_au);
        }
    }
    k26astro_ephem_close(e);
    if (fails > 0) {
        fprintf(stderr, "gen_de441_inner: %d/%d body queries failed\n",
                fails, 5);
        return 3;
    }

    /* Copy. */
    if (copy_file_(in_path, out_path) != 0) {
        return 4;
    }
    struct stat st;
    if (stat(out_path, &st) == 0) {
        fprintf(stderr,
            "gen_de441_inner: wrote %s (%lld bytes)\n",
            out_path, (long long)st.st_size);
    } else {
        fprintf(stderr, "gen_de441_inner: wrote %s\n", out_path);
    }
    return 0;
}
