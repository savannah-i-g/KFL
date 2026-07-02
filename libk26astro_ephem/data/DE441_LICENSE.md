# DE441 Inner-Planet Ephemeris — Source & Attribution

## What this file is

`de441_inner.spk` ships the JPL DE441 planetary ephemeris, restricted
to the inner solar-system bodies (Sun, Mercury, Venus, Earth, Mars).
It is a NAIF/SPICE Type-2 SPK kernel — Chebyshev-polynomial position
records relative to the solar-system barycentre, valid over the DE441
nominal coverage (1550..2650 CE for the inner-planet records).

The K26 build pipeline does **not** rebundle the entire DE441
distribution. Savannah supplies an upstream binary (the canonical
NASA NAIF distribution, or an offline-extracted subset); the
`tools/gen_de441_inner` host tool verifies that the inner-planet
queries succeed at J2000 TDB and copies the file byte-for-byte to
the K26 tree at `libk26astro_ephem/data/de441_inner.spk`. The K26
data apk (`libk26astro_ephem-data`) installs it to
`/usr/share/k26astro/ephem/de441_inner.spk` — the path that
`k26astro_ephem_load_default()` reads.

## Source

The canonical upstream is:

> NAIF / Solar System Dynamics Group, JPL
> https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/

The file commonly referred to as DE441 lives there as `de441.bsp`
(full, ~3 GB) and as smaller per-time-range subsets. Inner-planet
subsets are also distributed for ephemeris-only applications.

## Publication

Folkner W.M., Park R.S. (2014):
"The Planetary and Lunar Ephemerides DE430 and DE431",
IPN Progress Report 42-196.

Park R.S., Folkner W.M., Williams J.G., Boggs D.H. (2021):
"The JPL Planetary and Lunar Ephemerides DE440 and DE441",
The Astronomical Journal 161:105.

## License

The DE441 ephemerides are produced by the U.S. Government (Jet
Propulsion Laboratory, California Institute of Technology, under
contract with NASA) and are in the public domain under U.S. law.
Redistribution is permitted; attribution to NASA / JPL is requested
but not required by law.

K26's copy ships under the same public-domain status; downstream
users are encouraged to keep the Folkner / Park citations in any
publication that depends on these ephemerides.
