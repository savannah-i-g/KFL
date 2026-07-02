# Getting started

This guide walks through the three kinds of KFL program — compute, plot,
and simulation — building and running each one. It assumes `kflc` is on
your `PATH` (see the [README](../README.md) for installation) and that a
C++ compiler is available, since `kflc` invokes it to build the program it
emits.

## Compiling a program

`kflc program.kfl -o program` produces a native executable. A program that
uses only the language core links against the C math library alone:

```
kflc program.kfl -o program            # KFLC_LDLIBS defaults to -lm
./program
```

A program that uses a KFL library names the libraries to link through the
`KFLC_LDLIBS` environment variable, and — if the headers and archives are
not on the default search path — adds `-I`/`-L` paths through
`KFLC_CFLAGS`:

```
KFLC_LDLIBS="-lk26plot -lk26compute -lk26m3d -lcairo -lm" \
    kflc program.kfl -o program
```

Every example under `kflc/examples/` states its exact build line in its
header comment.

## 1. A compute program

Compute programs are pure functions plus a `run` entry point that prints
results. Save this as `kepler.kfl`:

```
form KEPLER
    # Kepler's third law in solar units (mu = 4 pi^2): T = sqrt(a^3).
    fn double period_years(double a_au)
        return sqrt(a_au * a_au * a_au)
    end

    fn void run()
        print "planet   a[AU]    T[yr]"
        print "Earth    1.0000   ", period_years(1.0000)
        print "Mars     1.5237   ", period_years(1.5237)
        print "Jupiter  5.2029   ", period_years(5.2029)
    end
end
```

```
kflc kepler.kfl -o kepler && ./kepler
planet   a[AU]    T[yr]
Earth    1.0000   1
Mars     1.5237   1.88083
Jupiter  5.2029   11.8677
```

`print` takes comma-separated arguments; string literals are emitted
verbatim and numeric expressions are evaluated and formatted. The language
provides scalar math (`sqrt`, `sin`, `exp`, `pow`, ...), `vector` and
`matrix` types with element indexing, and `while`/`if` control flow. See
[`kflc/GRAMMAR.md`](../kflc/GRAMMAR.md) for the complete surface.

A fuller compute example is `kflc/examples/compute_kepler.kfl`.

## 2. A plotting program

Plotting programs pair `fn data` producers with top-level `plot`
declarations. A producer fills `vector`s and hands each to a `series_*`
statement; the `plot` renders one PNG and one SVG:

```
form CURVE
    fn data parabola
        let xs: vector = linspace(-3, 3, 200)
        let ys: vector = zeros(200)
        let i: int = 0
        while i < 200
            ys[i] = xs[i] * xs[i]
            i = i + 1
        end
        series_line "y = x^2" xs ys
    end

    plot p_parabola data parabola
        title    "A parabola"
        x_label  "x"
        y_label  "y"
end
```

```
KFLC_LDLIBS="-lk26plot -lk26compute -lk26m3d -lcairo -lm" \
    kflc curve.kfl -o curve
./curve                                # writes p_parabola.png, p_parabola.svg
```

Series kinds are `series_line`, `series_scatter`, `series_errorbar`,
`series_histogram`, `series_box`, and `series_heatmap` (which takes a
single `matrix` and renders a 2-D color field). The worked examples
`kflc/examples/orbital_mechanics.kfl` (heliocentric orbits + a Hohmann
transfer) and `kflc/examples/porkchop.kfl` (an Earth-to-Mars launch-window
heatmap, solving Lambert's problem across a departure/time-of-flight grid)
show the plotting surface against real physics.

## 3. A simulation program

Simulation programs build a physical world in one or more `fn world`
blocks, advance it, and print observations. The runtime integrates an
N-body gravitational world and can report an observer's line-of-sight
range and direction to a body. This shape draws on the astrodynamics
libraries, so it links the astro stack:

```
KFLC_LDLIBS="-lk26astro_rt -lk26astro_vehicle -lk26astro_grav \
    -lk26astro_conics -lk26astro_atmos -lk26astro_body -lk26astro_ephem \
    -lk26astro_core -lk26compute -lk26tick -lk26m3d -lgfortran -lm" \
    kflc kflc/examples/astro_apophis_2029.kfl -o apophis
./apophis
```

The shipped simulations —
`kflc/examples/astro_apophis_2029.kfl` (the 2029 Apophis Earth flyby) and
`kflc/examples/astro_juno_flyby.kfl` — are the reference programs for the
`fn world` surface: they declare bodies with `astro_body`, advance the
world, and print `observe` results. Programs that call the astrodynamics
libraries also need each library's builtin manifest (`.kflbi`) available,
which the corresponding `-dev` package installs to
`/usr/share/kflc/builtins/`; override the search location with the
`K26_KFL_BUILTIN_PATH` environment variable for an uninstalled tree.

## Where to go next

- [`kflc/GRAMMAR.md`](../kflc/GRAMMAR.md) — the complete language reference.
- [`docs/libraries.md`](libraries.md) — the scientific library suite and
  its KFL builtin surface.
- [`kflc/examples/`](../kflc/examples/) — runnable programs for every
  language feature.
