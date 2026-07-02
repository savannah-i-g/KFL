# KFL — the KFL programming language

KFL is a small, block-structured language for writing **artifact-producing
batch programs**: numerical computations, plots, and physical simulations
that run to completion, print their results, and write image files. It is
not a UI language — a KFL program has no window, no widgets, and no event
loop.

A source file (extension `.kfl`) is compiled by `kflc`:

```
kflc program.kfl -o program        # emit + compile to a native binary
./program                          # run it
```

`kflc` translates the `.kfl` into standalone C++, then invokes the system
C++ toolchain to produce a native executable. The generated program runs
any simulations, renders any plots, calls the program entry function, and
exits. Link libraries are supplied through the `KFLC_LDLIBS` environment
variable (see *Compilation and the generated program*).

KFL is hand-authored, block-structured, and indentation-tolerant. It is
not JSON, not YAML, not INI.

---

## Lexical structure

- **Comments** — `#` to end of line. No block-comment form.
- **Strings** — double-quoted, with C-style escapes (`\n`, `\t`, `\"`,
  `\\`). An embedded NUL is forbidden. UTF-8 inside a string is preserved
  verbatim.
- **Identifiers** — `[a-zA-Z_][a-zA-Z0-9_]*`, case-sensitive.
- **Integers** — optional leading `-`, then decimal digits. No hex,
  octal, or binary.
- **Floats** — integer `.` integer, e.g. `0.5`, `3.141592653589793`.
  Scientific notation (`1.32712440018e20`) is accepted for the
  named-argument values of simulation statements.

Whitespace separates tokens; newlines terminate statements. Indentation
is cosmetic and ignored by the parser. The lexer also recognises a few
legacy token forms (`WxH` sizes, `#rrggbb` colours, `Ctrl+Q`-style
shortcuts) retained for backward compatibility; the current language has
no construct that consumes them.

---

## Program structure

Every file is exactly one `form` block:

```
form <NAME>
    <form-level construct>...
end
```

`<NAME>` is an identifier that names the program. There is **no window
header** — the language carries no title, size, or configuration and
produces no graphical surface. (`title`, `size`, and `cfg` lines are still
accepted so that older files parse, but they are ignored.)

The constructs allowed directly inside a `form` are:

| Construct                    | Purpose                                                            |
|------------------------------|--------------------------------------------------------------------|
| `fn <type> <name>(<args>)`   | A typed function. `fn void run()` is the program entry point.       |
| `fn data <name>`             | A plot-data producer (builds series for a `plot`).                 |
| `fn world <name>`            | A simulation callback run against a fresh world.                  |
| `plot <name> data <fn>`      | A figure declaration; renders to `<name>.png` and `<name>.svg`.    |
| `arena <name> capacity ...`  | A form-lifetime bump allocator (see *Memory model*).              |
| `arg <name> [default <v>]`   | A command-line-settable global.                                   |
| `frame <name> ...`           | A named reference frame for the simulation surface.               |
| `epoch <name> "<ISO>" <sc>`  | A compile-time epoch constant.                                    |
| `tick <handler> interval_ms N` | A periodic simulation callback.                                 |

Order is free. A program need only contain what it uses: a pure-compute
program is one or more `fn`s and a `fn void run()`; a plotting program adds
`fn data` producers and `plot` declarations; a simulation adds a
`fn world`.

---

## Functions

```
fn <return-type> <name>(<arg>, <arg>, ...)
    <statement>...
end
```

Return and argument types are `void`, `double`, `int`, `bool`, `string`,
`vector`, `matrix`, or a registered **opaque** type (see *Simulation
surface*). Each argument is `<type> <name>`, optionally prefixed with an
ownership qualifier (`own` / `borrow` / `ptr`; see *Memory model*).

```
fn double period_years(double a_au)
    return sqrt(a_au * a_au * a_au)
end
```

### The `run` entry point

A parameterless `fn void run()` is the conventional entry point. After the
program has executed any `fn world` simulations and rendered any `plot`
figures, the generated `main()` calls `run()` — it is where a compute
program does its work and prints results.

```
fn void run()
    print "Earth period T = ", period_years(1.0), " yr"
end
```

### Statements

| Statement                         | Form                                                     |
|-----------------------------------|----------------------------------------------------------|
| Binding                           | `let <name>: <type> = <expr>`                             |
| Constant binding                  | `const <name>: <type> = <expr>`                          |
| Assignment                        | `<name> = <expr>`                                         |
| Indexed assignment                | `<name>[<expr>] = <expr>`                                 |
| Conditional                       | `if <expr>` … [`else`] … `end`                            |
| Loop                              | `while <expr>` … `end`                                    |
| Return                            | `return [<expr>]`                                         |
| Print                             | `print <arg>, <arg>, ...`                                 |
| Expression statement              | `<expr>`                                                  |

The type annotation on `let` / `const` is written after a colon; when
omitted the binding defaults to `double`. Blocks (`if`, `else`, `while`)
are closed with `end`.

```
let n: int = 240
let mu: double = 4 * 3.141592653589793 * 3.141592653589793
let i: int = 0
while i < n
    x[i] = cos(i)
    i = i + 1
end
```

### Expressions

The expression sub-language covers the usual arithmetic and comparison
operators over `int` and `double`, calls to user `fn`s, and a whitelist of
built-ins:

- **Scalar maths** — `sqrt`, `sin`, `cos`, `tan`, `atan2`, `exp`, `log`,
  `log10`, `pow`, `abs`, `fabs`, `floor`, `ceil`, `min`, `max`, `fmod`.
- **String helpers** — `strlen`, `streq`, `starts_with`, `ends_with`,
  `concat`.
- **Vector / matrix helpers** — statistics and linear algebra such as
  `mean`, `std`, `dot`, `norm`, `mat_mul`, `mat_vec`, `transpose`,
  `solve`. Vector/matrix arguments must be in-scope bindings referenced by
  name, not anonymous temporaries.

### Vectors

A `vector` binding is initialised with a constructor — `zeros(N)`,
`ones(N)`, `linspace(a, b, N)`, `arange(a, b, step)`, or the literal form
`[ ... ]` — then read and written by index (indexing is bounds-checked at
run time):

```
let xs: vector = linspace(0.3, 35, 200)
let ys: vector = zeros(200)
let i: int = 0
while i < 200
    ys[i] = sqrt(1.0 / xs[i])
    i = i + 1
end
```

### Matrices

A `matrix` binding is a row-major 2-D array of doubles. It is initialised
either with the runtime constructor `zeros(rows, cols)` or with a nested
literal `[[ ... ], [ ... ]]`, and its cells are read and written with two
indices (both bounds-checked at run time):

```
let grid: matrix = zeros(300, 240)
grid[i][j] = value
```

Matrices back the `series_heatmap` plotting statement (see *Plotting*).

---

## The `print` statement

`print` writes one line to standard output. Its arguments are
comma-separated; each argument is either a **string literal** (printed
verbatim) or a **numeric expression** (evaluated and printed as a number):

```
print "planet      a[AU]   T[yr]"
print "Mercury  ", 0.3871, "   ", period_years(0.3871)
```

Commas that occur inside a string literal, or inside parentheses or
brackets of a call such as `f(g(a, b))`, do **not** split arguments — only
top-level commas separate them. Each `print` statement emits exactly one
line, terminated by a newline.

---

## Plotting

A plotting program has two parts: `fn data` producers that fill series,
and top-level `plot` declarations that render them.

A `fn data <name>` producer builds series and hands each to a `series_*`
statement (a full producer appears in the plotting example below). Most
kinds take a pair of x/y vectors — `series_<kind> "<label>" <xs> <ys>`,
where `<xs>` and `<ys>` name in-scope `vector` bindings — and cover
`series_line`, `series_scatter`, `series_errorbar`, `series_histogram`,
and `series_box`. The `series_heatmap "<label>" <matrix>` statement
instead takes a single `matrix` binding and renders it as a 2-D color
field over the index-space extent `[0, cols] x [0, rows]`, with the color
scale auto-fitted to the data range.

A top-level `plot` binds a producer to a figure and its labels:

```
plot p_speed data speed_curve
    title    "Circular orbital speed vs distance"
    x_label  "r (AU)"
    y_label  "v (km/s)"
```

Each `plot <name>` renders to `<name>.png` and `<name>.svg` in the working
directory via `libk26plot`. A plotting program must therefore link the
plotting libraries, for example:

```
KFLC_LDLIBS="-lk26plot -lk26compute -lk26m3d -lcairo -lm" \
    kflc orbital_mechanics.kfl -o orbital
```

---

## Memory model

KFL bindings and function arguments carry optional ownership qualifiers,
enforced by the compiler:

| Qualifier | Meaning                                                                    |
|-----------|----------------------------------------------------------------------------|
| `own`     | Single ownership. Assigning one `own` binding to another requires `move()`. |
| `borrow`  | Read-only reference; emitted as `const`. It may not outlive its source.     |
| `ptr`     | Raw passthrough; user-managed lifetime, no compile-time lifetime check.     |

`move(x)` makes an ownership transfer explicit and invalidates the source;
reading a moved-from binding within the same block is a compile error, and
returning a `borrow` of a function-local binding is a compile error.

```
fn double demo_move(own body src)
    let dst: own body = move(src)   # ownership transfers; src invalidated
    return 1.0
end

fn double demo_borrow(borrow body src)
    let view: borrow body = src     # read-only view for the fn lifetime
    return 2.0
end
```

### Arenas

A form-level arena is a bump allocator that lives for the program's
lifetime:

```
arena world_arena  capacity 64 MB
arena scratch      capacity 16 KB
arena small_buffer capacity 4096 reset_mode manual
```

`capacity` is a byte count with an optional `KB` / `MB` / `GB` suffix. A
function opts into an arena with an `allocator` binding in its prologue;
allocations inside the function then draw from that arena:

```
fn double demo()
    allocator = scratch
    return 1.0
end
```

The default reset mode, `fn_exit`, resets the arena automatically on every
return path from a function that binds it. `reset_mode manual` leaves the
caller responsible for resetting the arena when its data becomes stale.

---

## Simulation surface

The compiler ships no physics of its own. Numerical libraries publish
opaque types and function bindings through `.kflbi` manifest files, which
`kflc` loads at start-up from its builtins directory (default
`/usr/share/kflc/builtins/`, overridable with `KFLC_BUILTINS_PATH`). Each
manifest registers:

- **Opaque types** — nominal handle types such as `world`, `body`, and
  `epoch`. The compiler tracks them by name and rejects cross-type
  confusion. Once registered, an opaque name may be used in any type
  position (`let e: body = ...`, `fn world`, function arguments).
- **Builtins** — named functions with a fixed arity, callable from KFL
  expressions and simulation statements.

With the astronomy libraries installed, a program can drive a physical
simulation.

### `fn world` and simulation statements

A `fn world <name>` is a callback executed against a freshly created world.
Its body admits the ordinary statements plus a declarative simulation
surface:

- `astro_body <name> <key>=<value> ...` — declare a body from named
  parameters (`gm`, `mass`, `radius`, `parent`, …).
- `step <dt-expr>` — advance the whole world by `dt` seconds.
- `propagate <body> for <dt-expr>` — advance a single body.
- `for_each <name> in <world>` … `end` — read-only iteration over bodies.
- `observe <target> from <observer> [key=value ...]` — compute an
  observation and print it to standard output.

```
form APOPHIS2029

    frame ICRF inertial
    epoch J2025 "2025-01-01T00:00:00" UTC

    fn world apophis_demo
        astro_body sun   gm=1.32712440018e20 mass=1.989e30 radius=6.957e8
        astro_body earth gm=3.986004418e14   mass=5.972e24 radius=6.371e6 parent=sun
        astro_body apophis mass=2.7e10 gm=1.8e0 radius=185.0 parent=sun

        step 3600.0
        observe apophis from earth mode=apparent
    end

end
```

*(This is the shape of `examples/astro_apophis_2029.kfl`.)*

### Frames, epochs, ticks

- `frame <name> inertial` or `frame <name> body_fixed body <ident>`
  declares a named reference frame used by the runtime's coordinate
  transforms.
- `epoch <name> "<ISO-8601>" <scale>` binds a compile-time epoch constant.
  `<scale>` is one of `TAI`, `UTC`, `UT1`, `TT`, `TDB`. The ISO-8601 string
  is parsed at compile time, never at run time.
- `tick <handler> interval_ms <N>` registers a `fn` as a periodic callback,
  invoked every `N` milliseconds while a simulation runs.

A simulation program links the astronomy libraries via `KFLC_LDLIBS`, for
example:

```
KFLC_LDLIBS="-lk26astro_rt -lk26astro_grav -lk26astro_conics \
    -lk26astro_body -lk26astro_ephem -lk26astro_core \
    -lk26compute -lk26tick -lk26m3d -lgfortran -lm" \
    kflc astro_apophis_2029.kfl -o apophis
```

---

## Command-line arguments

`arg <name> [default <value>] [mutable | readonly]` declares a global that
the program reads by name and the user sets on the command line. The type
is inferred from the default: a boolean default exposes `--<name>` /
`--no-<name>`, while a numeric or string default is set with
`--<name> <value>`. A `readonly` arg keeps its default and exposes no flag.

---

## Compilation and the generated program

`kflc program.kfl -o program` emits a self-contained C++ translation unit
and compiles it to a native binary. Useful modes:

| Mode          | Effect                                              |
|---------------|-----------------------------------------------------|
| `-o <file>`   | Emit, compile, and link a binary at `<file>`.       |
| `-c`          | Also keep the emitted `.kflc.cc` beside the binary. |
| `--emit`      | Write the generated C++ to stdout; do not compile.  |
| `--dump`      | Dump the parsed AST to stdout.                       |
| `--check`     | Parse and emit silently; exit non-zero on error.    |

The emitted program is an ordinary `main()`. It:

1. parses declared `arg`s from `argv`;
2. initialises any form-level arenas;
3. runs each `fn world` against a fresh world (in a portable,
   bit-reproducible mode), whose `observe` statements print to stdout;
4. renders each `plot` to `<name>.png` and `<name>.svg`;
5. calls `run()` if the program defines one;
6. returns 0.

Compilation is controlled by three environment variables: `CXX` (the C++
compiler, default `c++`), `KFLC_CFLAGS` (compile flags), and `KFLC_LDLIBS`
(link libraries, default `-lm`). A program that uses a numerical library
must extend `KFLC_LDLIBS` with the matching `-l` flags.

---

## Worked example — a compute program

A headless program that prints Keplerian quantities. It uses only scalar
maths, so the default `KFLC_LDLIBS=-lm` suffices. This is the shape of
`examples/compute_kepler.kfl`:

```
form COMPUTEKEPLER

    # Orbital period in years for a circular orbit of radius a (AU).
    fn double period_years(double a_au)
        return sqrt(a_au * a_au * a_au)
    end

    # Circular orbital speed (AU/yr) at radius a (AU).
    fn double circular_speed(double a_au)
        let mu: double = 4 * 3.141592653589793 * 3.141592653589793
        return sqrt(mu / a_au)
    end

    fn void run()
        print "planet      a[AU]   T[yr]     v[AU/yr]"
        print "Mercury  ", 0.3871, "   ", period_years(0.3871), "  ", circular_speed(0.3871)
        print "Earth    ", 1.0000, "   ", period_years(1.0000), "  ", circular_speed(1.0000)
        print "Mars     ", 1.5237, "   ", period_years(1.5237), "  ", circular_speed(1.5237)
    end

end
```

```
kflc compute_kepler.kfl -o kepler
./kepler
```

---

## Worked example — a plotting program

A program that samples a curve, marks a few points, and writes a figure.
It follows the shape of `examples/orbital_mechanics.kfl` (which combines
several such producers and plots with a `run()`):

```
form ORBITAL

    fn double circular_speed(double a_au)
        let mu: double = 4 * 3.141592653589793 * 3.141592653589793
        return sqrt(mu / a_au)
    end

    # Producer: a continuous curve plus discrete planet markers.
    fn data speed_curve
        let rs: vector = linspace(0.3, 35, 200)
        let vs: vector = zeros(200)
        let mu: double = 4 * 3.141592653589793 * 3.141592653589793
        let i: int = 0
        while i < 200
            vs[i] = sqrt(mu / rs[i])
            i = i + 1
        end
        series_line "v_circ" rs vs

        let pr: vector = zeros(3)
        let pv: vector = zeros(3)
        pr[0] = 0.3871
        pr[1] = 1.0000
        pr[2] = 1.5237
        let j: int = 0
        while j < 3
            pv[j] = sqrt(mu / pr[j])
            j = j + 1
        end
        series_scatter "planets" pr pv
    end

    plot p_speed data speed_curve
        title    "Circular orbital speed vs heliocentric distance"
        x_label  "r (AU)"
        y_label  "v (AU/yr)"

    fn void run()
        print "Earth circular speed = ", circular_speed(1.0), " AU/yr"
    end

end
```

```
KFLC_LDLIBS="-lk26plot -lk26compute -lk26m3d -lcairo -lm" \
    kflc orbital.kfl -o orbital
./orbital           # writes p_speed.png and p_speed.svg
```
