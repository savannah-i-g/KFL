# libk26util

General-purpose C utilities for the KFL stack: string helpers
exposed to KFL's expression layer, POSIX filesystem primitives, and
atomic file I/O.

## Public surface

### `<k26_str.h>` string helpers

KFL-callable string utilities backed by the kflc compiler's
expression builtins (`strlen`, `streq`, `starts_with`, `ends_with`,
`concat`). Each function tolerates NULL inputs by returning the
empty answer for its return type, so KFL form-arguments may be
unset without crashing.

### `<k26_fsops.h>` filesystem mutation primitives

Pure POSIX filesystem operations: recursive delete, copy, and move;
size and item walk for progress estimation; human-readable byte
formatting; collision-avoiding unique-name generation; trash routing
through `~/.local/share/k26/trash` with sidecar `.trashinfo`
metadata. Recursion uses `opendir`, `readdir`, and `lstat` rather
than `nftw`, and never follows symlinks into a target tree. No FLTK
or X11 dependency.

### `<k26_atomic_io.h>` atomic file I/O

Atomic write-then-rename primitives for crash-safe configuration and
small-file updates.

## Dependencies

C compiler (GCC 9 or later, Clang 12 or later), GNU make, libc. No
external libraries.

## Build and install

Built as part of the KFL_Stack distribution by the top-level
Makefile (`bedrock` phase). Direct builds:

```
make
sudo make install                # default PREFIX=/usr/local
```

Installed footprint: three headers under `$PREFIX/include/`, a
static archive `libk26util.a`, and a shared object with conventional
soname versioning (`libk26util.so.0.1.0` plus `.so.0` and `.so`
symlinks).

## License

MIT. See the top-level `LICENSE` file in the KFL_Stack distribution.
