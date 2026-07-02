# Homebrew packaging (scaffold)

This directory is reserved for a Homebrew Formula and tap configuration.
Implementation is deferred — Ubuntu `.deb` is the only target with a
confirmed consumer in v0.1.

When Homebrew support lands, a single `kfl-stack.rb` formula that
builds the whole tree via `make all && make install PREFIX=#{prefix}`
is the expected shape. Per-library decomposition is not how Homebrew
formulae are typically structured.

Fortran-on-macOS notes: GNU gfortran availability differs between
Intel and Apple Silicon hosts; the formula will need to depend on
`gcc` (which provides gfortran on macOS) and may need conditional
build flags for the four Fortran-backed libraries
(libk26astro_quad / ode / geomag / atmos).

The top-level `make brew` target prints this README's pointer and exits
without action.
