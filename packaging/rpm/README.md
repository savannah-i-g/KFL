# RPM packaging (scaffold)

This directory is reserved for `.spec` files and RPM build tooling.
Implementation is deferred — Ubuntu `.deb` is the only target with a
confirmed consumer in v0.1.

When RPM packaging lands, it should follow the same per-library
decomposition as the Debian source (`../debian/control`): one source
RPM (`kfl-stack`) producing parallel `libk26<name>0` / `libk26<name>-devel`
binary RPMs.

The top-level `make rpm` target prints this README's pointer and exits
without action.
