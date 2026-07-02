# KFL_Stack — host-native orchestrator.
#
# Builds the KFL compiler and the portable K26 library stack in
# topological order. Host build (glibc + gfortran), no podman sidecar.
# Fortran tier (QUADPACK / ODEPACK / IGRF / NRLMSISE) is conditional:
# skipped cleanly when gfortran is absent.
#
# Per-lib Makefiles are the canonical build units; this file just
# orchestrates them. Override PREFIX / LIBDIR / DESTDIR as needed:
#   make install PREFIX=/usr/local                 # host install
#   make install PREFIX=/usr LIBDIR=/usr/lib/x86_64-linux-gnu DESTDIR=...
# The .deb build (packaging/debian/rules) passes the multi-arch
# LIBDIR; manual builds default to $PREFIX/lib.

PREFIX  ?= /usr/local
DESTDIR ?=
LIBDIR  ?= $(PREFIX)/lib

VERSION         := $(shell cat VERSION 2>/dev/null || echo 0.0.0)
KFL_GRAMMAR_VER := $(shell cat KFL_GRAMMAR_VERSION 2>/dev/null || echo 0.0)
GIT_SHA         := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)

# ----- Build dependency ordering (topological) ---------------------

BEDROCK     = libk26util libk26m3d libk26compute libk26tick
AUX_ALWAYS  = libk26plot libk26geo
AUX_CURL    = libk26http   # needs libcurl4-openssl-dev; skipped if absent
AUX         = $(AUX_ALWAYS) $(AUX_CURL)
ASTRO_CORE  = libk26astro_core

# Topological order — each lib's .so links cross-K26-lib .so deps, so
# producers must build before consumers. Ephem before body; conics +
# grav after body; fit standalone (pure-C, no deps).
ASTRO_TIER1 = libk26astro_ephem libk26astro_body libk26astro_conics \
              libk26astro_grav libk26astro_fit libk26astro_vehicle

# Fortran-backed libs — atmos goes last (links to astro_core).
ASTRO_FORTRAN = libk26astro_quad libk26astro_ode libk26astro_geomag \
                libk26astro_atmos

# Runtime manager — links astro_grav + atmos + conics + body + ephem +
# vehicle + core + compute + tick + m3d + libgfortran. Must build
# last, and requires the Fortran tier (atmos.so) to exist.
ASTRO_RT    = libk26astro_rt

# The mission-design and spacecraft-subsystem tier — propulsion,
# guidance/navigation/control, thermal, sensors, trajectory design,
# mission timelines, electrical power, communications, and
# atmospheric-flight aerodynamics — lives in the sibling KFL_Missions
# repository. Nothing in this stack depends on it.

# Every C library that builds unconditionally (no Fortran tier).
ALL_C_LIBS = $(BEDROCK) $(AUX) $(ASTRO_CORE) $(ASTRO_TIER1)

# Optional out-of-tree library overlay. A sibling project (e.g. the
# defense subsystem stack) may drop a kfl_overlay.mk fragment here to
# append extra library directories onto a tier. The fragment is not
# tracked and is absent from release tarballs; when present it is
# expected to add its dirs to ASTRO_TIER1_EXTRA so they build after
# the core composition tier they depend on.
ASTRO_TIER1_EXTRA =
-include kfl_overlay.mk
ASTRO_TIER1 += $(ASTRO_TIER1_EXTRA)

EXPORTS = PREFIX=$(PREFIX) LIBDIR=$(LIBDIR) DESTDIR=$(DESTDIR)

.PHONY: all bedrock aux astro compiler tools test test-tick test-astro test-compiler \
        install uninstall \
        deb deb-lintian deb-clean rpm brew dist \
        clean clean-host-artifacts distclean \
        kfl-version kfl-grammar-version help \
        check-gfortran

# ----- Top-level ---------------------------------------------------

all: bedrock aux astro compiler
	@echo ""
	@echo "KFL_Stack $(VERSION) (grammar v$(KFL_GRAMMAR_VER), git $(GIT_SHA))"
	@echo "build complete. Run 'make test' to verify."

help:
	@echo "KFL_Stack v$(VERSION) (grammar v$(KFL_GRAMMAR_VER))"
	@echo ""
	@echo "Common targets:"
	@echo "  make                build everything in topological order"
	@echo "  make test           run unit + integration tests"
	@echo "  make install        install to \$$DESTDIR\$$PREFIX (default /usr/local)"
	@echo "  make uninstall      remove installed files (manifest-driven)"
	@echo "  make deb            build .deb packages (see packaging/debian/)"
	@echo "  make rpm            [stub — see packaging/rpm/README.md]"
	@echo "  make brew           [stub — see packaging/brew/README.md]"
	@echo "  make dist           produce KFL_Stack-\$$VERSION.tar.gz"
	@echo "  make clean          remove build outputs per-lib"
	@echo "  make distclean      clean + remove dist/"
	@echo "  make kfl-version    print version"
	@echo ""
	@echo "Component phases: bedrock aux astro compiler"

# ----- Build phases ------------------------------------------------

bedrock:
	@for l in $(BEDROCK); do \
	    echo "==> $$l"; $(MAKE) -C $$l $(EXPORTS) || exit 1; \
	done

aux: bedrock
	@for l in $(AUX_ALWAYS); do \
	    echo "==> $$l"; $(MAKE) -C $$l $(EXPORTS) || exit 1; \
	done
	@if echo '#include <curl/curl.h>' | $(CC) -E -x c - >/dev/null 2>&1; then \
	    for l in $(AUX_CURL); do \
	        echo "==> $$l (libcurl)"; \
	        $(MAKE) -C $$l $(EXPORTS) || exit 1; \
	    done; \
	else \
	    echo ""; \
	    echo "==> SKIP $(AUX_CURL) — libcurl headers not found"; \
	    echo "    install with: sudo apt-get install libcurl4-openssl-dev"; \
	    echo ""; \
	fi

astro: bedrock aux
	@echo "==> $(ASTRO_CORE)"
	@$(MAKE) -C $(ASTRO_CORE) $(EXPORTS)
	@for l in $(ASTRO_TIER1); do \
	    echo "==> $$l"; $(MAKE) -C $$l $(EXPORTS) || exit 1; \
	done
	@if command -v gfortran >/dev/null 2>&1; then \
	    for l in $(ASTRO_FORTRAN); do \
	        echo "==> $$l (Fortran)"; \
	        $(MAKE) -C $$l $(EXPORTS) || exit 1; \
	    done; \
	    echo "==> $(ASTRO_RT)"; \
	    $(MAKE) -C $(ASTRO_RT) $(EXPORTS) || exit 1; \
	else \
	    echo ""; \
	    echo "==> SKIP $(ASTRO_FORTRAN) $(ASTRO_RT) — gfortran not found"; \
	    echo "    rt is skipped because it links atmos which needs Fortran"; \
	    echo ""; \
	fi

compiler:
	@echo "==> kflc"
	$(MAKE) -C kflc

# Host-side tooling. k26astro_replay is the REFERENCED-mode op-log
# inspector consumed by KFL_Verify replay claims and by
# tests/astro/test_referenced_replay.c (which expects the binary at
# tools/k26astro_replay/k26astro_replay relative to the test). Requires
# the full astro stack to be linkable, so build after astro.
tools: astro
	@echo "==> tools/k26astro_replay"
	@$(MAKE) -C tools/k26astro_replay

# ----- Tests -------------------------------------------------------

test: test-tick test-compiler test-astro

# Per-lib unit tests for libs that have a 'test' target.
test-units: bedrock aux $(ASTRO_CORE)
	@for l in libk26astro_core libk26astro_body libk26astro_ephem \
	          libk26astro_conics libk26astro_grav libk26astro_rt \
	          libk26astro_fit libk26astro_vehicle; do \
	    if [ -f $$l/Makefile ] && grep -q '^test:' $$l/Makefile; then \
	        echo "==> $$l test"; \
	        $(MAKE) -C $$l test $(EXPORTS) || exit 1; \
	    fi; \
	done

test-tick: bedrock
	@if [ -d tests ]; then $(MAKE) -C tests tick; fi

test-astro: astro
	@if [ -d tests/astro ]; then $(MAKE) -C tests astro; fi
	@if command -v gfortran >/dev/null 2>&1; then \
	    for l in $(ASTRO_FORTRAN); do \
	        if [ -f $$l/Makefile ] && grep -q '^test:' $$l/Makefile; then \
	            echo "==> $$l test"; \
	            $(MAKE) -C $$l test $(EXPORTS) || exit 1; \
	        fi; \
	    done; \
	fi
	@if [ -d tests/integration ]; then $(MAKE) -C tests integration; fi

test-compiler: compiler
	$(MAKE) -C kflc test

# ----- Install / uninstall -----------------------------------------

install: all
	@for l in $(ALL_C_LIBS); do \
	    if grep -q '^install:' $$l/Makefile 2>/dev/null; then \
	        $(MAKE) -C $$l install $(EXPORTS) || exit 1; \
	    fi; \
	done
	@if command -v gfortran >/dev/null 2>&1; then \
	    for l in $(ASTRO_FORTRAN); do \
	        $(MAKE) -C $$l install $(EXPORTS) || exit 1; \
	    done; \
	    $(MAKE) -C $(ASTRO_RT) install $(EXPORTS); \
	fi
	@$(MAKE) -C kflc install $(EXPORTS)

uninstall:
	@echo "Manual uninstall: rm files listed in packaging/debian/*/install"
	@echo "(per-lib Makefiles do not provide reverse-install targets)"

# ----- Packaging ---------------------------------------------------

deb:
	@command -v dpkg-buildpackage >/dev/null 2>&1 || \
	    { echo "FAIL: install build deps:"; \
	      echo "  sudo apt install debhelper devscripts dpkg-dev fakeroot"; exit 1; }
	@command -v lintian >/dev/null 2>&1 || \
	    echo "warning: lintian not installed; install with: sudo apt install lintian"
	dpkg-buildpackage -uc -us -b -rfakeroot
	@mkdir -p dist/debs
	@mv ../kfl-stack_*_amd64.deb \
	    ../kfl-stack-core_*_all.deb \
	    ../kfl-stack-aux_*_all.deb \
	    ../kfl-stack-astro_*_all.deb \
	    ../kflc_*_amd64.deb \
	    ../libk26*_*_amd64.deb \
	    dist/debs/ 2>/dev/null || true
	@ls dist/debs/*.deb 2>/dev/null | wc -l | xargs -I{} echo "Built {} .deb files in dist/debs/"

deb-lintian: deb
	lintian -i dist/debs/*.deb || true

deb-clean:
	rm -rf build dist/debs debian/.debhelper debian/tmp debian/files \
	       debian/*.substvars debian/*.debhelper.log \
	       $(addprefix ../kfl-stack,_*.deb _*.changes _*.buildinfo) \
	       $(addprefix ../kflc,_*.deb) \
	       ../libk26*_*.deb ../libk26*_*.changes ../libk26*_*.buildinfo

rpm:
	@echo "rpm packaging scaffold; see packaging/rpm/README.md"
	@exit 1

brew:
	@echo "brew formula scaffold; see packaging/brew/README.md"
	@exit 1

dist:
	@mkdir -p dist
	@if [ -d .git ]; then \
	    git archive --format=tar.gz \
	        --prefix=KFL_Stack-$(VERSION)/ \
	        -o dist/KFL_Stack-$(VERSION).tar.gz HEAD; \
	    echo "wrote dist/KFL_Stack-$(VERSION).tar.gz"; \
	else \
	    echo "no .git here; cannot 'git archive'"; exit 1; \
	fi

# ----- Clean -------------------------------------------------------

clean:
	@for l in $(ALL_C_LIBS) $(ASTRO_FORTRAN) kflc; do \
	    [ -d $$l ] && $(MAKE) -C $$l clean 2>/dev/null || true; \
	done
	@if [ -d tests ]; then $(MAKE) -C tests clean 2>/dev/null || true; fi

# Wipe host-built .a / .o / .so / .mod across every library. Run this
# before a fresh or cross build so host toolchain artifacts (a given
# libc / gcc version) do not leak into a differently-configured build
# and trigger LTO bytecode version mismatches.
clean-host-artifacts:
	@for lib in $(ALL_C_LIBS) $(ASTRO_FORTRAN) $(ASTRO_RT); do \
	    rm -f $$lib/*.a $$lib/*.so $$lib/*.so.* \
	          $$lib/src/*.o $$lib/src/*.mod 2>/dev/null; \
	done
	-rm -f libk26astro_quad/src/upstream/quadpack/*.o \
	       libk26astro_quad/src/upstream/quadpack/*.mod
	-rm -f libk26astro_ode/src/upstream/odepack/*.o \
	       libk26astro_ode/src/upstream/odepack/*.mod
	-rm -f libk26astro_fit/src/upstream/cminpack/*.o
	-rm -f libk26astro_geomag/src/upstream/igrf/*.o
	-rm -f libk26astro_atmos/src/upstream/nrlmsise00/*.o
	-rm -f kflc/bin/kflc kflc/libkflc.a kflc/src/*.o

distclean: clean clean-host-artifacts
	rm -rf dist build

# ----- Version inspection ------------------------------------------

kfl-version:
	@echo "KFL_Stack     : $(VERSION)"
	@echo "KFL grammar   : $(KFL_GRAMMAR_VER)"
	@echo "git SHA       : $(GIT_SHA)"
	@if [ -x kflc/bin/kflc ]; then \
	    printf "kflc          : built\n"; \
	else \
	    printf "kflc          : (not built — run 'make compiler')\n"; \
	fi

kfl-grammar-version:
	@cat KFL_GRAMMAR_VERSION

check-gfortran:
	@command -v gfortran >/dev/null 2>&1 \
	    && echo "gfortran: $$(gfortran --version | head -1)" \
	    || { echo "gfortran missing — Fortran libs will be skipped"; exit 1; }
