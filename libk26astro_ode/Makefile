CC      ?= gcc
FC      = gfortran
AR      ?= ar

# Initial-value ODE solver via ODEPACK / DLSODA (jacobwilliams modernized
# fork; vendored under src/upstream/odepack/).

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude
CFLAGS  += -ffp-contract=off -fexcess-precision=standard -frounding-math
CFLAGS  += -fPIC

FFLAGS  ?= -O2 -g -fimplicit-none
FFLAGS  += -ffp-contract=off
FFLAGS  += -ffree-line-length-none
FFLAGS  += -Wall
FFLAGS  += -fPIC

UPSTREAM_FFLAGS = $(FFLAGS) \
                  -fallow-argument-mismatch \
                  -Wno-unused-dummy-argument \
                  -Wno-maybe-uninitialized \
                  -Wno-uninitialized \
                  -Wno-integer-division \
                  -Wno-conversion \
                  -Wno-compare-reals \
                  -Wno-character-truncation \
                  -Wno-unused-variable \
                  -Wno-unused-label \
                  -Wno-unused-function

C_SRC = $(wildcard src/*.c)
C_OBJ = $(C_SRC:.c=.o)

F_UPSTREAM = src/upstream/odepack/M_odepack.f90
F_WRAPPER  = src/k26astro_ode_iface.f90
F_OBJ      = $(F_UPSTREAM:.f90=.o) $(F_WRAPPER:.f90=.o)

LIB = libk26astro_ode.a

PUBLIC_HDRS = $(wildcard include/k26astro_ode/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=


SONAME_MAJOR = 0
LIB_VER      = 0.4.0
LIB_BASE     = libk26astro_ode
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lgfortran -lm

TEST_SRC  = $(wildcard tests/test_*.c)
TEST_BIN  = $(TEST_SRC:.c=)
TEST_LINK = $(LIB) -lgfortran -lm

.PHONY: all test install clean
all: $(LIB) $(SHLIB_FULL)

$(LIB): $(C_OBJ) $(F_OBJ)
	$(AR) rcs $@ $^

$(SHLIB_FULL): $(C_OBJ) $(F_OBJ)
	$(CC) -shared -Wl,-soname,$(SHLIB_SONAME) -Wl,--version-script=../k26astro_exports.map -o $@ $^ $(SHLIB_LDLIBS)

src/%.o: src/%.c $(PUBLIC_HDRS)
	$(CC) $(CFLAGS) -c -o $@ $<

src/k26astro_ode_iface.o: src/k26astro_ode_iface.f90 \
                          src/upstream/odepack/M_odepack.o
	$(FC) $(FFLAGS) -Isrc/upstream/odepack -Jsrc -c -o $@ $<

src/upstream/odepack/M_odepack.o: \
        src/upstream/odepack/M_odepack.f90
	$(FC) $(UPSTREAM_FFLAGS) \
	    -Jsrc/upstream/odepack \
	    -c -o $@ $<

test: $(LIB) $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

tests/test_%: tests/test_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_LINK)

install: $(LIB) $(SHLIB_FULL)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_ode
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -d $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_ode/odepack
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	install -m644 $(PUBLIC_HDRS) \
	    $(DESTDIR)$(PREFIX)/include/k26astro_ode/
	install -m644 kfl/k26astro_ode.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/
	install -m644 src/upstream/odepack/LICENSE \
	    $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_ode/odepack/LICENSE
	install -m644 src/upstream/odepack/UPSTREAM.md \
	    $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_ode/odepack/UPSTREAM.md

clean:
	rm -f src/*.o src/*.mod
	rm -f src/upstream/odepack/*.o src/upstream/odepack/*.mod
	rm -f $(LIB) $(SHLIB_FULL) $(TEST_BIN)
