CC      ?= gcc
FC      = gfortran
AR      ?= ar

# Planetary atmospheric model (refraction + single-scatter Rayleigh +
# HG-Mie + NRLMSISE-00 thermosphere). NRL FORTRAN under public domain
# (US Government work).
#
# Depends on: libk26m3d (K26V3) + libk26astro_core (constants).

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26astro_core/include
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard -frounding-math
CFLAGS  += -fPIC

FFLAGS  ?= -O2 -g
FFLAGS  += -ffp-contract=off
FFLAGS  += -fdefault-real-8 -fdefault-double-8
FFLAGS  += -ffree-line-length-none
FFLAGS  += -Wall
FFLAGS  += -fPIC

UPSTREAM_FFLAGS = $(FFLAGS) \
                  -ffixed-form \
                  -Wno-unused-dummy-argument \
                  -Wno-maybe-uninitialized \
                  -Wno-integer-division \
                  -Wno-conversion \
                  -Wno-compare-reals \
                  -Wno-conversion-extra \
                  -std=legacy \
                  -fdec

C_SRC = $(wildcard src/*.c)
C_OBJ = $(C_SRC:.c=.o)

F_UPSTREAM = src/upstream/nrlmsise00/NRLMSISE-00_model.FOR
F_UPSTREAM_OBJ = src/upstream/nrlmsise00/NRLMSISE-00_model.o

F_WRAPPER = src/k26astro_atmos_nrlmsise_iface.f90
F_WRAPPER_OBJ = src/k26astro_atmos_nrlmsise_iface.o

F_OBJ = $(F_UPSTREAM_OBJ) $(F_WRAPPER_OBJ)

LIB = libk26astro_atmos.a

PUBLIC_HDRS  = $(wildcard include/k26astro_atmos/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=


SONAME_MAJOR = 0
LIB_VER      = 0.4.1
LIB_BASE     = libk26astro_atmos
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lgfortran -lm

TEST_SRC = $(wildcard tests/test_*.c)
TEST_BIN = $(TEST_SRC:.c=)
TEST_LINK = $(LIB) \
    ../libk26astro_core/libk26astro_core.a \
    ../libk26m3d/libk26m3d.a \
    -lgfortran -lm

.PHONY: all test install clean
all: $(LIB) $(SHLIB_FULL)

$(LIB): $(C_OBJ) $(F_OBJ)
	$(AR) rcs $@ $^

$(SHLIB_FULL): $(C_OBJ) $(F_OBJ)
	$(CC) -shared -Wl,-soname,$(SHLIB_SONAME) -Wl,--version-script=../k26astro_exports.map -o $@ $^ $(SHLIB_LDLIBS)

src/%.o: src/%.c $(PUBLIC_HDRS)
	$(CC) $(CFLAGS) -c -o $@ $<

src/k26astro_atmos_nrlmsise_iface.o: src/k26astro_atmos_nrlmsise_iface.f90
	$(FC) $(FFLAGS) -Jsrc -c -o $@ $<

src/upstream/nrlmsise00/NRLMSISE-00_model.o: src/upstream/nrlmsise00/NRLMSISE-00_model.FOR
	$(FC) $(UPSTREAM_FFLAGS) -c -o $@ $<

test: $(LIB) $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

tests/test_%: tests/test_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_LINK)

install: $(LIB) $(SHLIB_FULL)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_atmos
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	install -m644 $(PUBLIC_HDRS) \
	    $(DESTDIR)$(PREFIX)/include/k26astro_atmos/
	install -m644 kfl/k26astro_atmos.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/
	install -Dm644 src/upstream/nrlmsise00/LICENSE \
	    $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_atmos/nrlmsise00/LICENSE

clean:
	rm -f src/*.o $(F_UPSTREAM_OBJ) $(F_WRAPPER_OBJ) src/*.mod \
	      $(LIB) $(SHLIB_FULL) $(TEST_BIN)
