CC      ?= gcc
AR      ?= ar

# Nonlinear least-squares via CMINPACK (BSD-3, vendored). Pure C, no
# Fortran involvement.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -D_GNU_SOURCE
CFLAGS  += -Iinclude
CFLAGS  += -Isrc/upstream/cminpack
CFLAGS  += -ffp-contract=off -fexcess-precision=standard -frounding-math
CFLAGS  += -fPIC

UPSTREAM_CFLAGS = $(CFLAGS) \
                  -Wno-unused-parameter \
                  -Wno-unused-variable \
                  -Wno-unused-but-set-variable \
                  -Wno-sign-compare

C_SRC = src/k26astro_fit.c
C_OBJ = src/k26astro_fit.o

UPSTREAM_SRC = \
    src/upstream/cminpack/dpmpar.c \
    src/upstream/cminpack/enorm.c \
    src/upstream/cminpack/fdjac2.c \
    src/upstream/cminpack/lmdif.c \
    src/upstream/cminpack/lmdif1.c \
    src/upstream/cminpack/lmder.c \
    src/upstream/cminpack/lmder1.c \
    src/upstream/cminpack/lmpar.c \
    src/upstream/cminpack/qrfac.c \
    src/upstream/cminpack/qrsolv.c

UPSTREAM_OBJ = $(UPSTREAM_SRC:.c=.o)

LIB = libk26astro_fit.a

PUBLIC_HDRS = $(wildcard include/k26astro_fit/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.1
LIB_BASE     = libk26astro_fit
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lm

# Export only the k26astro_* public surface; keep the vendored CMINPACK
# symbols (lmdif, enorm, ...) private so they cannot collide in a
# consumer's link.
VERSION_SCRIPT = ../k26astro_exports.map

TEST_SRC  = $(wildcard tests/test_*.c)
TEST_BIN  = $(TEST_SRC:.c=)
TEST_LINK = $(LIB) -lm

.PHONY: all test install clean
all: $(LIB) $(SHLIB_FULL)

$(LIB): $(C_OBJ) $(UPSTREAM_OBJ)
	$(AR) rcs $@ $^

$(SHLIB_FULL): $(C_OBJ) $(UPSTREAM_OBJ)
	$(CC) -shared -Wl,-soname,$(SHLIB_SONAME) \
	    -Wl,--version-script=$(VERSION_SCRIPT) -o $@ $^ $(SHLIB_LDLIBS)

src/k26astro_fit.o: src/k26astro_fit.c $(PUBLIC_HDRS) \
                   src/upstream/cminpack/cminpack.h
	$(CC) $(CFLAGS) -c -o $@ $<

src/upstream/cminpack/%.o: src/upstream/cminpack/%.c \
                          src/upstream/cminpack/cminpack.h \
                          src/upstream/cminpack/cminpackP.h
	$(CC) $(UPSTREAM_CFLAGS) -c -o $@ $<

test: $(LIB) $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

tests/test_%: tests/test_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_LINK)

install: $(LIB) $(SHLIB_FULL)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_fit
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -d $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_fit/cminpack
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	install -m644 $(PUBLIC_HDRS) \
	    $(DESTDIR)$(PREFIX)/include/k26astro_fit/
	install -m644 kfl/k26astro_fit.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/
	install -m644 src/upstream/cminpack/LICENSE \
	    $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_fit/cminpack/LICENSE
	install -m644 src/upstream/cminpack/UPSTREAM.md \
	    $(DESTDIR)$(PREFIX)/share/licenses/libk26astro_fit/cminpack/UPSTREAM.md

clean:
	rm -f src/*.o
	rm -f src/upstream/cminpack/*.o
	rm -f $(LIB) $(SHLIB_FULL) $(TEST_BIN)
