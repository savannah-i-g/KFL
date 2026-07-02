CC      ?= gcc
AR      ?= ar

# Runtime manager. Composes K26AstroGravState + ephem/frame/conics
# into K26AstroWorld. Multi-rate scheduler, MERCURIUS hybrid, observer
# light-time/aberration, versioned snapshot format.
#
# Depends on: libk26astro_{core,body,ephem,conics,grav,atmos},
#             libk26tick, libk26m3d, libk26compute.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26astro_core/include
CFLAGS  += -I../libk26astro_body/include
CFLAGS  += -I../libk26astro_ephem/include
CFLAGS  += -I../libk26astro_conics/include
CFLAGS  += -I../libk26astro_grav/include
CFLAGS  += -I../libk26astro_atmos/include
CFLAGS  += -I../libk26astro_vehicle/include
CFLAGS  += -I../libk26tick/include
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard -frounding-math
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26astro_rt.a

PUBLIC_HDRS  = $(wildcard include/k26astro_rt/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.0
LIB_BASE     = libk26astro_rt
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lgfortran -lm

.PHONY: all install clean test
all: $(LIB) $(SHLIB_FULL)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(SHLIB_FULL): $(OBJ)
	$(CC) -shared -Wl,-soname,$(SHLIB_SONAME) -o $@ $^ $(SHLIB_LDLIBS)

src/%.o: src/%.c $(PUBLIC_HDRS) $(PRIVATE_HDRS)
	$(CC) $(CFLAGS) -c -o $@ $<

install: $(LIB) $(SHLIB_FULL)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_rt
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	for h in $(PUBLIC_HDRS); do \
	    install -m644 $$h $(DESTDIR)$(PREFIX)/include/k26astro_rt/; \
	done
	install -d $(DESTDIR)$(PREFIX)/share/doc/libk26astro_rt/examples
	for e in $(wildcard examples/*.kfl); do \
	    install -m644 $$e $(DESTDIR)$(PREFIX)/share/doc/libk26astro_rt/examples/; \
	done
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -m644 kfl/k26astro_rt.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/k26astro_rt.kflbi

TEST_LINK = $(LIB) \
            ../libk26astro_vehicle/libk26astro_vehicle.a \
            ../libk26astro_grav/libk26astro_grav.a \
            ../libk26astro_conics/libk26astro_conics.a \
            ../libk26astro_atmos/libk26astro_atmos.a \
            ../libk26astro_body/libk26astro_body.a \
            ../libk26astro_ephem/libk26astro_ephem.a \
            ../libk26astro_core/libk26astro_core.a \
            ../libk26compute/libk26compute.a \
            ../libk26tick/libk26tick.a \
            ../libk26m3d/libk26m3d.a -lgfortran -lm

TEST_SRC = $(wildcard tests/test_*.c)
TEST_BIN = $(TEST_SRC:.c=)

test: $(LIB) $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

tests/test_%: tests/test_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(TEST_LINK)

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) $(TEST_BIN)
