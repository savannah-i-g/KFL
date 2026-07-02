CC      ?= gcc
AR      ?= ar

# Substrate for the K26 astro suite. Provides: sector-grid position
# (K26AstroPos) with Q64.64 determinism fallback, time-scale arithmetic
# (K26AstroEpoch + TAI/UTC/UT1/TT/TDB conversions), inertial frame
# registry + transforms, SI unit conversions, IAU 2015 / CODATA 2018
# constant catalogue.
#
# Depends on libm only; the determinism-mode contract requires
# -ffp-contract=off + -fexcess-precision=standard on every consumer.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26astro_core.a

PUBLIC_HDRS  = $(wildcard include/k26astro_core/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.1
LIB_BASE     = libk26astro_core
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lm

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
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_core
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	for h in $(PUBLIC_HDRS); do \
	    install -m644 $$h $(DESTDIR)$(PREFIX)/include/k26astro_core/; \
	done
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -m644 kfl/k26astro_core.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/k26astro_core.kflbi

TEST_LINK = $(LIB) ../libk26m3d/libk26m3d.a -lm

test: $(LIB) tests/test_pos.c tests/test_epoch.c tests/test_frame.c tests/test_units.c tests/test_sum.c
	$(CC) $(CFLAGS) -o tests/test_pos    tests/test_pos.c    $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_epoch  tests/test_epoch.c  $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_frame  tests/test_frame.c  $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_units  tests/test_units.c  $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_sum    tests/test_sum.c    $(TEST_LINK)
	./tests/test_pos && ./tests/test_epoch && ./tests/test_frame && \
	    ./tests/test_units && ./tests/test_sum

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) tests/test_pos tests/test_epoch \
	    tests/test_frame tests/test_units tests/test_sum
