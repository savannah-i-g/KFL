CC      ?= gcc
AR      ?= ar

# JPL Development Ephemeris ingest + query.
#
# Depends on: libk26astro_core, libk26m3d.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26astro_core/include
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26astro_ephem.a

PUBLIC_HDRS  = $(wildcard include/k26astro_ephem/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.1
LIB_BASE     = libk26astro_ephem
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lm

.PHONY: all install clean test tools
all: $(LIB) $(SHLIB_FULL)

tools: tools/gen_de441_inner

tools/gen_de441_inner: tools/gen_de441_inner.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) ../libk26astro_core/libk26astro_core.a ../libk26m3d/libk26m3d.a -lm

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(SHLIB_FULL): $(OBJ)
	$(CC) -shared -Wl,-soname,$(SHLIB_SONAME) -o $@ $^ $(SHLIB_LDLIBS)

src/%.o: src/%.c $(PUBLIC_HDRS) $(PRIVATE_HDRS)
	$(CC) $(CFLAGS) -c -o $@ $<

install: $(LIB) $(SHLIB_FULL)
	install -d $(DESTDIR)$(LIBDIR)
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_ephem
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	for h in $(PUBLIC_HDRS); do \
	    install -m644 $$h $(DESTDIR)$(PREFIX)/include/k26astro_ephem/; \
	done
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -m644 kfl/k26astro_ephem.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/k26astro_ephem.kflbi

TEST_LINK = $(LIB) ../libk26astro_core/libk26astro_core.a ../libk26m3d/libk26m3d.a -lm

test: $(LIB) tests/test_cheby.c tests/test_spk.c tests/test_light_time.c tests/test_chain_walker.c tests/test_outer_planets.c
	$(CC) $(CFLAGS) -o tests/test_cheby          tests/test_cheby.c          $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_spk            tests/test_spk.c            $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_light_time     tests/test_light_time.c     $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_chain_walker   tests/test_chain_walker.c   $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_outer_planets  tests/test_outer_planets.c  $(TEST_LINK)
	./tests/test_cheby && ./tests/test_spk && ./tests/test_light_time && \
	    ./tests/test_chain_walker && ./tests/test_outer_planets

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) tests/test_cheby tests/test_spk \
	    tests/test_light_time tests/test_chain_walker tests/test_outer_planets \
	    tests/fixtures/*.spk \
	    tools/gen_de441_inner
