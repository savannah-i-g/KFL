CC      ?= gcc
AR      ?= ar

# Body taxonomy, orbital elements, attitude state, IAU rotation models,
# Earth orientation (GCRS↔ECEF full IAU 2006/2000A pipeline), star
# catalogue ingest.
#
# Depends on: libk26astro_core, libk26astro_ephem, libk26m3d.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26astro_core/include
CFLAGS  += -I../libk26astro_ephem/include
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26astro_body.a

PUBLIC_HDRS  = $(wildcard include/k26astro_body/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.0
LIB_BASE     = libk26astro_body
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
	install -d $(DESTDIR)$(PREFIX)/include/k26astro_body
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	for h in $(PUBLIC_HDRS); do \
	    install -m644 $$h $(DESTDIR)$(PREFIX)/include/k26astro_body/; \
	done
	install -d $(DESTDIR)$(PREFIX)/share/kflc/builtins
	install -m644 kfl/k26astro_body.kflbi \
	    $(DESTDIR)$(PREFIX)/share/kflc/builtins/k26astro_body.kflbi

TEST_LINK = $(LIB) \
            ../libk26astro_ephem/libk26astro_ephem.a \
            ../libk26astro_core/libk26astro_core.a \
            ../libk26m3d/libk26m3d.a -lm

test: $(LIB) tests/test_elements.c tests/test_attitude.c tests/test_rotation.c tests/test_earth.c tests/test_star.c tests/test_attitude_torque_registry.c tests/test_attitude_step_ext.c tests/test_body_set_mass.c
	$(CC) $(CFLAGS) -o tests/test_elements tests/test_elements.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_attitude tests/test_attitude.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_rotation tests/test_rotation.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_earth    tests/test_earth.c    $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_star     tests/test_star.c     $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_attitude_torque_registry \
	    tests/test_attitude_torque_registry.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_attitude_step_ext \
	    tests/test_attitude_step_ext.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_body_set_mass \
	    tests/test_body_set_mass.c $(TEST_LINK)
	./tests/test_elements && ./tests/test_attitude && ./tests/test_rotation && \
	    ./tests/test_earth && ./tests/test_star && \
	    ./tests/test_attitude_torque_registry && \
	    ./tests/test_attitude_step_ext && \
	    ./tests/test_body_set_mass

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) tests/test_elements tests/test_attitude \
	    tests/test_rotation tests/test_earth tests/test_star \
	    tests/test_attitude_torque_registry tests/test_attitude_step_ext \
	    tests/test_body_set_mass
