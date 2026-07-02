CC      ?= gcc
AR      ?= ar

# Geographic primitives. Depends on libk26m3d (math types), libexpat
# (XML SAX), libz, libm. cJSON not used here.

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -I../libk26m3d/include
CFLAGS  += -ffp-contract=off -fexcess-precision=standard
CFLAGS  += $(shell pkg-config --cflags expat zlib 2>/dev/null)
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26geo.a

PUBLIC_HDRS  = $(wildcard include/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.0
LIB_BASE     = libk26geo
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = $(shell pkg-config --libs expat zlib 2>/dev/null) -lm

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
	install -d $(DESTDIR)$(PREFIX)/include
	install -m755 $(SHLIB_FULL) $(DESTDIR)$(LIBDIR)/$(SHLIB_FULL)
	ln -sf $(SHLIB_FULL)   $(DESTDIR)$(LIBDIR)/$(SHLIB_SONAME)
	ln -sf $(SHLIB_SONAME) $(DESTDIR)$(LIBDIR)/$(SHLIB_LINK)
	install -m644 $(LIB) $(DESTDIR)$(LIBDIR)/$(LIB)
	for h in $(PUBLIC_HDRS); do install -m644 $$h $(DESTDIR)$(PREFIX)/include/; done

TEST_LINK = $(LIB) ../libk26m3d/libk26m3d.a -lexpat -lz -lm
TESTS = tests/test_mercator tests/test_rdp tests/test_osmtags tests/test_spgrid tests/test_osm_xml

test: $(LIB) $(TESTS:=.c)
	$(CC) $(CFLAGS) -o tests/test_mercator tests/test_mercator.c $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_rdp      tests/test_rdp.c      $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_osmtags  tests/test_osmtags.c  $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_spgrid   tests/test_spgrid.c   $(TEST_LINK)
	$(CC) $(CFLAGS) -o tests/test_osm_xml  tests/test_osm_xml.c  $(TEST_LINK)
	./tests/test_mercator && ./tests/test_rdp && ./tests/test_osmtags && \
	    ./tests/test_spgrid && ./tests/test_osm_xml

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) $(TESTS)
