CC      ?= gcc
AR      ?= ar

# Fixed-step accumulator + named rate channels. Pure C, no external
# deps beyond libm. Builds host (glibc) and sidecar (musl).

CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -Iinclude -Isrc
CFLAGS  += -ffp-contract=off -fexcess-precision=standard
CFLAGS  += -fPIC

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
LIB = libk26tick.a

PUBLIC_HDRS  = $(wildcard include/*.h)
PRIVATE_HDRS = $(wildcard src/*.h)

PREFIX  ?= /usr
LIBDIR  ?= $(PREFIX)/lib
DESTDIR ?=

SONAME_MAJOR = 0
LIB_VER      = 0.4.1
LIB_BASE     = libk26tick
SHLIB_FULL   = $(LIB_BASE).so.$(LIB_VER)
SHLIB_SONAME = $(LIB_BASE).so.$(SONAME_MAJOR)
SHLIB_LINK   = $(LIB_BASE).so
SHLIB_LDLIBS = -lm

TEST_SRC = $(wildcard tests/test_*.c)
TEST_BIN = $(TEST_SRC:.c=)

.PHONY: all test install clean
all: $(LIB) $(SHLIB_FULL)

test: $(LIB) $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

tests/test_%: tests/test_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) -lm

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

clean:
	rm -f src/*.o $(LIB) $(SHLIB_FULL) $(TEST_BIN)
