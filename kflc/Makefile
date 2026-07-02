CC      ?= cc
PREFIX  ?= /usr
DESTDIR ?=

CFLAGS  ?= -O2 -g
override CFLAGS += -std=c11 -Wall -Wextra -Wpedantic \
                   -Iinclude -Isrc

BIN      = bin/kflc
LIB      = libkflc.a
# Library objects: parser + lexer + dump + arena/diag util. Linked by
# both the kflc binary and any embedding tool (k26form) that needs to
# load .kfl files.
LIB_OBJS = src/util.o src/lexer.o src/parser.o src/dump.o src/serialize.o \
           src/expr.o src/stmt.o src/opaque.o src/builtin_register_astro.o \
           src/builtin_loader.o src/runtime.o
# Driver objects: emitter + CLI driver. Only the kflc binary needs these.
DRV_OBJS = src/emit.o src/main.o

.PHONY: all install clean test
all: $(BIN) $(LIB)

# Round-trip property test. Builds tests/round_trip against libkflc.a
# and runs it over every example .kfl, asserting parse / serialize /
# parse yields a structurally identical AST (the §3.3 semantic-
# equivalence promise). Standalone shell wrapper at tests/run.sh.
TEST_BIN = tests/round_trip
TEST_EXAMPLES = $(sort $(wildcard examples/*.kfl))

$(TEST_BIN): tests/round_trip.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

MANIFEST_TEST = tests/test_manifest_loader
$(MANIFEST_TEST): tests/test_manifest_loader.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

# Emit-side test for parse-time body bindings. Compiles a fixture
# via `kflc --emit` and greps the output for `_kfl_body_<NAME>_idx`
# declarations + idx-deref uses, and the absence of runtime
# `find_body` calls for known names.
BODY_BINDING_TEST = tests/test_body_binding_emit
$(BODY_BINDING_TEST): tests/test_body_binding_emit.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

# Emit-side test for the KFL memory model. Verifies form-level
# arena declarations + init factory + fn-body allocator binding +
# conditional kflc_runtime.h include + slim include surface for
# forms that do not use the memory keywords.
W8_ARENA_TEST = tests/test_w8_arena_emit
$(W8_ARENA_TEST): tests/test_w8_arena_emit.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

# Gate for compile-time enforcement of own / borrow / ptr / arena
# qualifiers. Writes small .kfl fixtures to /tmp, runs the built
# kflc binary against them, asserts that positive cases emit the
# expected C++ shape and negative cases trigger the right diag error.
MEMORY_QUAL_TEST = tests/test_memory_qualifiers
$(MEMORY_QUAL_TEST): tests/test_memory_qualifiers.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

test: $(TEST_BIN) $(MANIFEST_TEST) $(BODY_BINDING_TEST) $(W8_ARENA_TEST) $(MEMORY_QUAL_TEST) $(BIN)
	./$(TEST_BIN) $(TEST_EXAMPLES)
	./$(TEST_BIN) integration_tests/astro_w8.kfl
	./$(MANIFEST_TEST)
	./$(BODY_BINDING_TEST)
	./$(W8_ARENA_TEST)
	./$(MEMORY_QUAL_TEST)

$(LIB): $(LIB_OBJS)
	$(AR) rcs $@ $^

$(BIN): $(DRV_OBJS) $(LIB)
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(DRV_OBJS) $(LIB) $(LDFLAGS)

src/%.o: src/%.c include/kflc.h src/internal.h
	$(CC) $(CFLAGS) -c -o $@ $<

# dump.c and serialize.c only need the public header.
src/dump.o: src/dump.c include/kflc.h
	$(CC) $(CFLAGS) -c -o $@ $<

src/serialize.o: src/serialize.c include/kflc.h
	$(CC) $(CFLAGS) -c -o $@ $<

AR ?= ar

install: $(BIN) $(LIB)
	install -Dm755 $(BIN)                $(DESTDIR)$(PREFIX)/bin/kflc
	install -Dm644 $(LIB)                $(DESTDIR)$(PREFIX)/lib/$(LIB)
	install -Dm644 include/kflc.h        $(DESTDIR)$(PREFIX)/include/kflc.h
	install -Dm644 include/kflc_runtime.h \
	    $(DESTDIR)$(PREFIX)/include/kflc_runtime.h
	install -Dm644 GRAMMAR.md            $(DESTDIR)$(PREFIX)/share/doc/kflc/GRAMMAR.md
	install -Dm644 share/templates/template.kfl \
	    $(DESTDIR)$(PREFIX)/share/doc/kflc/templates/template.kfl
	@for f in $(wildcard examples/*.kfl); do \
	    install -Dm644 $$f \
	        $(DESTDIR)$(PREFIX)/share/doc/kflc/examples/$$(basename $$f); \
	done

clean:
	rm -f src/*.o $(BIN) $(LIB) $(TEST_BIN) $(MANIFEST_TEST) \
	      $(BODY_BINDING_TEST) $(W8_ARENA_TEST) $(MEMORY_QUAL_TEST) \
	      tests/bench_parse
	rm -rf bin
