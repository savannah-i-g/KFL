# KFL_Stack/tests/Makefile — top-level test harness for carved libs.
#
# Only subdirs whose parent libraries are part of the KFL_Stack carve
# live here: astro/ (12-lib astro gate suite, has its own Makefile)
# and tick/ (single-file libk26tick test). Other K26 test subdirs
# (cam/, fx/, model/, pylon/) live with their parent
# libs.

CC       ?= gcc
CFLAGS   ?= -std=c11 -O0 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE

INC      = -I../libk26m3d/include -I../libk26compute/include \
           -I../libk26tick/include

LIB_TICK = ../libk26tick/libk26tick.a

BIN_DIR  = bin
TESTS    = $(BIN_DIR)/test_tick

.PHONY: all astro tick clean
all: tick astro

# tick: build + run
tick: $(BIN_DIR)/test_tick
	@printf "  %-25s" "test_tick"
	@if ./$(BIN_DIR)/test_tick > $(BIN_DIR)/test_tick.log 2>&1; then \
	    printf "  PASS\n"; \
	else \
	    printf "  FAIL\n"; cat $(BIN_DIR)/test_tick.log; exit 1; \
	fi

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/test_tick: tick/test_tick.c $(LIB_TICK) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIB_TICK) -lm

# astro: delegate to its own Makefile
astro:
	$(MAKE) -C astro

clean:
	rm -rf $(BIN_DIR)
	$(MAKE) -C astro clean 2>/dev/null || true
