# tests/astro/Makefile - astro test gates.
#
# Cross-cutting tests that exercise the full astro stack (core +
# body + ephem + conics + grav) end-to-end.

CC      ?= gcc
CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE
CFLAGS  += -ffp-contract=off -fexcess-precision=standard

INC = -I../../libk26astro_core/include \
      -I../../libk26astro_body/include \
      -I../../libk26astro_ephem/include \
      -I../../libk26astro_conics/include \
      -I../../libk26astro_grav/include \
      -I../../libk26astro_atmos/include \
      -I../../libk26astro_vehicle/include \
      -I../../libk26astro_rt/include \
      -I../../libk26tick/include \
      -I../../libk26compute/include \
      -I../../libk26m3d/include

LIBS = ../../libk26astro_rt/libk26astro_rt.a \
       ../../libk26astro_vehicle/libk26astro_vehicle.a \
       ../../libk26astro_atmos/libk26astro_atmos.a \
       ../../libk26astro_grav/libk26astro_grav.a \
       ../../libk26astro_conics/libk26astro_conics.a \
       ../../libk26astro_body/libk26astro_body.a \
       ../../libk26astro_ephem/libk26astro_ephem.a \
       ../../libk26astro_core/libk26astro_core.a \
       ../../libk26compute/libk26compute.a \
       ../../libk26tick/libk26tick.a \
       ../../libk26m3d/libk26m3d.a

LDLIBS = -lgfortran -lm

BIN_DIR = bin

# Gates: two_body_kepler, two_body_wh, burrau_ias15, mercury_gr_ppn1,
# inner_ss_100yr, coord_stress, frame_roundtrip, snapshot_roundtrip,
# test_referenced_replay, test_outer_ss_100jup, test_kozai_lidov,
# test_sj_s_10kyr. The light_time test lives in libk26astro_ephem/tests;
# the KFL snapshot round-trip lives in kflc/Makefile.
TESTS = $(BIN_DIR)/two_body_kepler $(BIN_DIR)/two_body_wh \
        $(BIN_DIR)/burrau_ias15 $(BIN_DIR)/mercury_gr_ppn1 \
        $(BIN_DIR)/inner_ss_100yr $(BIN_DIR)/coord_stress \
        $(BIN_DIR)/frame_roundtrip $(BIN_DIR)/snapshot_roundtrip \
        $(BIN_DIR)/test_referenced_replay \
        $(BIN_DIR)/test_outer_ss_100jup $(BIN_DIR)/test_kozai_lidov \
        $(BIN_DIR)/test_sj_s_10kyr

.PHONY: all run clean
all: $(TESTS)
	@$(MAKE) run

run: $(TESTS)
	@printf "\n==== running astro test gates ====\n"
	@passed=0; skipped=0; \
	for t in $(TESTS); do \
	    printf "  %-30s" "$$(basename $$t)"; \
	    $$t > $$t.log 2>&1; rc=$$?; \
	    if [ $$rc -eq 0 ]; then \
	        printf "  PASS\n"; cat $$t.log | sed 's/^/    /'; \
	        passed=$$((passed+1)); \
	    elif [ $$rc -eq 77 ]; then \
	        printf "  SKIP\n"; cat $$t.log | sed 's/^/    /'; \
	        skipped=$$((skipped+1)); \
	    else \
	        printf "  FAIL (rc=$$rc)\n"; cat $$t.log; exit 1; \
	    fi; \
	done; \
	printf "==== %d passed, %d skipped of %d gate(s) ====\n" \
	    $$passed $$skipped $(words $(TESTS))

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/two_body_kepler: two_body_kepler.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/two_body_wh: two_body_wh.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/burrau_ias15: burrau_ias15.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/mercury_gr_ppn1: mercury_gr_ppn1.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/inner_ss_100yr: inner_ss_100yr.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/coord_stress: coord_stress.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/frame_roundtrip: frame_roundtrip.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/snapshot_roundtrip: snapshot_roundtrip.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/test_referenced_replay: test_referenced_replay.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/test_outer_ss_100jup: test_outer_ss_100jup.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/test_kozai_lidov: test_kozai_lidov.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

$(BIN_DIR)/test_sj_s_10kyr: test_sj_s_10kyr.c $(LIBS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS) $(LDLIBS)

clean:
	rm -rf $(BIN_DIR)
