# k26astro_replay — REFERENCED-mode op-log inspector (host-only).
#
# Compares two logs byte-for-byte (--check) or dumps a human-
# readable view (--dump). The schema constants come from
# libk26astro_rt's referenced.h.

CC      ?= gcc
CFLAGS  ?= -std=c11 -O2 -g -Wall -Wextra -Wpedantic -D_GNU_SOURCE

INC = -I../../libk26astro_rt/include \
      -I../../libk26astro_grav/include \
      -I../../libk26astro_body/include \
      -I../../libk26astro_ephem/include \
      -I../../libk26astro_core/include \
      -I../../libk26astro_conics/include \
      -I../../libk26compute/include \
      -I../../libk26tick/include \
      -I../../libk26m3d/include

LIBS = ../../libk26astro_rt/libk26astro_rt.a \
       ../../libk26astro_grav/libk26astro_grav.a \
       ../../libk26astro_conics/libk26astro_conics.a \
       ../../libk26astro_body/libk26astro_body.a \
       ../../libk26astro_ephem/libk26astro_ephem.a \
       ../../libk26astro_core/libk26astro_core.a \
       ../../libk26compute/libk26compute.a \
       ../../libk26tick/libk26tick.a \
       ../../libk26m3d/libk26m3d.a -lm

.PHONY: all clean
all: k26astro_replay

k26astro_replay: main.c $(LIBS)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS)

clean:
	rm -f k26astro_replay
