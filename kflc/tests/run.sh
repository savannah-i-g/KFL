#!/bin/sh
# tests/run.sh — convenience wrapper for `make test`.
#
# Runs the round-trip property test (BL-6 / WP-B) and prints
# a one-line summary. Exit code mirrors the test binary: 0 if all
# examples round-trip cleanly, 1 otherwise.
#
# Usage:
#   ./tests/run.sh                # run all examples
#   ./tests/run.sh path/to.kfl    # run a single .kfl

set -eu

cd "$(dirname "$0")/.."

# Rebuild the test binary on demand. Skip the full `make all` so we
# don't relink the driver unnecessarily.
make tests/round_trip >/dev/null

if [ "$#" -gt 0 ]; then
    exec ./tests/round_trip "$@"
fi

exec ./tests/round_trip examples/*.kfl
