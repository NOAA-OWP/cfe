#!/bin/bash
#
# bisect_regression.sh — Git bisect helper for CFE regression testing
#
# Builds CFE, runs a regression test, and exits with the comparison result.
# Designed to be called by: git bisect run bash test/scripts/bisect_regression.sh [args]
#
# Usage:
#   git bisect start <bad> <good>
#   git bisect run bash test/scripts/bisect_regression.sh \
#       --test-dir /path/to/test --ngen /path/to/ngen [--tol 1e-6]
#
# Options:
#   --test-dir DIR   Regression test directory (contains realization.json, reference/)
#   --ngen PATH      Path to ngen binary
#   --tol FLOAT      Comparison tolerance (default: 1e-6)
#   --build-dir DIR  CFE build directory (default: ./build)
#   --jobs N          Parallel build jobs (default: 4)
#
# Exit codes (per git-bisect convention):
#   0   = good commit (test passes)
#   1   = bad commit (test fails)
#   125 = skip (build failed, untestable)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

TEST_DIR=""
NGEN_BIN=""
TOL="1e-6"
BUILD_DIR="$CFE_ROOT/build"
JOBS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test-dir)  TEST_DIR="$2"; shift 2 ;;
        --ngen)      NGEN_BIN="$2"; shift 2 ;;
        --tol)       TOL="$2"; shift 2 ;;
        --build-dir) BUILD_DIR="$2"; shift 2 ;;
        --jobs)      JOBS="$2"; shift 2 ;;
        *)           echo "Unknown option: $1" >&2; exit 125 ;;
    esac
done

[ -n "$TEST_DIR" ] || { echo "ERROR: --test-dir required" >&2; exit 125; }
[ -n "$NGEN_BIN" ] || { echo "ERROR: --ngen required" >&2; exit 125; }

COMMIT="$(git -C "$CFE_ROOT" rev-parse --short HEAD)"
echo "=== Bisect testing commit $COMMIT ==="

# Build CFE
echo "Building CFE..."
if ! cmake --build "$BUILD_DIR" --target cfebmi -j"$JOBS" 2>&1 | tail -3; then
    echo "Build failed — skipping commit $COMMIT"
    exit 125
fi

# Run regression test
echo "Running regression test..."
bash "$SCRIPT_DIR/run_ngen_regression.sh" "$TEST_DIR" --ngen "$NGEN_BIN" --tol "$TOL"
RESULT=$?

if [ $RESULT -eq 0 ]; then
    echo "=== Commit $COMMIT: GOOD ==="
else
    echo "=== Commit $COMMIT: BAD ==="
    RESULT=1
fi

exit $RESULT
