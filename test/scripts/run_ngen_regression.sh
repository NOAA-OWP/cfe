#!/bin/bash
#
# run_ngen_regression.sh — Run a single ngen regression test and compare output
#
# Usage:
#   run_ngen_regression.sh <test_dir> [options]
#
# Arguments:
#   test_dir    Test directory containing realization.json, catchment_data.geojson,
#               nexus_data.geojson, and reference/ subdirectory
#
# Options:
#   --ngen PATH         Path to ngen binary (default: auto-detect via $NGEN or ./cmake_build/ngen)
#   --tol FLOAT         Absolute tolerance for comparison (default: 1e-6)
#   --baseline          Generate baseline: run ngen and save output as new reference
#   --var COLUMN        Output variable column index, 0-based (default: 2 = discharge)
#   --ref FILE          Reference file name within reference/ (default: auto-detect)
#
# The script expects the test directory layout:
#   test_dir/
#     realization.json
#     catchment_data.geojson
#     nexus_data.geojson
#     reference/           <- reference output files
#     output/              <- created/cleared by this script
#
# Exit: 0 = PASS, 1 = FAIL, 2 = setup error

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE="$SCRIPT_DIR/compare_discharge.py"

# Defaults
NGEN_BIN="${NGEN:-}"
TOL="1e-6"
BASELINE=false
REF_FILE=""

# Parse args
TEST_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ngen)     NGEN_BIN="$2"; shift 2 ;;
        --tol)      TOL="$2"; shift 2 ;;
        --baseline) BASELINE=true; shift ;;
        --ref)      REF_FILE="$2"; shift 2 ;;
        -*)         echo "Unknown option: $1" >&2; exit 2 ;;
        *)          TEST_DIR="$1"; shift ;;
    esac
done

if [ -z "$TEST_DIR" ]; then
    echo "Usage: run_ngen_regression.sh <test_dir> [--ngen PATH] [--tol FLOAT] [--baseline]" >&2
    exit 2
fi

TEST_DIR="$(cd "$TEST_DIR" && pwd)"

# Find ngen
if [ -z "$NGEN_BIN" ]; then
    # Walk up looking for cmake_build/ngen
    d="$TEST_DIR"
    while [ "$d" != "/" ]; do
        if [ -x "$d/cmake_build/ngen" ]; then
            NGEN_BIN="$d/cmake_build/ngen"
            break
        fi
        d="$(dirname "$d")"
    done
fi

[ -x "$NGEN_BIN" ] || { echo "ERROR: ngen not found. Set --ngen or \$NGEN." >&2; exit 2; }

# Validate test directory
for f in realization.json catchment_data.geojson nexus_data.geojson; do
    [ -f "$TEST_DIR/$f" ] || { echo "ERROR: missing $TEST_DIR/$f" >&2; exit 2; }
done

# Determine ngen root (realization paths are relative to cwd)
# Walk up from test_dir to find a directory that contains the paths used in realization
NGEN_ROOT="$TEST_DIR"
while [ "$NGEN_ROOT" != "/" ]; do
    if [ -d "$NGEN_ROOT/extern" ] || [ -d "$NGEN_ROOT/cmake_build" ]; then
        break
    fi
    NGEN_ROOT="$(dirname "$NGEN_ROOT")"
done

TEST_NAME="$(basename "$TEST_DIR")"
echo "============================================================"
echo "Test: $TEST_NAME"
echo "ngen: $NGEN_BIN"
echo "dir:  $TEST_DIR"
echo "============================================================"

# Clear and create output directory
mkdir -p "$TEST_DIR/output"
rm -f "$TEST_DIR/output/"*.csv 2>/dev/null || true

# Run ngen from the ngen root so relative paths in realization.json resolve
cd "$NGEN_ROOT"

# Compute paths relative to ngen root
REL_TEST="$(python3 -c "import os; print(os.path.relpath('$TEST_DIR', '$NGEN_ROOT'))")"

"$NGEN_BIN" \
    "$REL_TEST/catchment_data.geojson" "" \
    "$REL_TEST/nexus_data.geojson" "" \
    "$REL_TEST/realization.json" 2>&1 | tail -5

# Find the output file (first CSV in output/)
OUTPUT_CSV="$(find "$TEST_DIR/output" -name '*.csv' -not -name 'nex-*' | head -1)"

if [ -z "$OUTPUT_CSV" ] || [ ! -f "$OUTPUT_CSV" ]; then
    echo "FAIL: no output CSV produced in $TEST_DIR/output/"
    exit 1
fi

echo "Output: $OUTPUT_CSV"

# Baseline mode: save output as reference
if $BASELINE; then
    mkdir -p "$TEST_DIR/reference"
    BASELINE_FILE="$TEST_DIR/reference/q_baseline.csv"
    cp "$OUTPUT_CSV" "$BASELINE_FILE"
    echo "Baseline saved: $BASELINE_FILE"
    echo "PASS (baseline generated)"
    exit 0
fi

# Find reference file
if [ -n "$REF_FILE" ]; then
    REF_PATH="$TEST_DIR/reference/$REF_FILE"
elif [ -f "$TEST_DIR/reference/q_baseline.csv" ]; then
    REF_PATH="$TEST_DIR/reference/q_baseline.csv"
elif [ -f "$TEST_DIR/reference/q.csv" ]; then
    REF_PATH="$TEST_DIR/reference/q.csv"
elif [ -f "$TEST_DIR/reference/q.out" ]; then
    REF_PATH="$TEST_DIR/reference/q.out"
else
    echo "FAIL: no reference file found in $TEST_DIR/reference/"
    exit 1
fi

echo "Reference: $REF_PATH"
echo ""

# Compare
python3 "$COMPARE" "$REF_PATH" "$OUTPUT_CSV" --tol "$TOL"
