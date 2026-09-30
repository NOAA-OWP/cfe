#!/bin/bash
#
# generate_standalone_baseline.sh — Generate reference output from cfe_bmi_driver
#
# Runs the standalone CFE BMI driver to produce full-precision reference output
# for regression testing. This output has no formatting truncation (unlike Fred's
# %.8f output), so ngen-vs-standalone comparisons can achieve exact match.
#
# Usage:
#   generate_standalone_baseline.sh --driver PATH --config FILE --forcing FILE --outdir DIR
#
# Options:
#   --driver PATH    Path to cfe_bmi_driver executable
#   --config FILE    CFE config file (standalone version, not BMI)
#   --forcing FILE   AORC forcing CSV
#   --outdir DIR     Output directory for reference files (created if needed)
#   --verbose N      Verbosity level (default: 0)
#
# Produces: q.out, fluxes.out, storage.out, thetas.out in outdir

set -euo pipefail

DRIVER=""
CONFIG=""
FORCING=""
OUTDIR=""
VERBOSE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --driver)   DRIVER="$2"; shift 2 ;;
        --config)   CONFIG="$2"; shift 2 ;;
        --forcing)  FORCING="$2"; shift 2 ;;
        --outdir)   OUTDIR="$2"; shift 2 ;;
        --verbose)  VERBOSE="$2"; shift 2 ;;
        *)          echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[ -x "$DRIVER" ]  || { echo "ERROR: driver not found or not executable: $DRIVER" >&2; exit 2; }
[ -f "$CONFIG" ]   || { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }
[ -f "$FORCING" ]  || { echo "ERROR: forcing not found: $FORCING" >&2; exit 2; }
[ -n "$OUTDIR" ]   || { echo "ERROR: --outdir required" >&2; exit 2; }

mkdir -p "$OUTDIR"

echo "Generating standalone baseline..."
echo "  Driver:  $DRIVER"
echo "  Config:  $CONFIG"
echo "  Forcing: $FORCING"
echo "  Output:  $OUTDIR"

"$DRIVER" \
    -c "$CONFIG" \
    -f "$FORCING" \
    -q "$OUTDIR/q.out" \
    -x "$OUTDIR/fluxes.out" \
    -s "$OUTDIR/storage.out" \
    -t "$OUTDIR/thetas.out" \
    -v "$VERBOSE"

echo ""
echo "Baseline files generated:"
for f in q.out fluxes.out storage.out thetas.out; do
    if [ -f "$OUTDIR/$f" ]; then
        lines=$(wc -l < "$OUTDIR/$f")
        echo "  $f: $lines lines"
    fi
done
echo "Done."
