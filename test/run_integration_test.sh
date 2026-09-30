#!/bin/bash
# run_integration_test.sh — Run cfe_bmi_driver with a config and compare against golden outputs
# Usage: run_integration_test.sh <driver> <config> <forcing> <golden_dir> [tolerance]
# Golden dir must contain q.out, fluxes.out, storage.out (and optionally thetas.out)
set -e
DRIVER="$1"; CONFIG="$2"; FORCING="$3"; GOLDEN="$4"; TOL="${5:-1e-10}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE="$SCRIPT_DIR/compare_outputs.sh"

[ -x "$DRIVER" ] || { echo "ERROR: driver not found: $DRIVER" >&2; exit 1; }
[ -f "$CONFIG" ] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
[ -f "$FORCING" ] || { echo "ERROR: forcing not found: $FORCING" >&2; exit 1; }
[ -d "$GOLDEN" ] || { echo "ERROR: golden dir not found: $GOLDEN" >&2; exit 1; }

TMPDIR=$(mktemp -d); trap 'rm -rf "$TMPDIR"' EXIT

ARGS="-c $CONFIG -f $FORCING -q $TMPDIR/q.out -x $TMPDIR/fluxes.out -s $TMPDIR/storage.out -v 0"
[ -f "$GOLDEN/thetas.out" ] && ARGS="$ARGS -t $TMPDIR/thetas.out"

"$DRIVER" $ARGS 2>/dev/null

PASS=0
bash "$COMPARE" "$GOLDEN/q.out"       "$TMPDIR/q.out"       "$TOL" "q.out"       || PASS=1
bash "$COMPARE" "$GOLDEN/fluxes.out"  "$TMPDIR/fluxes.out"  "$TOL" "fluxes.out"  || PASS=1
bash "$COMPARE" "$GOLDEN/storage.out" "$TMPDIR/storage.out" "$TOL" "storage.out" || PASS=1
[ -f "$GOLDEN/thetas.out" ] && { bash "$COMPARE" "$GOLDEN/thetas.out" "$TMPDIR/thetas.out" "$TOL" "thetas.out" || PASS=1; }

exit $PASS
