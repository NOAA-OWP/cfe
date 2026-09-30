#!/bin/bash
# run_restart_test.sh — Verify hotstart restart produces identical output to a continuous run
#
# Splits a forcing file at a given step, runs the standalone driver for:
#   1. A continuous reference run (full forcing)
#   2. First half with hotstart generation
#   3. Second half from hotstart config
# Compares the second-half output values against the corresponding
# portion of the reference run.
#
# Usage: run_restart_test.sh <driver> <config> <forcing> <split_at_step> [tolerance]
set -e

DRIVER="$1"; CONFIG="$2"; FORCING="$3"; SPLIT="$4"; TOL="${5:-1e-10}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE="$SCRIPT_DIR/compare_outputs.sh"

[ -x "$DRIVER" ]  || { echo "ERROR: driver not found: $DRIVER" >&2; exit 1; }
[ -f "$CONFIG" ]   || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
[ -f "$FORCING" ]  || { echo "ERROR: forcing not found: $FORCING" >&2; exit 1; }
[ -n "$SPLIT" ]    || { echo "ERROR: split_at_step required" >&2; exit 1; }

WORKDIR=$(mktemp -d); trap 'rm -rf "$WORKDIR"' EXIT

# ---- Split forcing file (header + data) ----
head -1 "$FORCING" > "$WORKDIR/header.csv"
cat "$WORKDIR/header.csv" > "$WORKDIR/forcing_1.csv"
tail -n +2 "$FORCING" | head -n "$SPLIT" >> "$WORKDIR/forcing_1.csv"
cat "$WORKDIR/header.csv" > "$WORKDIR/forcing_2.csv"
tail -n +$((SPLIT + 2)) "$FORCING" >> "$WORKDIR/forcing_2.csv"

# ---- Prepare configs ----
# Strip output/forcing/hotstart settings from base config, then add our own.
# The hotstart config inherits these, so we strip and re-add for it too.
strip_output_keys() {
    grep -v -E '^[[:space:]]*(output_path_name|output_internal_fluxes|output_internal_storages|output_volume_balance|output_soil_moisture_theta|output_discharge|output_total_discharge|output_new_config|output_time_standard|output_file_delimiter|output_value_format|control_input_forcing)' "$1"
}
add_output_keys() {
    cat >> "$1" << 'EOF'
control_input_forcing_filename=BMI
output_time_standard_format=timestep
output_file_delimiter=space
output_value_format="%.8e"
EOF
}

strip_output_keys "$CONFIG" > "$WORKDIR/base.cf3"
add_output_keys "$WORKDIR/base.cf3"

cp "$WORKDIR/base.cf3" "$WORKDIR/config_ref.cf3"
cp "$WORKDIR/base.cf3" "$WORKDIR/config_1.cf3"
echo "output_new_config_filename_prefix=\"$WORKDIR/hotstart\"" >> "$WORKDIR/config_1.cf3"

# ---- Run ----
OUTPUTS="q fluxes storage thetas"
run_driver() {
    local cfg="$1" forcing="$2" prefix="$3"
    "$DRIVER" -c "$cfg" -f "$forcing" \
        -q "$WORKDIR/${prefix}_q.out" \
        -x "$WORKDIR/${prefix}_fluxes.out" \
        -s "$WORKDIR/${prefix}_storage.out" \
        -t "$WORKDIR/${prefix}_thetas.out" \
        -v 0 2>/dev/null
}

# Reference (continuous)
run_driver "$WORKDIR/config_ref.cf3" "$FORCING" ref

# First half (with hotstart)
run_driver "$WORKDIR/config_1.cf3" "$WORKDIR/forcing_1.csv" half1

# Find and prepare hotstart config
HOTSTART=$(ls "$WORKDIR"/hotstart.*.cf3 2>/dev/null | head -1)
[ -f "$HOTSTART" ] || { echo "ERROR: hotstart config not generated" >&2; exit 1; }

strip_output_keys "$HOTSTART" > "$WORKDIR/config_2.cf3"
add_output_keys "$WORKDIR/config_2.cf3"

# Second half (from hotstart)
run_driver "$WORKDIR/config_2.cf3" "$WORKDIR/forcing_2.csv" half2

# ---- Compare second-half values against reference ----
# Output format: "step value1 value2 ..." with 2 comment lines at top.
# Reference second half starts at data line SPLIT+1 = file line SPLIT+3.
# Second-half output has its own 2 comment lines then data starting at 0.
PASS=0
for name in $OUTPUTS; do
    ref="$WORKDIR/ref_${name}.out"
    h2="$WORKDIR/half2_${name}.out"
    [ -f "$ref" ] && [ -f "$h2" ] || continue

    # Extract value columns only (skip timestep = field 1)
    tail -n +$((SPLIT + 3)) "$ref" | grep -v '^[[:space:]]*#' | cut -d' ' -f2- > "$WORKDIR/ref_vals"
    grep -v '^[[:space:]]*#' "$h2" | cut -d' ' -f2- > "$WORKDIR/h2_vals"

    bash "$COMPARE" "$WORKDIR/ref_vals" "$WORKDIR/h2_vals" "$TOL" "restart_${name}" || PASS=1
done

exit $PASS
