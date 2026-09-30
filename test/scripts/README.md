# test/scripts — Regression Testing Utilities

Portable scripts for CFE ngen integration regression testing. These are
independent of any specific test data location — pass paths as arguments.

## Scripts

### compare_discharge.py

Compare two discharge time series (any combination of ngen CSV, standalone
driver output, or Fred-format CSV). Auto-detects format. Reports exact
matches, truncation-noise matches, max abs/rel diffs, pass/fail.

```bash
# Compare ngen output against standalone driver baseline
python3 compare_discharge.py reference/q.out output/cat-17536.csv --tol 1e-6

# Explicit format override
python3 compare_discharge.py ref.csv test.csv --ref-fmt fred --test-fmt ngen
```

### run_ngen_regression.sh

Run a single ngen regression test. Expects the standard test directory layout:

```
test_dir/
  realization.json
  catchment_data.geojson
  nexus_data.geojson
  reference/          # contains q.csv, q.out, or q_baseline.csv
  output/             # cleared and repopulated each run
```

```bash
# Run test and compare against existing reference
bash run_ngen_regression.sh /path/to/test_dir --ngen /path/to/ngen --tol 1e-6

# Generate a new baseline from current build
bash run_ngen_regression.sh /path/to/test_dir --ngen /path/to/ngen --baseline
```

### generate_standalone_baseline.sh

Generate full-precision reference output from `cfe_bmi_driver` (the standalone
driver). Produces q.out, fluxes.out, storage.out, thetas.out.

```bash
bash generate_standalone_baseline.sh \
    --driver build/cfe_bmi_driver \
    --config test_dir/cfe_config/standalone.cf3 \
    --forcing test_dir/forcing/cat-17536.csv \
    --outdir test_dir/reference
```

### bisect_regression.sh

Git bisect helper. Builds CFE, runs a regression test, returns appropriate
exit codes for `git bisect run`.

```bash
git bisect start <bad-commit> <good-commit>
git bisect run bash test/scripts/bisect_regression.sh \
    --test-dir /path/to/test_dir \
    --ngen /path/to/ngen \
    --tol 1e-6
```

Exit codes: 0 = good, 1 = bad, 125 = skip (build failed).

## Config Templates (`configs/`)

Realization templates and CFE config files live in `configs/`. Environment-
specific paths use `$VAR` syntax — edit `env.sh` once, then render with
`envsubst`:

```bash
cp configs/env.sh my_env.sh && vi my_env.sh
source my_env.sh
envsubst < configs/realization_pt_enabled.json > my_realization.json
```

See `configs/README.md` for the full variable reference.

## Workflow: Establishing a Regression Baseline

1. Build CFE from the known-good commit
2. Render realization templates (see above)
3. Generate baselines:
   - For ngen tests: `run_ngen_regression.sh <test_dir> --baseline`
   - For standalone tests: `generate_standalone_baseline.sh --driver ... --outdir ...`
4. Subsequent builds compare against these baselines automatically

## Workflow: Comparing Against an External Reference

The comparator auto-detects three output formats, so you can drop in a reference
file from any source — a standalone driver binary, an earlier build, etc.

**Supported reference formats:**
- **ngen CSV**: header row + `step,time,value` columns
- **Standalone driver**: `#` comment lines + `timestep value` (space-delimited)
- **Fred-format CSV**: `#` comment lines + `datetime,value` (comma-delimited)

**Steps:**

```bash
# 1. Copy the external reference into the test's reference/ directory
cp /path/to/freds_q.csv my_test/reference/freds_q.csv

# 2. Run ngen and compare, pointing --ref at the external file
bash test/scripts/run_ngen_regression.sh my_test/ \
    --ngen /path/to/ngen --ref freds_q.csv --tol 1e-6
```

Use `--ref <filename>` to name the reference file explicitly (looked up inside
`reference/`). Without `--ref`, the runner auto-detects in priority order:
`q_baseline.csv` → `q.csv` → `q.out`. This lets you keep both a self-regression
baseline and an external reference side by side:

```bash
my_test/reference/
  q_baseline.csv      # self-regression (exact match expected)
  external_q.csv      # external (use --ref external_q.csv, expect ~1e-6 diffs)
```

You can also compare two files directly without running ngen:

```bash
python3 test/scripts/compare_discharge.py \
    /path/to/freds_q.csv /path/to/ngen_output.csv --tol 1e-6
```

**Tolerance guidance:**
- Self-regression (same build): `1e-10` — expect exact match
- Cross-binary (standalone driver vs ngen BMI): `1e-6` — FP accumulation diffs
- Cross-version (earlier build vs current): `1e-6` to `1e-4` — depends on code changes

## Existing Scripts (test/ root)

The `test/` directory also contains:

- **compare_outputs.sh** — awk-based numerical comparator for standalone driver
  output files (used by CTest integration tests)
- **run_integration_test.sh** — CTest-registered standalone driver golden tests
- **run_restart_test.sh** — Hotstart/restart validation for the standalone driver

Those scripts are wired into `CMakeLists.txt` and run via `ctest`. The scripts
in this directory are for ngen-level regression testing, which is outside the
CTest suite.
