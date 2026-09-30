# CFE v3 — Build, Test, and Run Instructions

CFE v3 supersedes all previous versions (1.x, 2.x). For legacy v2 needs, see
the [v2.1.0 tag](https://github.com/NOAA-OWP/cfe/tree/v2.1.0).

Building CFE requires a C compiler (GCC or Clang) and [CMake](https://cmake.org/) >= 3.10.

## Quick Start

```bash
git clone https://github.com/NOAA-OWP/cfe
cd cfe
cmake -B build -S .
cmake --build build
ctest --test-dir build
```

This builds the shared library (`libcfebmi`), the BMI driver executable
(`cfe_bmi_driver`), the config migration utility (`cfe_migrate_config`),
and all test executables, then runs the full test suite (unit tests +
integration tests with golden output comparison).

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `STANDALONE` | `OFF` | Build the standalone non-BMI driver (`cfe_main_driver`) |
| `NGEN` | `ON` | Accepted for compatibility (no effect — always builds) |
| `CMAKE_BUILD_TYPE` | (none) | Set to `Debug` for debug symbols and verbosity=1 |

Example with standalone driver:
```bash
cmake -B build -S . -DSTANDALONE=ON
cmake --build build
```

## Running Tests

### All tests
```bash
ctest --test-dir build
```

### Unit tests only
```bash
ctest --test-dir build -E integration
```

### Integration tests only
```bash
ctest --test-dir build -L integration
```

The integration tests run the `cfe_bmi_driver` with v3 configs (bucket and
DSBM) against golden reference outputs at 1e-10 tolerance (exact match).

### Verbose output on failure
```bash
ctest --test-dir build --output-on-failure
```

## Running the BMI Driver

The `cfe_bmi_driver` executable reads a config file and forcing data, runs
the model through the BMI interface, and writes output files.

> **Note:** Output directories are not created automatically. Create them
> before running the driver:
> ```bash
> mkdir -p output
> ```

### Example
```bash
build/cfe_bmi_driver \
    -c configs/bmi_config_cat87_v3.cf3 \
    -f forcings/cat87_01Dec2015.csv \
    -q output/q.out \
    -x output/fluxes.out \
    -s output/storage.out \
    -t output/thetas.out \
    -v 1
```

### Driver options
```
Usage:
  cfe_bmi_driver -c <config> -f <forcing> [OPTIONS]

Required:
  -c <file>     Configuration file (.cf3)

Optional:
  -f <file>     Forcing data file (overrides config)
  -q <file>     Discharge output file (m/timestep)
  -b <file>     Volume balance summary file
  -x <file>     Internal fluxes output file
  -s <file>     Internal storages output file
  -t <file>     Soil moisture theta output (DSBM only)
  -v <level>    Verbosity (0=quiet, 1=normal, 2=verbose)
  -dryrun       Run 120 timesteps without forcing data
```

## Config File Format

CFE v3 uses a keyword-based config format with `cfe_config_version=3.0`.
Comments with `#` or `//`, one keyword per line, optional unit annotations in
brackets (for human readability only — not parsed for conversion):
```
cfe_config_version=3.0[]
control_model_timestep_h=1.0[h]
control_input_forcing_filename=BMI
soil_depth_m=2.0[m]
soil_Clapp_Hornberger_exponent_b=4.05[]
soil_sat_hydraulic_conductivity_cm_per_h=1.2168[cm h-1]
control_soil_simulate_discrete_soil_moisture_true_false=TRUE
...
```

See `configs/clean_config.cf3` for a fully annotated template of all options.

### Migrating from v2 configs

CFE v3 does not support legacy v2 config files. Use the migration utility:
```bash
build/cfe_migrate_config old_config.txt new_config.cf3
```

This converts v2 key names, units, and storage representations to v3 format.
Configs that relied on Nash Cascade surface routing (removed in v3) will
produce an error — GIUH ordinates are required.

## ngen Framework Integration

CFE v3 builds a shared library (`libcfebmi.dylib` / `libcfebmi.so`) that
can be loaded by the [ngen](https://github.com/NOAA-OWP/ngen) framework.

### Building within ngen
```bash
# From ngen directory
git submodule update --remote extern/cfe/cfe
cmake -B extern/cfe/cfe/cmake_build -S extern/cfe/cfe/
make -C extern/cfe/cfe/cmake_build
```

The realization config should reference the shared library and a v3 config:

```json
{
    "model_type_name": "bmi_c",
    "library_file": "./extern/cfe/cfe/cmake_build/libcfebmi.so",
    "init_config": "./extern/cfe/cfe/configs/bmi_config_cat87_v3.cf3"
}
```

### BMI Variable Names

**Inputs:** `rainfall_depth_m` (m s-1), `et_potential_m` (m s-1)

**Key outputs:** `discharge_m`, `surface_runoff_m`, `lateral_flow_m`,
`baseflow_m`, `actual_et_m`, `potential_et_m`, `giuh_outflow_m`,
`soil_to_gw_percolation_flux_m` (all m/timestep)

**Calibration parameters** (16 total, accessible via `set_value` / `get_value_ptr`):
BMI parameter names use user-facing units (cm/h for conductivity, cm for
head), matching the config file conventions.
See [README.md](README.md#calibration-parameters) for the full config-to-BMI
mapping table, or `bmi_cfe.c` `param_var_names[]` for the authoritative list.

Note: `refkdt` is a constant (=3.0, per Schaake et al. 1996) and is not
calibratable.

### Notes

- `include/bmi.h` must be the canonical [CSDMS BMI-C](https://github.com/csdms/bmi-c)
  release header. ABI compatibility with ngen cannot be guaranteed if this file
  is modified or replaced with a non-standard version.
- The `-DNGEN=ON` flag is accepted for compatibility but is no longer
  required. A plain `cmake -B build -S .` builds everything.
