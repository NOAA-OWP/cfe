# Regression Test Config Templates

CFE config files and ngen realization templates for regression testing.
Environment-specific paths use `$VAR` syntax for `envsubst` rendering.

## Quick Start

```bash
# 1. Copy the env file and edit paths for your environment
cp env.sh my_env.sh
vi my_env.sh

# 2. Source it and render a realization
source my_env.sh
envsubst < realization_pt_disabled.json > my_realization.json

# 3. Copy the CFE config (no substitution needed)
cp cfe_pt_disabled.cf3 /path/to/test/cfe_config/cat-17536_config.cf3
```

## Files

| File | Tests | Description |
|------|-------|-------------|
| `env.sh` | all | Environment variables — copy and edit for your setup |
| `cfe_pt_enabled.cf3` | 2 | PT + soil evap enabled; internal day_of_year via `control_simulation_start_date` |
| `cfe_pt_disabled.cf3` | 3, 4, 6 | PT + soil evap disabled; NOM or external module provides PET |
| `realization_pt_enabled.json` | 2 | SLOTH + CFE with AORC mappings; internal day_of_year |
| `realization_pt_disabled.json` | 6 | SLOTH + CFE only; simplest ngen configuration |
| `realization_nom_coupled.json` | 3, 4 | SLOTH + NOM + CFE; add `model_params` for test 4 |

## Environment Variables (env.sh)

| Variable | Example | Used In |
|----------|---------|---------|
| `SLOTH_LIB` | `./extern/sloth/cmake_build/libslothmodel` | all realizations |
| `NOM_LIB` | `./extern/noah-owp-modular/cmake_build/libsurfacebmi` | nom_coupled |
| `NOM_CONFIG` | `./data/test/nom_config/{{id}}.namelist.input` | nom_coupled |
| `CFE_LIB` | `/abs/path/to/cfe/build/libcfebmi` | all realizations |
| `CFE_CONFIG` | `./data/test/cfe_config/{{id}}_config.cf3` | all realizations |
| `FORCING_DIR` | `./data/test/forcing/` | all realizations |
| `OUTPUT_DIR` | `./data/test/output/` | all realizations |
| `START_TIME` | `2012-10-01 00:00:00` | all realizations |
| `END_TIME` | `2021-08-31 00:00:00` | all realizations |

`{{id}}` in `CFE_CONFIG` and `NOM_CONFIG` is ngen's own catchment-ID
substitution — `envsubst` only expands `$VAR`/`${VAR}` syntax, so the
double-brace `{{id}}` passes through untouched.

## Standalone Driver Use

CFE configs default to `control_input_forcing_filename=BMI` for ngen.
For standalone driver use:

```bash
sed 's|control_input_forcing_filename=BMI|control_input_forcing_filename=/path/to/forcing.csv|' \
    cfe_pt_disabled.cf3 > standalone_config.cf3
```

## SLOTH Constants

| Parameter | Value | Notes |
|-----------|-------|-------|
| `sloth_ice_fraction` | 0.0 | No frozen soil |
| `sloth_smp` | 0.0 | No soil moisture profile coupling |
| `sloth_et_potential_m` | 0.0 | Placeholder; PT or NOM provides actual PET |
| `sloth_vegetated_fraction` | 0.70 | Must match `catchment_forested_fraction_0-1` in CFE config |
| `sloth_rsurf_exp` | 5.0 | Bare soil surface resistance exponent |

In NOM-coupled realizations, SLOTH only provides `ice_fraction` and `smp` — NOM
provides the rest through its own outputs (`FVEG`, `RSURF_EXP`, `EVAPOTRANS`).
