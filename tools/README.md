# cfe_migrate_config — v2 to v3 Config Migration

Converts a legacy CFE v2 config file to v3 format. For the v2 codebase, see
the [v2.1.0 tag](https://github.com/NOAA-OWP/cfe/tree/v2.1.0).

## Usage

```bash
cfe_migrate_config <input_v2_config> <output_v3_config>
```

The utility is built automatically with the rest of CFE:
```bash
cmake -B build -S .
cmake --build build
./build/cfe_migrate_config old_config.txt new_config.cf3
```

## What it does

1. Reads v2 `key=value[units]` format
2. Renames keys to v3 names
3. Converts units where needed
4. Writes a v3 config with `cfe_config_version=3.0`
5. Adds v3-only control keys with sensible defaults

If the v2 config uses Nash Cascade surface routing without GIUH ordinates,
the utility errors — GIUH is required in v3.

If the input already has a `cfe_config_version` key, it is detected as v3 and
skipped (exit 0, no output written).

## Key Mapping

| v2 Key | v3 Key | Notes |
|--------|--------|-------|
| `forcing_file` | `control_input_forcing_filename` | |
| `num_timesteps` | `control_total_num_simulation_timesteps` | |
| `verbosity` | `control_verbosity` | |
| `soil_params.depth` | `soil_depth_m` | |
| `soil_params.b` | `soil_Clapp_Hornberger_exponent_b` | |
| `soil_params.satdk` | `soil_sat_hydraulic_conductivity_cm_per_h` | **m/s → cm/h** (×360,000) |
| `soil_params.satpsi` | `soil_sat_capillary_head_cm` | **m → cm** (×100) |
| `soil_params.smcmax` | `soil_effective_porosity` | |
| `soil_params.wltsmc` | *(dropped)* | v3 auto-calculates from Clapp-Hornberger; v2 value is not migrated |
| `soil_params.slop` | `soil_to_gw_percolation_rate_limiter_0_to_1` | |
| `alpha_fc` | `soil_field_capacity_Pcap_over_Patm_0_1` | |
| `K_lf` | `soil_reservoir_rate_const_to_subsurface_lateral_flow` | |
| `soil_storage` | `state_soil_reservoir_init_storage_m` | **theta → meters** (×depth×porosity) |
| `max_gw_storage` | `gw_reservoir_max_storage_m` | |
| `Cgw` | `gw_discharge_coeff_m_per_timestep` | |
| `expon` | `gw_discharge_exponent` | |
| `gw_storage` | `state_gw_reservoir_init_storage_m` | |
| `K_nash_subsurface` | `subsurface_routing_nash_reservoir_time_constant_k` | |
| `nash_storage_subsurface` | `state_subsurface_routing_init_nash_cascade_storage_m` | |
| `giuh_ordinates` | `surface_routing_giuh_ordinates` | Count auto-detected |
| `surface_water_partitioning_scheme` | `partitioning_scheme_name` | |
| `a_Xinanjiang_inflection_point_parameter` | `partitioning_Xinanjiang_tension_water_inflection_point` | |
| `b_Xinanjiang_shape_parameter` | `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent` | |
| `x_Xinanjiang_shape_parameter` | `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent` | |
| `urban_decimal_fraction` | `catchment_impervious_fraction_0-1` | |

## Unit Conversions

| Parameter | v2 Units | v3 Units | Conversion |
|-----------|----------|----------|------------|
| `satdk` | m/s | cm/h | × 100 × 3600 |
| `satpsi` | m | cm | × 100 |
| `soil_storage` | dimensionless theta | m | × soil_depth × porosity |

## v3 Defaults Added

The following keys are added to the output with default values when not
present in the v2 config:

| Key | Default | Description |
|-----|---------|-------------|
| `cfe_config_version` | `3.0` | Always added |
| `control_model_timestep_h` | `1.0` | v2 default timestep |
| `control_soil_simulate_discrete_soil_moisture_true_false` | `FALSE` | DSBM off |
| `control_soil_simulate_freeze_thaw_true_false` | `FALSE` | SFT off |
| `control_ET_deepest_root_zone_discretization` | `4` | |
| `control_soil_use_lookup_table_num_points` | `0` | Analytic mode |
| `soil_ice_content_impervious_threshold` | `0.0` | |
| `state_surface_routing_init_giuh_convolution_queue_m` | all zeros | Matches ordinate count |

## Deprecated v2 Keys (skipped with a note)

These v2 keys are recognized and silently skipped — they have no v3 equivalent:

- `soil_params.expon`, `soil_params.expon_secondary`
- `refkdt`, `debug`
- `surface_runoff_scheme`, `N_nash_surface`, `K_nash_surface`
- `nash_storage_surface`, `nsubsteps_nash_surface`
- `Kinf_nash_surface`, `retention_depth_nash_surface`

## Wilting Point Auto-Calculation

If `soil_params.wltsmc` is not in the v2 config, the utility calculates it
from the Clapp-Hornberger relation at 15 atmospheres:

```
theta_wp = porosity × (satpsi_cm / psi_wp_cm) ^ (1/b)
```

where `psi_wp_cm = 15 atm` converted to cm of water head.
