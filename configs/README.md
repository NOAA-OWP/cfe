# CFE v3 Configuration Files

CFE v3 uses a keyword-based config format (`.cf3`). Legacy v2 configs are not
supported directly — use `cfe_migrate_config` to convert them. For the v2
codebase, see the [v2.1.0 tag](https://github.com/NOAA-OWP/cfe/tree/v2.1.0).

## Config Format

Comments with `#` or `//`, one keyword per line, optional unit annotations in
brackets. The bracketed units are for human readability only — the parser
strips them but does not interpret or convert based on them. Values must
always be in the units indicated by the key name:
```
cfe_config_version=3.0[]
control_model_timestep_h=1.0[h]
soil_depth_m=2.0[m]
soil_sat_hydraulic_conductivity_cm_per_h=1.2168[cm h-1]
control_soil_simulate_discrete_soil_moisture_true_false=TRUE
```

See `clean_config.cf3` for a fully annotated template of all options.

## Example Configs

| File | Partitioning | Soil | Notes |
|------|-------------|------|-------|
| `cfe_config_cat_87_pass.cf3` | Schaake | Linear reservoir | Used by unit tests |
| `migrated_legacy_cat87.cf3` | Xinanjiang | Linear reservoir | Integration test (migrated from v2) |
| `bmi_config_cat87_v3.cf3` | Schaake | DSBM (4-layer) | Integration test, DSBM enabled |
| `cfe_config_cat87_v3.cf3` | Schaake | DSBM (4-layer) | Standalone driver config |
| `laramie_dry_start_v3.cf3` | Schaake | DSBM (4-layer) | Integration test, dry start exercises wilting point ET reduction |
| `clean_config.cf3` | — | — | Annotated template |

## Migrating from v2

```bash
build/cfe_migrate_config old_config.txt new_config.cf3
```

The migration utility converts key names, units (satdk m/s to cm/h, satpsi m
to cm), and storage representations (theta to meters). Configs that relied on
Nash Cascade surface routing will error — GIUH ordinates are required in v3.

## Parameter Reference

Parameters marked with \* are calibratable via BMI `set_value` /
`get_value_ptr`. The legend for table columns:

- **Config Key**: the keyword used in `.cf3` config files
- **Config Unit**: the unit expected in the config file (embedded in the key name)
- **BMI Name**: the variable name used with BMI `set_value` / `get_value` / `get_value_ptr`
- **BMI Unit**: the unit reported by BMI `get_var_units` (internal model units)
- **Range**: calibration range from NWM experience (OWP/FLO)

### Config Units vs BMI Units

Config files use human-readable units so hydrologists can edit them directly
(e.g., `soil_sat_hydraulic_conductivity_cm_per_h=1.2168`). The expected unit
for each parameter is embedded in the key name — the `[unit]` brackets in
config files are annotations for human readers and are not parsed.

The BMI interface uses the same user-facing units as the config file for
all calibration parameters. The only exception is the groundwater
discharge coefficient, which uses per-timestep units in the config but
per-second units in BMI:

| Parameter | Config Unit | BMI Unit |
|-----------|-------------|----------|
| GW discharge coefficient | m timestep-1 | m s-1 |

### Soil Parameters

| Config Key | Config Unit | BMI Name | BMI Unit | Range | Description |
|-----------|-------------|----------|----------|-------|-------------|
| `soil_depth_m` | m | `param_soil_depth_m` | m | — | Soil column depth |
| `soil_Clapp_Hornberger_exponent_b`\* | - | `soil_Clapp_Hornberger_b` | - | 2.0–15.0 | Clapp-Hornberger exponent |
| `soil_sat_hydraulic_conductivity_cm_per_h`\* | cm h-1 | `soil_saturated_hydraulic_conductivity` | cm h-1 | 0.07–510.0 cm/h | Saturated hydraulic conductivity |
| `soil_sat_capillary_head_cm`\* | cm | `soil_saturated_capillary_head` | cm | — | Saturated capillary head |
| `soil_effective_porosity`\* | - | `soil_effective_porosity` | - | 0.16–0.58 | Effective porosity |
| `soil_to_gw_percolation_rate_limiter_0_to_1`\* | - | `soil_percolation_rate_limiter` | - | 0.0–1.0 | Percolation rate limiter |
| `soil_field_capacity_Pcap_over_Patm_0_1`\* | - | `soil_field_capacity_fraction` | - | 0.15–0.33 | Field capacity; typically not calibrated |
| `soil_reservoir_rate_const_to_subsurface_lateral_flow`\* | h-1 | `soil_lateral_flow_K` | h-1 | 0.0–1.0 | Lateral flow rate constant |
| `state_soil_reservoir_init_storage_m` | m | — | — | — | Initial soil storage |

Note: `refkdt` is a constant (=3.0, per Schaake et al. 1996) hardcoded in the
model. It is not a calibration parameter — calibrating both `refkdt` and `satdk`
is degenerate (a change in one can be fully compensated by the other).

### Groundwater Parameters

| Config Key | Config Unit | BMI Name | BMI Unit | Range | Description |
|-----------|-------------|----------|----------|-------|-------------|
| `gw_reservoir_max_storage_m`\* | m | `gw_max_storage_m` | m | 0.1–0.25 | Maximum GW storage |
| `gw_discharge_coeff_m_per_timestep`\* | m timestep-1 | `gw_discharge_coefficient` | m s-1 | 1.8e-6–1.8e-3 | GW discharge coefficient |
| `gw_discharge_exponent`\* | - | `gw_discharge_exponent` | - | 1.0–8.0 | GW discharge exponent |
| `state_gw_reservoir_init_storage_m` | m | — | — | — | Initial GW storage |

### GIUH Surface Routing

| Config Key | Units | Description |
|-----------|-------|-------------|
| `surface_routing_num_giuh_ordinates` | — | Number of GIUH ordinates (required) |
| `surface_routing_giuh_ordinates` | - | Comma-separated ordinates summing to 1.0 |
| `state_surface_routing_init_giuh_convolution_queue_m` | m | Initial convolution queue (one per ordinate) |

GIUH is the only surface routing method in CFE v3. Nash Cascade surface routing
was removed as it did not produce added model skill.

### Subsurface Routing (Nash Cascade)

| Config Key | Config Unit | BMI Name | BMI Unit | Range | Description |
|-----------|-------------|----------|----------|-------|-------------|
| `subsurface_routing_nash_reservoir_time_constant_k`\* | h-1 | `subsurface_nash_K` | h-1 | 0.0–1.0 | Subsurface Nash time constant |
| `state_subsurface_routing_init_nash_cascade_storage_m` | m | — | — | — | Initial Nash cascade storage (2 values) |

### Xinanjiang Parameters (when `partitioning_scheme_name=XINANJIANG`)

| Config Key | Config Unit | BMI Name | BMI Unit | Range |
|-----------|-------------|----------|----------|-------|
| `partitioning_Xinanjiang_tension_water_inflection_point`\* | - | `Xinanjiang_inflection_a` | - | -0.493–0.493 |
| `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent`\* | - | `Xinanjiang_shape_b` | - | 0.0–1.0 |
| `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent`\* | - | `Xinanjiang_shape_x` | - | 0.0–1.0 |

### Additional Calibration Parameters

| Config Key | Config Unit | BMI Name | BMI Unit | Range | Description |
|-----------|-------------|----------|----------|-------|-------------|
| `control_ET_simulate_Priestley_Taylor`\* | - | `Priestley_Taylor_alpha` | - | 0.75–1.6 | P-T alpha coefficient (0 disables PT PET) |
| `soil_ice_content_impervious_threshold`\* | - | `soil_ice_imperv_threshold` | - | — | Ice fraction threshold for Schaake |

### DSBM and v3 Feature Controls

| Config Key | Description |
|-----------|-------------|
| `control_soil_simulate_discrete_soil_moisture_true_false` | `TRUE` to enable DSBM |
| `control_soil_use_lookup_table_num_points` | LUT points (0 = analytic) |
| `control_ET_deepest_root_zone_discretization` | Deepest root zone layer (1–NDISC) |

### Catchment Metadata (debug only)

| Config Key | Unit | Description |
|-----------|------|-------------|
| `catchment_latitude_decimal_degree` | degrees | Latitude of catchment centroid |
| `catchment_longitude_decimal_degree` | degrees | Longitude of catchment centroid |
| `catchment_elevation` | m | Mean catchment elevation |

These fields are parsed and included in diagnostic output (`print_config`) but
are **not used by the model**. They are retained for provenance and debugging.

Calibration ranges are from NWM calibration experience (OWP/FLO). Parameters
marked `—` have no established calibration range or are not typically calibrated.

## Infiltration Excess Partitioning

1. **Schaake** (`partitioning_scheme_name=SCHAAKE`) — no additional parameters
2. **Xinanjiang** (`partitioning_scheme_name=XINANJIANG`) — requires the three Xinanjiang parameters above
