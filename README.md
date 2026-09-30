# Conceptual Functional Equivalent (CFE) Model

CFE is a simplified conceptual hydrological model designed to be functionally
equivalent to the stormflow generation components of the NOAA National Water
Model (versions 3.1 and earlier). Originally conceived by Fred L. Ogden
(NOAA/NWS Office of Water Prediction).

For the conceptual basis and hypotheses underpinning CFE, see
[MODEL.md](MODEL.md).

## Version 3 Highlights

- **Discrete Soil Moisture Balance Model (DSBM)**: optional 4-layer
  Noah-MP-style soil discretization with Darcy-Buckingham vertical fluxes
  and Clapp-Hornberger hydraulic properties.
- **GIUH surface routing** (Nash Cascade surface routing removed — did not
  add model skill). Subsurface lateral flow via 2-reservoir Nash cascade.
- **Priestley-Taylor PET** estimation from AORC radiation data (testing only).
- **v3 config format** (`.cf3`) with comments and explicit units. Legacy v2
  configs must be converted using `cfe_migrate_config` (see [INSTALL.md](INSTALL.md)).
- **BMI compliance** with the CSDMS BMI-C standard, including `get_value_ptr`
  for all variables and the ngen mass balance protocol.
- **16 calibration parameters** accessible via BMI `set_value`/`get_value_ptr`
  (see [Calibration Parameters](#calibration-parameters) below).

## Build and Run

```bash
cmake -B build -S .
cmake --build build
ctest --test-dir build
```

See [INSTALL.md](INSTALL.md) for detailed build options, driver usage, and
ngen framework integration instructions.

## Configuration

CFE v3 uses a keyword-based config format (`.cf3`). Example configs are in
the `configs/` directory. See [configs/README.md](configs/README.md) for the
full parameter reference, and `configs/clean_config.cf3` for an annotated
template.

## Calibration Parameters

CFE exposes 16 calibration parameters via BMI `set_value` / `get_value_ptr`.
Config files use descriptive key names with units in the name (e.g.,
`soil_sat_hydraulic_conductivity_cm_per_h`). The BMI interface uses shorter
canonical names with values in user-facing units (e.g.,
`soil_saturated_hydraulic_conductivity` in cm/h). Unit conversion to
internal SI happens automatically inside the BMI layer.

| BMI Parameter Name | BMI Unit | Config Key | Config Unit |
|----|------|-----|------|
| `soil_effective_porosity` | - | `soil_effective_porosity` | - |
| `soil_saturated_hydraulic_conductivity` | cm h-1 | `soil_sat_hydraulic_conductivity_cm_per_h` | cm h-1 |
| `soil_percolation_rate_limiter` | - | `soil_to_gw_percolation_rate_limiter_0_to_1` | - |
| `soil_Clapp_Hornberger_b` | - | `soil_Clapp_Hornberger_exponent_b` | - |
| `soil_lateral_flow_K` | h-1 | `soil_reservoir_rate_const_to_subsurface_lateral_flow` | h-1 |
| `subsurface_nash_K` | h-1 | `subsurface_routing_nash_reservoir_time_constant_k` | h-1 |
| `gw_discharge_coefficient` | m s-1 | `gw_discharge_coeff_m_per_timestep` | m timestep-1 |
| `gw_discharge_exponent` | - | `gw_discharge_exponent` | - |
| `gw_max_storage_m` | m | `gw_reservoir_max_storage_m` | m |
| `soil_saturated_capillary_head` | cm | `soil_sat_capillary_head_cm` | cm |
| `soil_field_capacity_fraction` | - | `soil_field_capacity_Pcap_over_Patm_0_1` | - |
| `Xinanjiang_inflection_a` | - | `partitioning_Xinanjiang_tension_water_inflection_point` | - |
| `Xinanjiang_shape_b` | - | `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent` | - |
| `Xinanjiang_shape_x` | - | `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent` | - |
| `Priestley_Taylor_alpha` | - | `control_ET_simulate_Priestley_Taylor` | - |
| `soil_ice_imperv_threshold` | - | `soil_ice_content_impervious_threshold` | - |

Both the BMI interface and config files use the same user-facing units (cm/h
for Ksat, cm for capillary head) so values are human-interpretable. Internal
model computations use SI; conversion happens inside the BMI layer.

Note: `refkdt` is a constant (=3.0, per Schaake et al. 1996) and is not
exposed as a calibration parameter.

## Testing

- **Unit tests**: CTest suite covering all BMI functions, calibration
  parameter round-trips, and mass balance protocol validation.
- **Integration tests**: golden output comparisons at 1e-10 tolerance.

```bash
ctest --test-dir build -L integration   # integration only
ctest --test-dir build -E integration   # unit tests only
```

See [test/README.md](test/README.md) for test code organization.

## Migration from v2

Legacy v2 configs are not supported. Use `cfe_migrate_config` to convert
them (see [INSTALL.md](INSTALL.md#migrating-from-v2-configs)). For the v2
codebase, see the [v2.1.0 tag](https://github.com/NOAA-OWP/cfe/tree/v2.1.0).
See [CHANGELOG.md](CHANGELOG.md) for the complete v2-to-v3 variable name
mapping and list of removed items.

## Getting Help

For questions, please open a GitHub Issue. See [CONTRIBUTING.md](CONTRIBUTING.md)
for development guidelines.
