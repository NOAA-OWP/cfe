# CFE ngen Example Configurations

Example ngen realizations and CFE configs showing the two main operating
modes. Copy and adapt to your catchment.

## Quick-Start: Which mode do I need?

| Mode | CFE computes PET? | AORC forcing needed? | Modules in realization |
|------|-------------------|---------------------|----------------------|
| **External PET** | No | No | SLOTH (or PET module) + CFE |
| **Internal PT** | Yes | Yes | SLOTH + CFE |

## External PET (simplest)

An upstream PET module (or SLOTH placeholder) provides `et_potential_m`.
CFE only needs rainfall, ET, and ice fraction — no AORC fields needed.

**Files:**
- `ngen_realization_cfe_external_pet.json` — ngen realization
- Any existing CFE config works as-is (leave PT and soil evap lines
  commented out or absent)

**What's different from a minimal config:**
- Nothing in the CFE config
- `variables_names_map` in the realization does not need AORC fields
  (CFE dynamically omits them from its BMI input list when PT is disabled)

## Internal Priestley-Taylor PET + Soil Evaporation

CFE computes PET internally using AORC radiation/temperature fields and
a dynamic soil skin temperature model.

**Files:**
- `ngen_realization_cfe_internal_pt.json` — ngen realization
- `cfe_config_internal_pt.cf3` — CFE config (annotated with `<<< NEW`)

**To enable internal PET on an existing config:**

In the CFE config (4 lines):
```
control_simulation_start_date=YYYY-MM-DD          # internal day_of_year
control_ET_simulate_Priestley_Taylor=1.26[]        # enables PT PET
control_soil_simulate_soil_evaporation=TRUE         # enables bare-soil evap
catchment_forested_fraction_0-1=0.70[A A-1]        # required for soil evap
```

In the ngen realization — add AORC mappings to `variables_names_map`:
```json
"TMP_2maboveground": "land_surface_air__temperature",
"DSWRF_surface": "land_surface_radiation~incoming~shortwave__energy_flux",
"DLWRF_surface": "land_surface_radiation~incoming~longwave__energy_flux",
"PRES_surface": "land_surface_air__pressure",
"SPFH_2maboveground": "atmosphere_air_water~vapor__relative_saturation",
"UGRD_10maboveground": "land_surface_wind__x_component_of_velocity",
"VGRD_10maboveground": "land_surface_wind__y_component_of_velocity"
```

**Note on `day_of_year`:** With `control_simulation_start_date` in the
config, CFE computes day_of_year internally and does not advertise it as
a BMI input — no external day_of_year source or mapping needed.

## Land-Cover Parameters

Both modes accept two optional BMI inputs for land-cover partitioning:

| BMI input | Purpose | Default source |
|-----------|---------|---------------|
| `param_catchment_vegetated_fraction` | Fraction of catchment with canopy | `catchment_forested_fraction_0-1` in config |
| `bare_soil_rsurf_exp` | Bare-soil surface resistance exponent | Hardcoded 3.5 |

These can be mapped from SLOTH (as shown in the examples) or from another
upstream module. If not mapped, CFE uses the config/default values.

## Regression Tests

See `test/scripts/` for portable regression test utilities, and
`test/scripts/configs/` for envsubst-based realization templates.
