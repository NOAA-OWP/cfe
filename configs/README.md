## Configuration File
Example configuration files are provided in this directory. To build and run the given examples see the instructions [here](https://github.com/NOAA-OWP/cfe/blob/master/INSTALL.md). A detailed description of the parameters for model configuration (i.e., initialize/setup) is provided below. The asterisk (*) denotes calibratable parameters.

| Variable | Datatype |  Limits  | Units | Role | Process | Description |
| -------- | -------- | ------ | ----- | ---- | ------- | ----------- |
| forcing_file | *char* | 256  |   | filename |   | path to forcing inputs csv; set to `BMI` if passed via `bmi.set_value*()`  |
| soil_params.depth | *double* |   | meters [m]| state |  | soil depth  |
| *soil_params.b | *double* |   |   | state |   | beta exponent on Clapp-Hornberger (1978) soil water relations  |
| *soil_params.satdk | *double* |   |  meters/second [m s-1] | state |  | saturated hydraulic conductivity  |
| *soil_params.satpsi  | *double* |   |  meters [m] | state |  | saturated capillary head  |
| *soil_params.slop   | *double* |   |  meters/meters [m/m]| state |  | this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow. |
| *soil_params.smcmax  | *double* |   |  meters/meters [m/m] | state |  | saturated soil moisture content  |
| soil_params.wltsmc | *double* |   |  meters/meters [m/m] | state |   | wilting point soil moisture content  |
| soil_params.expon  | *double* |   |  | parameter_adjustable |    | optional; defaults to `1.0`  |
| soil_params.expon_secondary  | *double* |  |   | parameter_adjustable |  | optional; defaults to `1.0` |
| *max_gw_storage | *double* |   |  meters [m] | parameter_adjustable |  | maximum storage in the conceptual reservoir |
| *Cgw | *double* |   |  meters/hour [m h-1] | parameter_adjustable |  | the primary outlet coefficient |
| *expon | *double* |   |   | parameter_adjustable |  | exponent parameter (1.0 for linear reservoir) |
| gw_storage | *double* |   |  meters/meters [m/m] | parameter_adjustable |  | initial condition for groundwater reservoir - it is the ground water as a decimal fraction of the maximum groundwater storage (max_gw_storage) for the initial timestep |
| alpha_fc | *double* |   |   | parameter_adjustable |  | field capacity |
| soil_storage| *double* |   | meters/meters [m/m] | parameter_adjustable |  | initial condition for soil reservoir - it is the water in the soil as a decimal fraction of maximum soil water storage (smcmax * depth) for the initial timestep |
| N_nash | *int* |   |   | parameter_adjustable |   | number of Nash lf reservoirs (optional, defaults to 2, ignored if storage values present)  |
| K_nash | *double* |   | 1/meters [m^-1]  | parameter_adjustable | subsurface runoff | Nash Config param for lateral subsurface runoff   |
| *K_lf | *double* |   |   | parameter_adjustable |  | Nash Config param - primary reservoir  |
| nash_storage | *double* |   |   | parameter_adjustable |  | Nash Config param - secondary reservoir   |
| giuh_ordinates   | *double* |   |   | parameter_adjustable |  | Giuh ordinates in dt time steps   |
| num_timesteps  | *int* |   |  | time_info |  | set to `1` if `forcing_file=BMI`   |
| verbosity | *int* | `0`-`3`  |   | option |   |  prints various debug and bmi info  |
| surface_partitioning_scheme | *char* | `Xinanjiang` or `Schaake`  |  | parameter_adjustable | direct runoff |    |
| surface_runoff_scheme | *char* | GIUH or NASH_CASCADE | | parameter_adjustable | surface runoff | also supports 1 for GIUH and 2 for NASH_CASCADE; default is GIUH |
| N_nash_surface | *int* |   |   | parameter_adjustable | surface runoff | number of Nash reservoirs for surface runoff   |
| K_nash_surface | *double* |   | 1/meters [m^-1]  | parameter_adjustable | surface runoff | Nash Config param for surface runoff   |
| nash_storage_surface | *double* |   | meters [m]  | parameter_adjustable | surface runoff | Nash Config param; reservoir surface storage; default is zero storage |
| *a_Xinanjiang_inflection_point_parameter | *double* |   |  | parameter_adjustable | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| *b_Xinanjiang_shape_parameter=1  | *double* |   |   | parameter_adjustable  | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| *x_Xinanjiang_shape_parameter=1  | *double* |   |   | parameter_adjustable | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| urban_decimal_fraction  | *double*  |  0.0 - 1.0 |   |  parameter_adjustable | direct runoff | when `surface_partitioning_scheme=Xinanjiang` |
| is_aet_rootzone                    | *boolean* | True, true or 1  |  | coupling parameter | `rootzone-based AET` | when `CFE coupled to SoilMoistureProfile` |
| max_rootzone_layer | *double* |  | meters [m] | parameter_adjustable | AET | layer of the soil that is the maximum root zone depth. That is, the depth of the layer where the AET is drawn from |
| soil_layer_depths | 1D array |  | meters [m] | parameter_adjustable | AET | an array of depths from the surface. Example, soil_layer_depths=0.1,0.4,1.0,2.0
| is_sft_coupled                   | *boolean* | True, true or 1  |  | coupling parameter | `ice_fraction-based runoff` | when `CFE coupled to SoilFreezeThaw`|

## Direct runoff options in CFE

The user has the option to pick a particular direct runoff (aka surface partitioning) method:

1. Schaake function (configuration: `surface_partitioning_scheme=Schaake`)
2. Xinanjiang function (configuration: `surface_partitioning_scheme=Xinanjiang`). When using this runoff method the user must also include four parameters.

If the **Xinanjiang** scheme is choosen, four parameters need to be included in the configuration file:
1. a_Xinanjiang_inflection_point_parameter
2. b_Xinanjiang_shape_parameter
3. x_Xinanjiang_shape_parameter
4. urban_decimal_fraction 

## Surface runoff options in CFE
The user has the option to pick a particular surface runoff (aka surface runoff scheme) method:

1. GIUH-based surface runoff (configuration: `surface_runoff_scheme=GIUH`). This is the default option.
2. Nash_Cascade-based surface runoff (configuration: `surface_runoff_scheme=NASH_CASCADE`). In this method, GIUH is used to derive Nash cascade parameters K and N.


## Rootzone-based Actual Evapotranspiration (AET)
The user has the option to turn ON and OFF rootzone-based AET, default option is OFF. To turn it ON, the following parameters need to be included in the configuration file.
1. `is_aet_rootzone=true`
2. `soil_layer_depths` 
3. `max_rootzone_layer`

## CFE coupled to Soil freeze-thaw model (SFT)
The Soil Freeze-Thaw (SFT) model is a standalone model.  For detailed information please refer to the [SFT repo](https://github.com/NOAA-OWP/SoilFreezeThaw). A few things to note when coupling CFE to SFT:
1. SFT model provides `ice fraction` to CFE runoff schemes (Schaake `ice_fraction_schaake` and Xinanjiang `ice_fraction_xinanjiang`)
2. To turn ON/OFF SFT set sft_coupled flag.
    * `is_sft_coupled` : (type boolean) if `true`, SFT is turned ON. (options: True, true, 1).
    * If the runoff scheme is Xinanjiang, no additional parameters are needed in the CFE config files.
    * If the runoff scheme is Schaake, the CFE config file will need an additional parameter, namely:
      * `ice_content_threshold` : (type double, unit m). This represent the ice content above which soil is impermeable.
  

**Note:** By defualt `is_sft_coupled` and `is_aet_rootzone` are set to `OFF`, that means these changes do not affect the basic functionality of CFE.
