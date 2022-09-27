# Conceptual Functional Equivalent (CFE) Model

This model is designed to be a simplified (conceptual) model of the National Water Model, which is functionally equivalent.  CFE is a model written by Fred Ogden and was not originally BMI compatible.  To see the original author code, please refer to the `original_author_code` directory.  For more information on the hypotheses and ideas underpinning the CFE model, see the [T-shirt Approximation of the National Water Model versions 1.2, 2.0, and 2.1](https://github.com/NOAA-OWP/cfe/edit/AET_rootzone/README.md#t-shirt-approximation-of-the-national-water-model-versions-12-20-and-21) section of this document.  The remainder of this document discusses the BMI enabled and expanded CFE model. 


There are multiple ways to run CFE:

1. Through the BMI commands. There are script options, detailed below, for running the model in standalone mode and BMI allows the Next Generation Water Resources Modeling Framework to run CFE and couple it to surface routines. See [Compiling and running CFE](#compiling-and-running-cfe)
2. As written by the original author. This includes a full program to read and process atmospheric forcing data, print the model output and check for mass balance closure. This code can be run from the [original_author_code](./original_author_code) directory. This code does not have a BMI implementation.

## Compiling and Running CFE
There are four examples for running CFE as described below. They assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. To build without cmake see section `Alternative: Compiling and Running CFE` below.

1. `Option BASE` : CFE reads local forcing data (standalone CFE BMI run; one BMI example)
2. `Option FORCING` : CFE uses an external module to read in a forcing file and pass those data using BMI (two BMIs example)
3. `Option FORCINGPET` : CFE uses external modules to read in a forcing file and calculate potential evapotranspiration, and pass those data using BMIs (three BMIs example)
4. `Option AETROOTZONE` : Option #3 running with rootzone-based actual evapotranspiration (AET). This version of CFE is coupled to the [SoilMoistureProfiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) module which provides soil moisture (1D array) that is used to estimate actual evapotranspiration (AET) from the deepest rootzone layer. This example requires the user to clone [SoilMoistureProfiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) repo. (four BMIs example)

````
git clone https://github.com/NOAA-OWP/cfe.git
cd cfe
git checkout ajk/sft_aet_giuh_merge (update after PR merge)
git clone https://github.com/NOAA-OWP/SoilMoistureProfiles.git smc_coupler (needed if AETROOTZONE=ON)
mkdir build && cd build
cmake ../ [-DBASE=ON,-DFORCING=ON,-DFORCINGPET=ON,-DAETROOTZONE=ON] (pick one option, e.g. `cmake ../ -DFORCING=ON`)
make
cd ..
run_cfe.sh [BASE, FORCING, FORCINGPET, AETROOTZONE] (pick one option) 
````

## Note for running cfe in the ngen framework
Follow general build instructions [here](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial) and then build [SLoTH](https://github.com/NOAA-OWP/SLoTH) using the following instructions
```
cd extern/sloth/ && git checkout latest 
git submodule update --init --recursive
cd ../..
cmake -B extern/sloth/cmake_build -S extern/sloth/
make -C extern/sloth/cmake_build
```

## Note for customized examples (examples different than the above four)
The configuration files must be passed in this order: (1) the CFE configuration file, (2) the forcing configuration file, (3) the potential evapotranspiration (PET) configuration file, and (4) the soil moisture profile configuration file


## Configuration File
A [configs/](./configs/) directory contains primiary configuration text files for three different catchments pertaining to each process identiified in   The table below details information for [catchment-87](./configs/cat_87_bmi_config_cfe.txt).
| Variable | Datatype |  Limits  | Units | Role | Process | Description |
| -------- | -------- | ------ | ----- | ---- | ------- | ----------- |
| forcing_file | *char* | 256  |   | filename |   | path to forcing inputs csv; set to `BMI` if passed via `bmi.set_value*()`  |
| soil_params.depth | *double* |   | meters [m]| state |  | soil depth  |
| soil_params.b | *double* |   |   | state |   | beta exponent on Clapp-Hornberger (1978) soil water relations  |
| soil_params.satdk | *double* |   |  meters/second [m s-1] | state |  | saturated hydraulic conductivity  |
| soil_params.satpsi  | *double* |   |  meters [m] | state |  | saturated capillary head  |
| soil_params.slop   | *double* |   |  meters/meters [m/m]| state |  | this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow. |
| soil_params.smcmax  | *double* |   |  meters/meters [m/m] | state |  | saturated soil moisture content  |
| soil_params.wltsmc | *double* |   |  meters/meters [m/m] | state |   | wilting point soil moisture content  |
| soil_params.expon  | *double* |   |  | parameter_adjustable |    | optional; defaults to `1.0`  |
| soil_params.expon_secondary  | *double* |  |   | parameter_adjustable |  | optional; defaults to `1.0` |
| max_gw_storage | *double* |   |  meters [m] | parameter_adjustable |  | maximum storage in the conceptual reservoir |
| Cgw | *double* |   |  meters/hour [m h-1] | parameter_adjustable |  | the primary outlet coefficient |
| expon | *double* |   |   | parameter_adjustable |  | exponent parameter (1.0 for linear reservoir) |
| gw_storage | *double* |   |  meters/meters [m/m] | parameter_adjustable |  | initial condition for groundwater reservoir - it is the ground water as a decimal fraction of the maximum groundwater storage (max_gw_storage) for the initial timestep |
| alpha_fc | *double* |   |   | parameter_adjustable |  | field capacity |
| soil_storage| *double* |   | meters/meters [m/m] | parameter_adjustable |  | initial condition for soil reservoir - it is the water in the soil as a decimal fraction of maximum soil water storage (smcmax * depth) for the initial timestep |
| K_nash | *int* |   |   | parameter_adjustable |   | number of Nash lf reservoirs (optional, defaults to 2, ignored if storage values present)  |
| K_lf | *double* |   |   | parameter_adjustable |  | Nash Config param - primary reservoir  |
| nash_storage | *double* |   |   | parameter_adjustable |  | Nash Config param - secondary reservoir   |
| giuh_ordinates   | *double* |   |   | parameter_adjustable |  | Giuh ordinates in dt time steps   |
| num_timesteps  | *int* |   |  | time_info |  | set to `1` if `forcing_file=BMI`   |
| verbosity | *int* | `0`-`3`  |   | option |   |  prints various debug and bmi info  |
| surface_partitioning_scheme | *char* | `Xinanjiang` or `Schaake`  |  | parameter_adjustable | direct runoff |    |
| a_Xinanjiang_inflection_point_parameter | *double* |   |  | parameter_adjustable | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| b_Xinanjiang_shape_parameter=1  | *double* |   |   | parameter_adjustable  | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| x_Xinanjiang_shape_parameter=1  | *double* |   |   | parameter_adjustable | direct runoff | when `surface_partitioning_scheme=Xinanjiang`   |
| aet_rootzone                    | *boolean* | True, true or 1  |  | coupling parameter | `rootzone-based AET` | when `CFE coupled to SoilMoistureProfile` |
| sft_coupled                     | *boolean* | True, true or 1  |  | coupling parameter | `ice-fraction based runoff` | when `CFE coupled to SoilFreezeThaw`|
## Alternative: Compiling and Running CFE
### 1. Read local forcing file
To compile and run CFE with locally read forcing data, run the following from the command line:

1. `gcc -lm ./src/main.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c -o cfe_base`. This will generate an executable called `cfe_base`.
2.  Then run the model with example forcing data: `./cfe_base ./configs/cat_58_bmi_config_cfe.txt`  


### 2. CFE Model gets forcings passed from BMI

CFE was designed to read its own forcing file, but we have added an option to get forcings passed in through BMI using its `set_value` functionality. To demonstrate this functionality we have included the BMI-enabled AORC forcing read module. Follow the steps below:  

1. `gcc -lm ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c  -o cfe_forcing`. This generates an executable called `cfe_forcing`. 
2. To run this executable you must pass the path to the corresponding configuration files for **BOTH** CFE and AORC (in that order): `./cfe_forcing ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt`

### 3. CFE Model gets forcings AND potential evapotranspiration passed from BMI
CFE can remove mass from the modeled system through evapotranspiration (directly from precipitation and from the soil using the Budyko function). Follow the steps below:  

1. `gcc -lm ./src/main_cfe_aorc_pet.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c -o cfe_forcingpet`. This generates an executable called `cfe_forcingpet`.
2. To run this executable you must pass the path to the corresponding configuration files for CFE, PET and AORC (in that order):  `./cfe_forcingpet ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt ./configs/cat_89_bmi_config_pet_pass.txt`

### 4. CFE rootzone-based example couples C and C++ modules and can't be built without cmake

**NOTE:** the configuration files must be passed in this order: (1) the CFE configuration file, (2) the forcing configuration file, (3) the potential evapotranspiration (PET) configuration file, and (4) the soil moisture profile configuration file

## Direct runoff options in CFE

The user has the option to pick a particular direct runoff (aka surface partitioning) method:

1. Schaake function (configuration: `surface_partitioning_scheme=Schaake`)
2. Xinanjiang function (configuration: `surface_partitioning_scheme=Xinanjiang`). When using this runoff method the user must also include three parameters.

If the **Xinanjiang** scheme is choosen, three parameters need to be included in the configuration file:
1. a_Xinanjiang_inflection_point_parameter
2. b_Xinanjiang_shape_parameter
3. x_Xinanjiang_shape_parameter 

## Rootzone-based Actual Evapotranspiration (AET)
The user has the option to turn ON and OFF rootzone-based AET, default option is OFF. To turn it ON, the following parameters need to be included in the configuration file.
1. aet_rootzone=true
2. soil_layer_depths=0.1,0.4,1.0,2.0 (this is just an example)
3. max_root_zone_layer=2

## Soil freeze-thaw model (SFT)
The Soil Freeze-Thaw (SFT) model is a standalone model.  For detailed information please refer to the [SFT repo](https://github.com/NOAA-OWP/SoilFreezeThaw). A few things to note when coupling SFT to CFE:
1. SFT model provides `ice fraction` to CFE runoff schemes (Schaake `ice_fraction_schaake` and Xinanjiang `ice_fraction_xinan`)
2. To turn ON/OFF SFT set sft_coupled flag.
    * `sft_coupled` : (type boolean) if `true`, SFT is turned ON. (options: True, true, 1).
    * If the runoff scheme is Xinanjiang, no additional parameters are needed in the CFE config files.
    * If the runoff scheme is Schaake, the CFE config file will need an additional parameter, namely:
      * `ice_content_threshold` : (type double, unit m). This represent the ice content above which soil is impermeable.
  

**Note:** By defualt `sft_coupled` and `aet_rootzone` are `OFF`, that means these changes do not effect the basic functionality of CFE.


# T-shirt Approximation of the National Water Model versions 1.2, 2.0, and 2.1
The CFE code is based on the "t-shirt" approximation of the National Water Model (versions 1.2, 2.0, and 2.1).  The t-shirt approximation was developed to test the hypothesis that the National Water Model runoff generation, vadose zone dynamics, and conceptual groundwater model can be greatly simplified by acknowledging that it is truly a conceptual model. The hypothesis is supported by a number of observations made during a 2017-2018 deep dive into the NWM code. These are:
1. Rainfall/throughfall/melt partitioning in the NWM is based on a simple curve-number like approach that was developed by Schaake et al. (1996) and which is very similar to the Probability Distributed Moisture (PDM) function by Moore, 1985.   The Schaake function is a single valued function of soil moisture deficit, predicts 100% runoff when the soil is saturated, like the curve-number method, and is fundamentally simple.
2. Run-on infiltration is strictly not calculated.  Overland flow routing applies the Schaake function repeatedly to predict this phenomenon, which violates the underlying assumption of the PDM method that only rainfall  inputs affect soil moisture.
3. The water-content based Richards' equation, applied using a coarse-discretization, can be replaced with a simple conceptual reservoir because it never allows saturation or infiltration-excess runoff unless deactivated by assuming no-flow lower boundary condition.  Since this form of Richards' equation cannot simulate heterogeneous soil layers, it can be replaced with a conceptual reservoir.
4. The lateral flow routing function in the NWM is purely conceptual.  It is activated whenever the soil water content in one or more of the four Richards-equation discretizations reaches the wilting point water content. This activation threshold is physically unrealistic, because in most soils lateral subsurface flow is not active until pore water pressures become positive at some point in the soil profile.  Furthermore, the lateral flow hydraulic conductivity is assumed to be the vertical hydraulic conductivity multiplied by a calibration factor "LKSATFAC" which is allowed to vary between 10 and 10,000 during calibration, resulting in an anisotropy ratio that varies over the same range, without correlation with physiographic characteristics or other support.

This code implements these assumptions using pure conceptualizations.  The formulation consists of the following:
1. Rainfall is partitioned into direct runoff and soil moisture using the Schaake function.
2. Rainfall that becomes direct runoff is routed to the catchment outlet using a geomorphological instantaneous unit hydrograph (GIUH) approach, eliminating the 250 m NWM routing grid, and the incorrect use of the Schaake function to simulate run-on infiltration.
3. Water partitioned by the Schaake function to be soil moisture is placed into a conceptual linear reservoir that consists of two outlets that apply a minimum storage activation threshold.   This activation threshold is identical for both outlets, and is based on an integral solution of the storage in the soil assuming Clapp-Hornberger parameters equal to those used in the NWM to determine that storage corresponding to a soil water content 0.5 m above the soil column bottom that produces a soil suction head equal to -1/3 atm, which is a commonly applied assumption used to estimate the field capacity water content. The first outlet calculates vertical percolation of water to deep groundwater using the saturated hydraulic conductivity of the soil multiplied by the NWM "slope" parameter, which when 1.0 indicates free drainage and when 0.0 indicates a no-flow lower boundary condition. The second outlet is used to calculate the flux to the soil lateral flow path, using a conceptual LKSATFAC-like calibration parameter.
4. The lateral flow is routed to the catchment outlet using a Nash-cascade of reservoirs to produce a mass-conserving delayed response, and eliminates the need for the 250 m lateral flow routing grid.
5. The groundwater contribution to base flow is modeled using either (a) an exponential nonlinear reservoir identical to the one in the NWM formulation, or (b) a nonlinear reservoir forumulation, which can also be made linear by assuming an exponent value equal to 1.0.

**The original author code was written entirely by Fred L. Ogden, May 22-24, 2020, in the service of the NOAA-NWS Office of Water Prediction, in Tuscaloosa, Alabama.**
