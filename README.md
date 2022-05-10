# Conceptual Functional Equivalent (CFE) Model
This model is designed to be a simplified (conceptual) model of the National Water Model, which is functionally equivalent.  CFE is a model written by Fred Ogden and was not originally BMI compatible.  To see the original author code, please refer to the `original_author_code` directory.  For more information on the hypotheses and ideas underpinning the CFE model, see the [T-shirt Approximation of the National Water Model versions 1.2, 2.0, and 2.1](https://github.com/NOAA-OWP/cfe/edit/AET_rootzone/README.md#t-shirt-approximation-of-the-national-water-model-versions-12-20-and-21) section of this document.  The remainder of this document discusses the BMI enabled and expanded CFE model. 

# Running CFE
This version of CFE is coupled to the [Soil Moisture Profiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) module which is used to estimate actual evapotranspiration (AET) from the deepest rootzone layer.  Therefore, you will need both repos to run the CFE model.  To run CFE:

````
git clone https://github.com/NOAA-OWP/cfe.git
cd cfe
git checkout AET_rootzone
git clone https://github.com/NOAA-OWP/SoilMoistureProfiles.git smc_coupler
mkdir build && cd build
cmake ../
make
cd ..
build/cfe_smp ./configs/laramie_bmi_config_cfe_pass.txt ./configs/laramie_bmi_config_aorc.txt ./configs/laramie_bmi_config_pet_pass.txt ./configs/laramie_bmi_config_smc_coupler.txt
````

**NOTE:** the configuration files must be passed in this order: (1) the CFE configuration file, (2) the forcing configuration file, (3) the potential evapotranspiration (PET) configuration file, and (4) the soil moisture profile configuration file

# Options in CFE
## Direct runoff 
The user has the option to pick a particular direct runoff (aka surface partitioning) method:
1. Schaake function (configuration: `surface_partitioning_scheme=Schaake`)
2. Xinanjiang function (configuration: `surface_partitioning_scheme=Xinanjiang`). When using this runoff method the user must also include three parameters.

If the **Xinanjiang** scheme is choosen, three parameters need to be included in the configuration file:
1. a_Xinanjiang_inflection_point_parameter
2. b_Xinanjiang_shape_parameter
3. x_Xinanjiang_shape_parameter 

## Soil freeze-thaw model (SFT)
The Soil Freeze-Thaw (SFT) model is a standalone model.  For detailed information please refer to the [SFT repo](https://github.com/NOAA-OWP/SoilFreezeThaw).  A few things to note when coupling SFT to CFE:
1. SFT model provides `ice fraction` to CFE runoff schemes (Schaake `ice_fraction_schaake` and Xinanjiang `ice_fraction_xinan`)
2. If the runoff scheme is Xinanjiang, no additional parameters are needed in the CFE config files. The framework will pass the parameter `ice_fraction_xina` to CFE.
3. If the runoff scheme is Schaake, the CFE config file will need additional parameters, namely:
    * `ice_fraction` : (type boolean) if coupled this parameter needs to be `true` (options: True, true, 1). 
    * `ice_content_threshold` : (type double, unit m), must be specified when `ice_fraction=1`.
  
**Note:** See SFT-enabled sample configuration file [here](./configs/cat_87_bmi_config_cfe.txt). By defualt `ice fraction` is zero, and these changes do not effect the basic functionality of CFE when uncoupled.  


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
