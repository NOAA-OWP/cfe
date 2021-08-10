# Conceptual Functional Equivalent Model
This model is designed to be a simplified model of the National Water Model, which is functionally equicalent.

# Compiling this code to run examples with a "pseudo" or "mini" framework.
## CFE Model reads its own forcing file
The BMI functionality was developed as a standalone module in C. To compile this code the developer used these steps:
1. `module load gnu/10.1.0`
2. `gcc -lm ./src/main.c ./src/cfe.c ./src/bmi_cfe.c -o run_bmi`
This should generate an executable called **run_cfe_bmi**. To run this executable you must pass the path to the corresponding configuration file: `./run_bmi ./configs/cat_58_bmi_config_cfe.txt`
Included in this repository is an environment file (env_cheyenne.sh), and a "make and run" file (make_and_run_bmi.sh), which will compile the code and run an example. If you are on the Cheyenne computer, or if you can modify these files to your machine, you can simply follow these two steps to run this code:
1. `source env_cheyenne.sh`
2. `./make_and_run_bmi.sh`
## CFE Model gets forcings passed from BMI
The CFE was designed to read its own forcing file, but we have added an option to get forcings passed in through BMI. To test this functionality we need to also include the AORC BMI model when compiling. The steps are very similar to the example above, but with just adding two additional src files, which come from the AORC BMI Module. Notice that below in step 2 the AORC files come from a different directory. 
1. `module load gnu/10.1.0`
2. `gcc -lm ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ../AORC/src/aorc.c ../AORC/src/bmi_aorc.c  -o run_cfe_bmi_pass_forcings`
This should generate an executable called **run_cfe_bmi_pass_forcings**. To run this executable you must pass the path to the corresponding configuration file for **BOTH** CFE and AORC (in that order): `./run_cfe_bmi_pass_forcings ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt`

# Making this code to run in the NGen Framework
make clean  
make (the resulting executable is moved to parent directory in the Makefile script)  
cd ..  
./run_cfe_bmi ./configs/cat_89_bmi_config_cfe.txt  


# The CFE was based on the t-short approximation of the National Water Model
# t-shirt approximation of the hydrologic routing funtionality of the National Water Model v 1.2, 2.0, and 2.1
This code was developed to test the hypothesis that the National Water Model runoff generation, vadose zone
dynamics, and conceptual groundwater model can be greatly simplified by acknowledging that it is truly a 
conceptual model. The hypothesis is supported by a number of observations made during a 2017-2018 deep dive
into the NWM code. Thesed are:
1. Rainfall/throughfall/melt partitioning in the NWM is based on a simple curve-number like approach that was developed by Schaake et al. (1996) and which is very similar to the Probability Distributed Moisture (PDM) function by Moore, 1985.   The Schaake function is a single valued function of soil moisture deficit, predicts 100% runoff when the soil is saturated, like the curve-number method, and is fundamentally simple.
2. Run-on infiltration is strictly not calculated.  Overland flow routing applies the Schaake function repeatedly to predict this phenomenon, which violates the underlying assumption of the PDM method that only rainfall  inputs affect soil moisture.
3. The water-content based Richards' equation, applied using a coarse-discretization, can be replaced with a simple conceptual reservoir because it never allows saturation or infiltration-excess runoff unless deactivated by assuming no-flow lower boundary condition.  Since this form of Richards' equation cannot simulate heterogeneous soil layers, it can be replaced with a conceptual reservoir.
4. The lateral flow routing function in the NWM is purely conceptual.  It is activated whenever the soil water content in one or more of the four Richards-equation discretizations reaches the wilting point water content. This activation threshold is physically unrealistic, because in most soils lateral subsurface flow is not active until pore water pressures become positive at some point in the soil profile.  Furthermore, the lateral flow hydraulic conductivity is assumed to be the vertical hydraulic conductivity multiplied by a calibration factor "LKSATFAC" which is allowed to vary between 10 and 10,000 during calibration, resulting in an anisotropy ratio that varies over the same range, without correlation with physiographic characteristics or other support.

This code implements these assumptions using pure conceptualizations.  The formulation consists of the following:
1. Rainfall is partitioned into direct runoff and soil moisture using the Schaake function.
2. Rainfall that becomes direct runoff is routed to the catchment outlet using a geomorphological instantanteous unit hydrograph (GIUH) approach, eliminating the 250 m NWM routing grid, and the incorrect use of the Schaake function to simulate run-on infiltration.
3. Water partitioned by the Schaake function to be soil moisture is placed into a conceptual linear reservoir that consists of two outlets that apply a minimum storage activation threshold.   This activation threshold is identical for both outlets, and is based on an integral solution of the storage in the soil assuming Clapp-Hornberger parameters equal to those used in the NWM to determine that storage corresponding to a soil water content 0.5 m above the soil column bottom that produces a soil suction head equal to -1/3 atm, which is a commonly applied assumption used to estimate the field capacity water content. The first outlet calculates vertical percolation of water to deep groundwater using the saturated hydraulic conductivity of the soil multiplied by the NWM "slope" parameter, which when 1.0 indicates free drainage and when 0.0 indicates a no-flow lower boundary condition. The second outlet is used to calculate the flux to the soil lateral flow path, using a conceptual LKSATFAC-like calibration parameter.
4. The lateral flow is routed to the catchment outlet using a Nash-cascade of reservoirs to produce a mass-conserving delayed response, and elminates the need for the 250 m lateral flow routing grid.
5. The groundwater contribution to base flow is modeled using either (a) an exponential nonlinear reservoir identical to the one in the NWM formulation, or (b) a nonlinear reservoir forumulation, which can also be made linear by assuming an exponent value equal to 1.0.

**This code was written entirely by Fred L. Ogden, May 22-24, 2020, in the service of the NOAA-NWS Office of Water Prediction, in Tuscaloosa, Alabama.**
