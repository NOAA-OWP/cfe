
# GIUH Generation

**Description**:  
Calculates the GIUH (Geomorphological Instantaneous Unit Hydrograph), extract the parameters for CFE from the hydrofabrics files and generates the CFE config files for the basins in the hydrofabrics

This code:

1) Download data from http://web.corral.tacc.utexas.edu/nfiedata/HAND/ for the desired HUC06 of interest
2) Allows the use of 10 or 30 meters DEM
3) Calculate a raster with the travel time. Three methods are available:
	3.1) 1 = constant velocity ; 
	3.2) 2= varying velocity Wong, 1997 - in this case the travel time is a function of rain rate, upstream area, roughness and slope. The CFE code is will need to be modified to use this function, since at this point only a fixed relationship is used, independently of the rain rate. The user can specify a fixed rain rate if this function is to be used with the current version of CFE. If this option is selected, and rain rate is not specified, 10 mm/hour is used. 
	3.3) 3= varying velocity for gully and overland flow (Wong, 1997) and constant velocity in the channel. Same as 2, but with constant velocity for the channel. 

To specify the parameters for gully, overland and channel velocity, as well as rainfall rate, see parameters in generate_travel_time_by_pixel.py. 
 
4) Calculates the travel time to the river network (flowpath file in the hydrofabrics)
5) For each sub-basin in the HUC06, calculate the GIUH which correspond to the histogram of travel time (bins are for each hour, as defined in CFE)
6) Extract calibrated soil and groundwater parameters from the NWM 2.1  
7) Generates the cat_XX_bmi_config_cfe.txt file needed to run CFE - the file name contains the ID of the sub-basin as per the hydrofabrics. 

The "others" folder contains functions to plot the outputs. 

# Dependencies

 This code was tested in linux

# Software Requirements:
1) TauDEM (which requires gdal, mpiexec,... see https://github.com/dtarb/TauDEM)
2) Python is required to generate TWI histogram and the width function per basin. Running this model requires python and the libraries listed in the environment file: environment.yml. To install anaconda see https://docs.anaconda.com/. To create an environment in anaconda, open the anaconda console and run: 

conda env create --file environment.yml
conda activate params

 	
Anaconda is not required, as long as all requirements listed in environment.yml are available in the python installation. 
 	
3) Curl to download the HAND DEM data

## Usage
1) Edit the workflow_hand_twi_giuh.env file. This file contains all parameters for the runs, including the HUC06 number, path to the hydrofabrics, and environmental variabels for TAUDEM and gdal
2) Run: 
	./workflow_hand_twi_giuh.sh 

# Data Requirements:
hydrofabrics (catchments.geojson,flowpaths.geojson,gwbucket-params-fullrouting.csv,soil-properties-fullrouting.csv). 
if methods 3.2 and 3.3 are to be used, a raster file with manning information is required. This file is not provided here due to the size (16GB)
The HUC06 of the area covered by the watershed (include in the "workflow_hand_twi_giuh.env" file)
If the polygons in the hydrofabric covers multiple HUC06, all HUC06 can be specified in the "workflow_hand_twi_giuh.env" file  

## Open source licensing info


## Credits and references


