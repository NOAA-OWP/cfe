/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1 
 *                and earlier
 *
 * Originally conceived and developed by: 
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS 
 *         Office of Water Prediction, Tuscaloosa, AL
 *
 */


#ifndef _NASH_H
#define _NASH_H

#include <math.h>

#define MAX_NUM_SURFACE_NASH_CASCADE  6  
#define MAX_NUM_SUBSURFACE_NASH_CASCADE 2 

// this data structure describes runoff using Nash Cascade model in the subsurface and on the surface
struct NASH_CASCADE_PARAMETERS_STRUCTURE {
  int    N_nash;                // Number of Nash cascade reservoirs; [-]
  double K_nash;                // Fraction of storage per hour that moves from one reservoir to the next (time constant); [1/hour]
  int    nsubsteps;             // the number of substeps that each dt is divided into
  double *nash_storage;         // storage array nash cascade reservoirs [m]

  int    is_riparian_gw;        // flag to turn on/off riparian groundwater (currently used in CASAM only)
};

double nash_cascade_routing(double runoff_m, double soil_storage_deficit_m,
			    struct NASH_CASCADE_PARAMETERS_STRUCTURE *nash_params);

#endif
