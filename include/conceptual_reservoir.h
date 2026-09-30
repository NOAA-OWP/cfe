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
 

#ifndef _CONCEPTUAL_RESERVOIR_H
#define _CONCEPTUAL_RESERVOIR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0



struct CONCEPTUAL_RESERVOIR_STRUCTURE {
  // this data structure describes a nonlinear reservoir having two outlets, one primary with an activation
  // threshold that may be zero, and a secondary outlet with a threshold that may be zero
  // this will also simulate a linear reservoir by setting the exponent parameter to 1.0 iff is_exponential==FALSE
  // iff is_exponential==TRUE, then it uses the exponential discharge function from the NWM V2.0 forumulation
  // as the primary discharge with a zero threshold, and does not calculate a secondary discharge.
  //--------------------------------------------------------------------------------------------------
  int    is_exponential;                // set this true TRUE to use the exponential form of the discharge equation
  double gw_storage;                    // Initial Storage - LKC: added since I need to keep track of it when changing parameters
  double storage_max_m;                 // maximum storage in this reservoir
  double storage_m;                     // state variable.
  double storage_change_m;              // storage change in the current step
  double coeff_primary;                 // the primary outlet
  double exponent_primary;
  double storage_threshold_primary_m;
  double storage_threshold_secondary_m;
  double coeff_secondary;
  double exponent_secondary;
  double ice_fraction_schaake;
  double ice_fraction_xinanjiang;
  int    is_sft_coupled;                // boolean - true if SFT is ON otherwise OFF (default is OFF)
  
  //---Root zone adjusted AET development -rlm -ajk -------------
  double *smc_profile;                  //soil moisture content profile
  int    n_soil_discs;                  // number of soil discretizations
  double *soil_disc_depths_m;           // soil discrete depths defined in the config file in units of [m]
  int    is_aet_rootzone;               // boolean - true if aet_root_zone is ON otherwise OFF (default is OFF)
  int    max_rootzone_disc;             // largest (deepest) disc containing roots
  double *delta_soil_disc_depth_m;      // used to calculate the total soil moisture in each discretization (disc)
  double soil_water_content_field_capacity;  // water content [m/m] at field capacity.  Used in AET routine 
  
  //---------------------------------------------------------------
};

extern void conceptual_reservoir_flux_calc(struct CONCEPTUAL_RESERVOIR_STRUCTURE *da_reservoir,
                                           double *primary_flux_m, double *secondary_flux_m);

#endif
