#ifndef _CONCEPTUAL_RESERVOIR_H
#define _CONCEPTUAL_RESERVOIR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

struct conceptual_reservoir {
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
  int    n_soil_layers;                 // number of soil layers
  double *soil_layer_depths_m;          // soil layer depths defined in the config file in units of [m]
  int    is_aet_rootzone;               // boolean - true if aet_root_zone is ON otherwise OFF (default is OFF)
  int    max_rootzone_layer;            // maximum root zone layer is used to identify the maximum layer to remove water from for AET
  double *delta_soil_layer_depth_m;     // used to calculate the total soil moisture in each layer
  double soil_water_content_field_capacity;  // water content [m/m] at field capacity.  Used in AET routine 
  
  //---------------------------------------------------------------
};

extern void conceptual_reservoir_flux_calc(struct conceptual_reservoir *da_reservoir,
                                           double *primary_flux_m, double *secondary_flux_m);

#endif
