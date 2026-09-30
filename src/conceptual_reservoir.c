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
 
 
#include "conceptual_reservoir.h"
#include <math.h>


int is_epsilon_less_than(double a, double eps) {
    if(fabs(a)<eps) return 1;
    else            return 0;
}

//##############################################################
//########## SINGLE OUTLET EXPONENTIAL RESERVOIR ###############
//##########                -or-                 ###############
//#####    TWO OUTLET LINEAR OR NONLINEAR RESERVOIR   ##########
//##############################################################
// This function calculates the flux from a linear, or nonlinear 
// conceptual reservoir with one or two outlets, or from an
// exponential nonlinear conceptual reservoir with only one outlet.
// In the non-exponential instance, each outlet can have its own
// activation storage threshold.  Flow from the second outlet is 
// turned off by setting the discharge coeff. to 0.0.
//################################################################

extern void conceptual_reservoir_flux_calc(struct CONCEPTUAL_RESERVOIR_STRUCTURE *da_reservoir,
                                           double *primary_flux_m, double *secondary_flux_m)
{
  //struct conceptual_reservoir  <<<<INCLUDED HERE FOR REFERENCE.>>>>
  //{
  // int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
  // double storage_max_m;
  // double storage_m;
  // double coeff_primary;
  // double exponent_primary;        // Fixed: was exponent_secondary in comment
  // double storage_threshold_primary_m;
  // double storage_threshold_secondary_m;
  // double coeff_secondary;
  // double exponent_secondary;
  // };
  
  // THIS FUNCTION CALCULATES THE FLUXES FROM A CONCEPTUAL NON-LINEAR (OR LINEAR) RESERVOIR WITH TWO OUTLETS
  // all fluxes calculated by this routine are instantaneous with units of the coefficient.
  
  // Define epsilon constants for better maintainability
  const double EPSILON_SMALL = 1.0e-06;
  const double EPSILON_TINY = 1.0e-07;
  const double EPSILON_EXPONENT = 1.0e-03;
  
  // Initialize output parameters
  *primary_flux_m = 0.0;
  *secondary_flux_m = 0.0;
  
  // *****************************************************************************
  // ------------------ Conceptual Ground Water Reservoir -------------------------
  // single outlet reservoir like the NWM V1.2 exponential conceptual gw reservoir
  
  if (da_reservoir->is_exponential == TRUE && is_epsilon_less_than(da_reservoir->coeff_secondary, EPSILON_SMALL)) {
    // ------------------------------------------------------------------------------------------  
    // The code goes here to mimic the exponential conceptual gw nonlinear reservoir in WRF-Hydro
    // It calculates the flux to streamflow and returns.
    // ------------------------------------------------------------------------------------------
    double exp_term = exp(da_reservoir->exponent_primary * da_reservoir->storage_m / da_reservoir->storage_max_m);
    
    *primary_flux_m = da_reservoir->coeff_primary * (exp_term - 1.0);
    
    return;
  }
  
  // *****************************************************************************
  
  // code goes past here iff it is not a single outlet exponential deep groundwater reservoir of the NWM variety
  // The vertical outlet is assumed to be primary (percolation) and satisfied first.
  
  // ----------------------------------------- use shorter variable names to clean appearance of the calculations
  double S = da_reservoir->storage_m;
  double S_max = da_reservoir->storage_max_m;
  double T_pri = da_reservoir->storage_threshold_primary_m;
  double T_sec = da_reservoir->storage_threshold_secondary_m;
  
  // Check if any outflow is possible from either outlet using the minimum threshold
  double min_threshold;
  if (T_pri < T_sec) {
    min_threshold = T_pri;
  } else {
    min_threshold = T_sec;
  }
  
  if (S < min_threshold || is_epsilon_less_than(S - min_threshold, EPSILON_SMALL)) {
    return;  // no outflow possible from either outlet
  }
  
  // Calculate outflow from primary outlet
  double S_excess = S - T_pri;
  
  if (S_excess > EPSILON_SMALL) {  // primary outlet can flow
    if (is_epsilon_less_than(fabs(da_reservoir->exponent_primary - 1.0), EPSILON_EXPONENT)) {
      // it's essentially 1.0, so treat it like a linear reservoir
      *primary_flux_m = da_reservoir->coeff_primary * S_excess / (S_max - T_pri);
    } else {
      // it's nonlinear
      *primary_flux_m = da_reservoir->coeff_primary * pow((S_excess / (S_max - T_pri)), da_reservoir->exponent_primary);
    }
    
    // Limit primary flux to available storage above threshold
    if (*primary_flux_m > S_excess) {
      *primary_flux_m = S_excess;
    }
  }
  
  // Check if secondary flux calculation is needed
  if (is_epsilon_less_than(da_reservoir->coeff_secondary, EPSILON_TINY) || S < T_sec) {
    return;  // no secondary flux to calculate
  }
  
  // *****************************************************************************
  // Calculate secondary outlet flux
  // *****************************************************************************
  
  S_excess = S - T_sec;
  
  // Ensure we have storage above secondary threshold
  if (is_epsilon_less_than(S_excess, EPSILON_SMALL)) {
    return;
  }
  
  if (is_epsilon_less_than(fabs(da_reservoir->exponent_secondary - 1.0), EPSILON_EXPONENT)) {
    // it's essentially 1.0, so treat it like a linear reservoir
    *secondary_flux_m = da_reservoir->coeff_secondary * S_excess / (S_max - T_sec);
  } else {
    // it's nonlinear
    *secondary_flux_m = da_reservoir->coeff_secondary * pow((S_excess / (S_max - T_sec)), da_reservoir->exponent_secondary);
  }
  
  // Handle flux limiting - PRIMARY ALWAYS GETS PRIORITY!
  // The critical case is when thresholds are equal (common at field capacity)
  if (is_epsilon_less_than(T_pri - T_sec, EPSILON_SMALL)) {
    // Thresholds are essentially equal - both outlets compete for same storage pool
    // PRIMARY IS SATISFIED FIRST! Secondary gets remainder.
    // Use the lower threshold to be safe (in case of tiny numerical differences)
    
    double lower_thresh;
    if (T_pri < T_sec) {
      lower_thresh = T_pri;
    } else {
      lower_thresh = T_sec;
    }
    double total_available = S - lower_thresh;
    
    // Primary flux gets first claim on available storage
    if (*primary_flux_m > total_available) {
      if (total_available > 0.0) {
        *primary_flux_m = total_available;
      } else {
        *primary_flux_m = 0.0;
      }
    }
    
    // Secondary gets whatever is left after primary is satisfied
    double remaining = total_available - (*primary_flux_m);
    if (*secondary_flux_m > remaining) {
      if (remaining > 0.0) {
        *secondary_flux_m = remaining;
      } else {
        *secondary_flux_m = 0.0;
      }
    }
    
  } else if (T_sec < T_pri) {
    // Secondary threshold is lower - but PRIMARY STILL GETS PRIORITY
    double total_available = S - T_sec;
    
    // Primary flux gets priority
    if (*primary_flux_m > total_available) {
      if (total_available > 0.0) {
        *primary_flux_m = total_available;
      } else {
        *primary_flux_m = 0.0;
      }
    }
    
    // Secondary gets remaining storage after primary is satisfied
    double remaining = total_available - (*primary_flux_m);
    if (*secondary_flux_m > remaining) {
      if (remaining > 0.0) {
        *secondary_flux_m = remaining;
      } else {
        *secondary_flux_m = 0.0;
      }
    }
    
  } else {
    // Primary threshold is lower - secondary gets remaining storage after primary
    double remaining = S_excess - (*primary_flux_m);
    if (*secondary_flux_m > remaining) {
      if (remaining > 0.0) {
        *secondary_flux_m = remaining;
      } else {
        *secondary_flux_m = 0.0;
      }
    }
  }
  
  return;
}


