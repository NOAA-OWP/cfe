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
 
#include <stdio.h>
#include "nash_cascade.h"


//##############################################################
//#################  NASH CASCADE ROUTING   ####################  Note: was named nash_cascade_surface() before refactor
//##############################################################        by FLO 6/25 to make surface and subsurface routing both use this function
double nash_cascade_routing(double runoff_m, double soil_storage_deficit_m,
			    struct NASH_CASCADE_PARAMETERS_STRUCTURE *nash_params)
{
  //##############################################################
  // Solve ffor the flow through the Nash cascade to delay the
  // arrival of the lateral flow into the channel
  //##############################################################

  int nsubsteps  = nash_params->nsubsteps;
  int N_nash     = nash_params->N_nash;
  double K_nash  = nash_params->K_nash;
  

  // local vars
  double dt_h  = 1.0;             // model timestep [hour]
  double subdt = dt_h/nsubsteps;
  double S     = 0.0;
  double dS    = 0.0;            // change in reservoir storage
  double Q_r;                    // discharge from reservoir
  double Q_out = 0.0;            // discharge at the outlet (the last reservoir) per subtimestep
  double Q_to_channel_m = 0.0;   // total outflow to channel per timestep
  
  nash_params->nash_storage[0] += runoff_m;

  // Loop through number of sub-timesteps
  for (int ts = 0; ts < nsubsteps; ts++) {

    //Loop through reservoirs (N_nash bounded by caller; no inner bounds check needed)
    for(int i = 0; i < N_nash; i++) {

      // if storage of ith reservoir is zero, move to the next reservoir
      if (nash_params->nash_storage[i] == 0.0)
	continue;

      // Route water through Nash reservoirs
      S = nash_params->nash_storage[i];


      Q_r = K_nash * S;                    // flow from reservoir i to i+1
      dS  = fmin(Q_r * subdt, S);          // storage change in reservoir i per subtimestep
      nash_params->nash_storage[i] -= dS;  // updated storage in reservoir i

      if(i < (N_nash-1))
        nash_params->nash_storage[i+1] += dS;
      else
        Q_out = Q_r;
      
    }

   Q_to_channel_m += Q_out * subdt;  // Q_r at the end of N_nash loop is the discharge at the outlet

  }

  // Return the flow output
  return (Q_to_channel_m);

}

