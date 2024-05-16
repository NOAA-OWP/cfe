#ifndef _NASH_C
#define _NASH_C

#include "nash_cascade.h"

//##############################################################
//###################   NASH CASCADE   #########################
//##############################################################
extern double nash_cascade(double flux_lat_m,int num_lateral_flow_nash_reservoirs,
                           double K_nash,double *nash_storage_arr)
{
  //##############################################################
  // Solve for the flow through the Nash cascade to delay the
  // arrival of the lateral flow into the channel
  //##############################################################
  // local vars
  int i;
  double outflow_m;
  static double Q[MAX_NUM_NASH_CASCADE];

  //Loop through reservoirs
  for(i = 0; i < num_lateral_flow_nash_reservoirs; i++)
    {
      Q[i] = K_nash*nash_storage_arr[i];
      nash_storage_arr[i]  -= Q[i];

      if (i==0) nash_storage_arr[i] += flux_lat_m;
      else      nash_storage_arr[i] +=  Q[i-1];
    }

  /*  Get Qout */
  outflow_m = Q[num_lateral_flow_nash_reservoirs-1];

  //Return the flow output
  return (outflow_m);

}


//##############################################################
//##############   NASH CASCADE SURFACE RUNOFF  ################
//##############################################################
double nash_cascade_surface_runoff(double runoff_m, double soil_storage_deficit_m,
				   struct nash_cascade_parameters *nash_params)
{
  //##############################################################
  // Solve for the flow through the Nash cascade to delay the
  // arrival of the lateral flow into the channel
  //##############################################################

  int nsubsteps = nash_params->nsubsteps;
  int N_nash    = nash_params->N_nash;
  double K_nash = nash_params->K_nash;
  double depth_r = nash_params->retention_depth;
  double K_infil = nash_params->K_infiltration;
  
  // local vars
  double dt_h  = 1.0;             // model timestep [hour]
  double subdt = dt_h/nsubsteps;
  double S     = 0.0;
  double dS    = 0.0;            // change in reservoir storage
  double dS_infil = 0.0;         // change in reservoir storage due to infiltration
  double Q_r;                    // discharge from reservoir
  double Q_out;                  // discharge at the outlet (the last reservoir) per subtimestep
  double Q_infil;                // discharge from reservoirs to soil
  double Q_to_channel_m = 0.0;  // total outflow to channel per timestep
  double Q_to_soil_m    = 0.0;  // runon infiltration (losses from surface runoff to soil)
  double soil_deficit_m = soil_storage_deficit_m; // local variable to track the soil storage deficit
  double infil_m;
  
  nash_params->nash_storage[0] += runoff_m;

  // Loop through number of sub-timesteps
  for (int ts = 0; ts < nsubsteps; ts++) {

    //Loop through reservoirs
    for(int i = 0; i < N_nash; i++) {

      // First: infiltration capacity should be satisfied before routing water through Nash reservoirs

      // compute runon infiltration
      Q_infil = K_infil * nash_params->nash_storage[i];
      infil_m = Q_infil * subdt;

      // case 1: soil_deficit greater than infil_m, all water can infiltrate, and soil_deficit decreases
      // case 2: soil_deficit less than infil_m, portion of water can infiltrate, soil_deficit gets zero
      // case 3: soil_deficit is zero, no infiltration

      // determine the right amount of infiltrated water based on soil deficit
      if (soil_deficit_m > 0.0 && soil_deficit_m > infil_m) {
	soil_deficit_m -= infil_m;           // update soil deficit for the next subtimstep iteration
      }
      else if (soil_deficit_m > 0.0) {
	infil_m = soil_deficit_m;
	soil_deficit_m = 0.0; 
      }
      else {
	infil_m = 0.0; // if got here, then soil is fully saturated, so no infiltration
      }

      // update nash storage (subtract infiltrated water)
      if  (nash_params->nash_storage[i] >= infil_m)
        nash_params->nash_storage[i] = nash_params->nash_storage[i] - infil_m;
      else {
        infil_m = nash_params->nash_storage[i];
        nash_params->nash_storage[i] = 0.0;
      }


      Q_to_soil_m += infil_m;  // water volume that infiltrates per subtimestep from each reservoir

      /*=========================================================================================*/
      // Second: time to route water through Nash reservoirs
      S = nash_params->nash_storage[i];

      // determine the amount of surface water available for routing from the first reservoir
      if (i == 0 &&  S > depth_r)
        S -= depth_r;
      else if (i == 0)
        S = 0.0;

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

  nash_params->runon_infiltration = Q_to_soil_m;
  
  // Return the flow output
  return (Q_to_channel_m);

}

#endif
