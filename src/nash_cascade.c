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
double nash_cascade_surface_runoff(double runoff_m, struct nash_cascade_parameters *nash_params)
{
  //##############################################################
  // Solve for the flow through the Nash cascade to delay the
  // arrival of the lateral flow into the channel
  //##############################################################

  int nsubsteps = nash_params->nsubsteps;
  int N_nash    = nash_params->N_nash;
  double K_nash = nash_params->K_nash;
  double depth_r = nash_params->retention_depth;
  //double K_infil = nash_params->K_infiltration;
  double c0 = nash_params->Kinf_c0;
  double c1 = nash_params->Kinf_c1;
  
  // local vars
  double dt_h  = 1.0;             // model timestep [hour]
  double subdt = dt_h/nsubsteps;
  double S     = 0.0;
  double dS    = 0.0;            // change in reservoir storage
  double dS_infil = 0.0;         // change in reservoir storage due to infiltration
  double Q_r;                    // discharge from reservoir
  double Q_out;                  // discharge at the outlet (the last reservoir) per subtimestep
  double Q_infil;                // discharge from reservoirs to soil
  double K_infil;                  // K infiltration = c0 + c1 * S
  double Q_to_channel_m = 0.0;  // total outflow to channel per timestep
  double Q_to_soil_m    = 0.0;  // runon infiltration (losses from surface runoff to soil)

  nash_params->nash_storage[0] += runoff_m;

  // for (int i=0; i<2; i++)
  //    printf("before: %lf \n",  nash_params->nash_storage[i]);
  // Loop through number of sub-timesteps
  for (int ts = 0; ts < nsubsteps; ts++) {

    //Loop through reservoirs
    for(int i = 0; i < N_nash; i++) {

      S = nash_params->nash_storage[i];

      // first reservoir
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

      //compute runon infiltration
      K_infil = c0 + c1 * nash_params->nash_storage[i];
      printf("N1: %lf %lf %lf %lf \n", c0, c1, K_infil, nash_params->nash_storage[i]);
      Q_infil = K_infil * nash_params->nash_storage[i];
      dS_infil = Q_infil * subdt;

      if  (nash_params->nash_storage[i] >= dS_infil)
        nash_params->nash_storage[i] = nash_params->nash_storage[i] - dS_infil;
      else {
        dS_infil = nash_params->nash_storage[i];
        nash_params->nash_storage[i] = 0.0;
      }

      Q_to_soil_m += dS_infil;  // water volume that infiltrates per subtimestep from each reservoir

    }

   Q_to_channel_m += Q_out * subdt;  // Q_r at the end of N_nash loop is the discharge at the outlet

  }

  // for (int i=0; i<2; i++)
  //    printf("after: %lf \n",  nash_params->nash_storage[i]);
  nash_params->runon_infiltration = Q_to_soil_m;
  //printf("Nash end = %lf %lf \n", nash_params->runon_infiltration, Q_to_channel_m);
  
  // Return the flow output
  return (Q_to_channel_m);

}

#endif
