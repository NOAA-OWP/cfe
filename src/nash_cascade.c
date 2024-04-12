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
  // local vars
  double dt_h = 1.0; // model timestep [hour]
  double subdt = dt_h/nsubsteps;
  double dS = 0.0; // change in reservoir storage
  double Q_r; // discharge from reservoir
  double Q_out; // discharge at the outlet (the last reservoir)
  
  double outflow_m = 0.0;
  nash_params->nash_storage[0] += runoff_m;
  
  //for (int i=0; i<2; i++)
  //    printf("before: %lf \n",  nash_params->nash_storage[i]);
  // Loop through number of sub-timesteps
  for (int ts = 0; ts < nsubsteps; ts++) {

    //Loop through reservoirs
    for(int i = 0; i < N_nash; i++) {

      Q_r = K_nash * nash_params->nash_storage[i]; // flow from reservoir i to i+1
      dS  = Q_r * subdt;                       // storage change in reservoir i per subtimestep
      nash_params->nash_storage[i] -= dS;          // updated storage in reservoir i

      if(i < (N_nash-1))
	nash_params->nash_storage[i+1] += dS;
      else
	Q_out = Q_r;
      
    }
    
    outflow_m += Q_out * subdt;  // Q_r at the end of N_nash loop is the discharge at the outlet
    
  }

  // Return the flow output
  return (outflow_m);

}

#endif
