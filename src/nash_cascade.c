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


#endif
