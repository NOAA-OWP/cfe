#ifndef _GIUH_C
#define _GIUH_C

#include "giuh.h"


//##############################################################
//############### GIUH CONVOLUTION INTEGRAL   ##################
//##############################################################
extern double giuh_convolution_integral(double runoff_m,int num_giuh_ordinates, 
					double *giuh_ordinates, double *runoff_queue_m_per_timestep)
{
  //##############################################################
  // This function solves the convolution integral involving N
  //  GIUH ordinates.
  //##############################################################
  double runoff_m_current_timestep;
  int N,i;
  
  N = num_giuh_ordinates;
  runoff_queue_m_per_timestep[N] = 0.0;
  
  for(i=0;i<N;i++)
    {
      runoff_queue_m_per_timestep[i] += giuh_ordinates[i]*runoff_m;
    }
  
  runoff_m_current_timestep = runoff_queue_m_per_timestep[0];
  
  for(i=1;i<=N;i++)  // shift all the entries in preperation for the next timestep
    {
      runoff_queue_m_per_timestep[i-1] = runoff_queue_m_per_timestep[i];
    }
  //runoff_queue_m_per_timestep[N-1]=0.0;
  
  return runoff_m_current_timestep;
}


#endif
