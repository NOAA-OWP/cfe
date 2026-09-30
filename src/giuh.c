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
 

#include "giuh.h"


//##############################################################
//############### GIUH CONVOLUTION INTEGRAL   ##################
//##############################################################
extern double giuh_convolution_integral(double runoff_m,int num_giuh_ordinates, 
					double *giuh_ordinates, double *giuh_runoff_queue_m_per_timestep)
{
  //##############################################################
  // This function solves the convolution integral involving N
  //  GIUH ordinates.
  //##############################################################
  double runoff_m_current_timestep;
  int N,i;
  
  N = num_giuh_ordinates;
  giuh_runoff_queue_m_per_timestep[N] = 0.0;
  
  for(i=0;i<N;i++)
    {
      giuh_runoff_queue_m_per_timestep[i] += giuh_ordinates[i]*runoff_m;
    }
  
  runoff_m_current_timestep = giuh_runoff_queue_m_per_timestep[0];
  
  for(i=1;i<=N;i++)  // shift all the entries in preperation for the next timestep
    {
      giuh_runoff_queue_m_per_timestep[i-1] = giuh_runoff_queue_m_per_timestep[i];
    }
  
  return runoff_m_current_timestep;
}
