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


#ifndef _GIUH_H
#define _GIUH_H

#ifndef MAX_NUM_GIUH_ORDINATES
#define MAX_NUM_GIUH_ORDINATES 11   // Ordinarily there should never be more than this on NextGen hydrofabric catchments
#endif

extern double giuh_convolution_integral(double runoff_m, int num_giuh_ordinates, 
					double *giuh_ordinates, double *giuh_runoff_queue_m_per_timestep);

#endif
