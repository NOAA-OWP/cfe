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
 
 
#ifndef CFE_SOIL_DISCRETE_H
#define CFE_SOIL_DISCRETE_H

// Discrete soil moisture constants
#ifndef NDISC
#define NDISC 4
#endif

#ifndef NUM_LOOKUP_TABLE_VALS
#define NUM_LOOKUP_TABLE_VALS 5
#endif

#ifndef LOOKUP_TABLE_THETA_MIN
#define LOOKUP_TABLE_THETA_MIN 3.0e-02
#endif 

#ifndef MAX_LOOKUP_TABLE_POINTS
#define MAX_LOOKUP_TABLE_POINTS 500
#endif

#endif // CFE_SOIL_DISCRETE_H
