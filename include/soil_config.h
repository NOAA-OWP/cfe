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


#ifndef SOIL_CONFIG_H
#define SOIL_CONFIG_H

// Fictitious lower bound on water content (m3/m3).
#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

// Tiny epsilon to avoid divide-by-zero.
#ifndef SOIL_EPS
#define SOIL_EPS 1.0e-12
#endif

// Inline helper macro
#ifndef SOIL_INLINE
#define SOIL_INLINE static inline
#endif

#endif // SOIL_CONFIG_H
