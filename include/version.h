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

#ifndef _VERSION_H
#define _VERSION_H

/* Integer-based semantic version for preprocessor comparisons.
 * Use CFE_VERSION_STRING for display; use CFE_VERSION_MAJOR etc. for #if guards. */
#define CFE_VERSION_MAJOR 3
#define CFE_VERSION_MINOR 0
#define CFE_VERSION_PATCH 0
#define CFE_VERSION_STRING "3.0.0-beta"

/* Float version — used by config file version detection and hotstart output.
 * Prefer the integer macros above for preprocessor guards. */
#define CFE_VERSION 3.00
#define CFE_SUBVERSION_STRING "beta"

#endif
