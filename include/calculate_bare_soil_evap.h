/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1
 *                and earlier
 *
 * Originally conceived and developed by:
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS
 *         Office of Water Prediction, Tuscaloosa, AL
 */
#ifndef CALCULATE_BARE_SOIL_EVAP_H
#define CALCULATE_BARE_SOIL_EVAP_H

#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

double calculate_bare_soil_evap(
    const cfe_forcing_struct *forcing,
    const cfe_parameters_struct *parameters,
    const cfe_state_struct *state,
    int dt_seconds,
    double aerodynamic_resistance_s_per_m,
    double rsurf_exp,
    double *surface_resistance_s_per_m_out,
    double *surface_relative_humidity_out);

#ifdef __cplusplus
}
#endif

#endif
