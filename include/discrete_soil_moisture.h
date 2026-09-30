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


#ifndef SOIL_KERNEL_STATELESS_H
#define SOIL_KERNEL_STATELESS_H

#include <stdio.h>
#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void et_from_soil_discrete
    (
    const SoilControl*    soil_control,
    const SoilGeometry*   soil_geometry, 
    const SoilParameters* soil_parameters, 
    SoilStateIn*    soil_state,
    struct EVAPOTRANSPIRATION_STRUCTURE* et_struct
    );
    
int DSBM_step_one_hour_stateless(
    const SoilControl        *ctrl,
    const SoilGeometry       *geom,
    const SoilParameters     *par,
    const SoilLookupTables   *lut,                     // may be NULL to use analytic CH
    const SoilStateIn        *sin,
    struct EVAPOTRANSPIRATION_STRUCTURE *evap_struct,  // Pass the actual ET struct
    const SoilForcing        *forcing,
    SoilStateOut             *sout,
    SoilFluxes               *flux,
    TimestepSoilVolbal       *mb,
    FILE                     *debug_fptr);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOIL_KERNEL_STATELESS_H
