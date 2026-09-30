#ifndef CFE_PET_PRIESTLEY_TAYLOR_H
#define CFE_PET_PRIESTLEY_TAYLOR_H

#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Calculate potential evapotranspiration using Priestley-Taylor.
 *
 * The soil temperature state must be advanced before this function is called.
 * This function reads forcing and state but does not modify hydrologic or
 * temperature state.  PET is returned in meters per timestep.
 */
double calculate_pet_priestley_taylor(
    const cfe_forcing_struct *forcing,
    int dt_seconds,
    double alpha_pt,
    const cfe_state_struct *state);

#ifdef __cplusplus
}
#endif

#endif /* CFE_PET_PRIESTLEY_TAYLOR_H */
