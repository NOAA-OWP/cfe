#ifndef CFE_SOIL_TEMPERATURE_H
#define CFE_SOIL_TEMPERATURE_H

#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void initialize_soil_temperature_state(
    cfe_pet_temperature_state_struct *temperature_state);

int update_soil_skin_temperature_state(
    const cfe_forcing_struct *forcing,
    int dt_seconds,
    int day_of_year,
    cfe_state_struct *state);

#ifdef __cplusplus
}
#endif

#endif /* CFE_SOIL_TEMPERATURE_H */
