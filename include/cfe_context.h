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
 

#ifndef CFE_CONTEXT_H
#define CFE_CONTEXT_H

#include "cfe_types.h"

//NOTE IF YOU ARE LOOKING FOR THE DEFINITION OF CFE_Context.h  it is in cfe_types.h

/* Create context from config (allocates context) */
int cfe_context_create_from_config(const char* cfg_path, CFE_Model_Context** out_ctx);

/* Destroy context (frees state allocations and the context) */
void cfe_context_destroy(CFE_Model_Context* ctx);

/* Run one step using the currently set forcing */
int cfe_context_update(CFE_Model_Context* ctx);

/* Per-field setters/getters for BMI-friendly mapping */
int cfe_context_set_rainfall_depth_m(CFE_Model_Context* ctx, double depth_m);
int cfe_context_set_et_potential_m(CFE_Model_Context* ctx, double et_m);
int cfe_context_get_outputs(const CFE_Model_Context* ctx, cfe_outputs_struct* out);
int cfe_context_get_time_step_seconds(const CFE_Model_Context* ctx, int* out_dt_s);
int cfe_context_get_current_step(const CFE_Model_Context* ctx, int* out_step);

/* Advance skin temperature, compute PET (if P-T enabled), and compute
 * bare-soil evaporation.  Called from both BMI and standalone paths. */
int cfe_update_pet_and_bare_soil_evap(
    cfe_forcing_struct *forcing,
    const cfe_parameters_struct *parameters,
    cfe_state_struct *state,
    const cfe_options_struct *options);

/* Sum all storage compartments (soil + gw + routing) */
double calculate_total_storage(const CFE_Model_Context* ctx);

#endif
