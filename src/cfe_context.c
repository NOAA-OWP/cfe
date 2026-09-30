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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfe_context.h"
#include "cfe_helpers.h"
#include "cfe_driver_utils.h"  // needed for cfe_initialize_volume_balance()
#include "cfe_soil_skin_temperature.h"
#include "cfe_pet_priestley_taylor.h"
#include "calculate_bare_soil_evap.h"

// This file contains the initial context setup to create a BMI model definition.  FLO 9/2025 with alot of help from claude.ai


int cfe_context_create_from_config(const char* cfg_path, CFE_Model_Context** out_ctx)
{
    if (cfg_path == NULL || out_ctx == NULL) return -1;

    CFE_Model_Context* ctx = (CFE_Model_Context*)calloc(1, sizeof(CFE_Model_Context));
    if (ctx == NULL) return -1;

    ctx->forcing.DLWRF_surface = CFE_DLWRF_SURFACE_UNINITIALIZED_SENTINEL;

    double version = read_cfe_config_version(cfg_path);
    if (fabs(version) < 1.0e-04) {
        fprintf(stderr,
            "WARNING: No cfe_config_version found in: %s\n"
            "         Add cfe_config_version=3.0 to the config file to suppress this warning.\n"
            "         If this is a legacy v2 config, use cfe_migrate_config to convert it.\n",
            cfg_path);
        version = 3.0;
    }

    if (parse_config_driver(cfg_path, version, &ctx->config, &ctx->parameters, &ctx->options) != 0) {
        free(ctx);
        return -1;
    }
    if (cfe_initialize(&ctx->parameters, &ctx->options, &ctx->state) != 0) {
        free(ctx);
        return -1;
    }

    if (ctx->config.state_pet_initialized) {
        ctx->state.pet_temperature_state.skin_temperature_k =
            ctx->config.state_skin_temperature_k;
        ctx->state.pet_temperature_state.upper_soil_temperature_k =
            ctx->config.state_upper_soil_temperature_k;
        ctx->state.pet_temperature_state.estimated_annual_air_temperature_k =
            ctx->config.state_estimated_annual_air_temperature_k;
        ctx->state.pet_temperature_state.air_temperature_time_integral_k_s =
            ctx->config.state_air_temperature_time_integral_k_s;
        ctx->state.pet_temperature_state.accumulated_time_s =
            ctx->config.state_accumulated_time_s;
        ctx->state.pet_temperature_state.initialized =
            ctx->config.state_pet_initialized;
    }

    // this statement allows the volbal structure to persist
    cfe_initialize_volume_balance(&ctx->parameters, &ctx->options, &ctx->state, &ctx->volbal);
    
    *out_ctx = ctx;
    return 0;
}

//######################
void cfe_context_destroy(CFE_Model_Context* ctx)
{
    if (ctx == NULL) return;
    if (ctx->serialized_state != NULL) {
        free(ctx->serialized_state);
        ctx->serialized_state = NULL;
        ctx->serialized_size = 0;
    }
    cfe_finalize(&ctx->state);
    free(ctx);
}

//################################### <- helper function needed in context_update()
double calculate_total_storage(const CFE_Model_Context* ctx) {
    double total = ctx->state.gw_storage_m;
    
    if (ctx->options.simulate_discrete_soil_moisture) {
        double dz[4] = {0.1, 0.3, 0.6, 1.0};
        for (int i = 0; i < 4; i++) {
            total += ctx->state.soil_discrete_storage_theta[i] * dz[i];
        }
    } else {
        total += ctx->state.soil_storage_m;
    }
    
    if (ctx->options.surface_routing_scheme == SURF_ROUTE_GIUH) {
        for (int i = 0; i < ctx->parameters.giuh_num_ordinates; i++) {
            total += ctx->state.giuh_queue_m[i];
        }
    }
    
    // Subsurface routing storage (Nash cascade — always MAX_NUM_SUBSURFACE_NASH_CASCADE = 2)
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        total += ctx->state.nash_subsurface_storage_m[i];
    }
    
    return total;
}
//####################
int cfe_update_pet_and_bare_soil_evap(
    cfe_forcing_struct *forcing,
    const cfe_parameters_struct *parameters,
    cfe_state_struct *state,
    const cfe_options_struct *options)
{
    if (forcing == NULL || parameters == NULL || state == NULL || options == NULL)
        return -1;

    if (options->enable_ET_Priestley_Taylor == TRUE ||
        options->simulate_soil_evaporation == TRUE) {
        update_soil_skin_temperature_state(
            forcing,
            options->time_step_seconds,
            forcing->day_of_year,
            state);
    }

    if (options->enable_ET_Priestley_Taylor == TRUE) {
        forcing->et_potential_m =
            calculate_pet_priestley_taylor(
                forcing,
                options->time_step_seconds,
                parameters->alpha_pt,
                state);
    }

    forcing->forest_pet_m = forcing->et_potential_m;
    forcing->bare_soil_aet_m = 0.0;

    if (options->simulate_discrete_soil_moisture == TRUE &&
        options->simulate_soil_evaporation == TRUE) {

        double local_bare_soil_evaporation_m;

        forcing->forest_pet_m =
            forcing->et_potential_m *
            parameters->catchment_vegetated_fraction;

        local_bare_soil_evaporation_m =
            calculate_bare_soil_evap(
                forcing,
                parameters,
                state,
                options->time_step_seconds,
                CFE_AERODYNAMIC_RESISTANCE_S_PER_M,
                parameters->bare_soil_rsurf_exp,
                NULL,
                NULL);

        forcing->bare_soil_aet_m =
            local_bare_soil_evaporation_m *
            parameters->catchment_bare_soil_fraction;
    }

    return 0;
}

int cfe_context_update(CFE_Model_Context* ctx)
{
    if (ctx == NULL) return -1;

    if (ctx->params_dirty) {
        cfe_resync_derived_params(&ctx->parameters, &ctx->options, &ctx->state);
        ctx->params_dirty = 0;
    }

    // Calculate storage at start of timestep
    ctx->timestep_storage_start_m = calculate_total_storage(ctx);
    
    // Store inputs ffor this timestep
    ctx->timestep_input_m = ctx->forcing.rainfall_depth_m;
    
    double dt = (double)ctx->options.time_step_seconds;

    if (!ctx->forcing.day_of_year_set_externally &&
        ctx->options.epoch_start_seconds > 0.0) {
        double current_epoch = ctx->options.epoch_start_seconds +
                               (double)ctx->state.current_time_step * dt;
        ctx->forcing.day_of_year = day_of_year_from_epoch(current_epoch);
    }

    cfe_update_pet_and_bare_soil_evap(
        &ctx->forcing, &ctx->parameters, &ctx->state, &ctx->options);

    if (cfe_step(&ctx->parameters, &ctx->options, &ctx->state, &ctx->forcing, dt,
        &ctx->last_outputs, &ctx->volbal) != 0) {
       return -1;
    }
    
    // Calculate storage at end of timestep
    ctx->timestep_storage_end_m = calculate_total_storage(ctx);
    
    // Calculate outputs ffor this timestep (ET + discharge)
    ctx->timestep_output_m = ctx->last_outputs.qout_m + 
                            (ctx->last_outputs.actual_et_m);
    
    return 0;
}


/* Setters/getters */

int cfe_context_set_rainfall_depth_m(CFE_Model_Context* ctx, double depth_m)
{
    if (ctx == NULL) return -1;
    ctx->forcing.rainfall_depth_m = depth_m;
    return 0;
}

int cfe_context_set_et_potential_m(CFE_Model_Context* ctx, double et_m)
{
    if (ctx == NULL) return -1;
    ctx->forcing.et_potential_m = et_m;
    return 0;
}

int cfe_context_get_outputs(const CFE_Model_Context* ctx, cfe_outputs_struct* out)
{
    if (ctx == NULL || out == NULL) return -1;
    *out = ctx->last_outputs;
    return 0;
}

int cfe_context_get_time_step_seconds(const CFE_Model_Context* ctx, int* out_dt_s)
{
    if (ctx == NULL || out_dt_s == NULL) return -1;
    *out_dt_s = ctx->options.time_step_seconds;
    return 0;
}

int cfe_context_get_current_step(const CFE_Model_Context* ctx, int* out_step)
{
    if (ctx == NULL || out_step == NULL) return -1;
    *out_step = ctx->state.current_time_step;
    return 0;
}
