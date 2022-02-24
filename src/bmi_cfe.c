#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"
#include <time.h>
#include <float.h>
#ifndef WATER_SPECIFIC_WEIGHT
#define WATER_SPECIFIC_WEIGHT 9810
#endif

#define CFE_DEGUG 0

#define INPUT_VAR_NAME_COUNT 2
#define OUTPUT_VAR_NAME_COUNT 6
#define STATE_VAR_NAME_COUNT 89   // must match var_info array size

#define PARAM_VAR_NAME_COUNT 10
static const char *param_var_names[PARAM_VAR_NAME_COUNT] = {
    "maxsmc", "satdk", "slope", "b", "multiplier", "Klf", 
    "Kn", "Cgw", "expon", "max_gw_storage"
};

static const char *param_var_types[PARAM_VAR_NAME_COUNT] = {
    "double", "double", "double", "double", "double", "double",
    "double", "double", "double", "double"
};
//----------------------------------------------
// Put variable info into a struct to simplify
// BMI implementation and avoid errors.
//---------------------------------------------
// Should we add "/0" after each string here?
// Everything works without it.
//---------------------------------------------
Variable var_info[] = {
	{ 0,  "timestep_rainfall_input_m",        "double", 1 },
	{ 1,  "soil_reservoir_storage_deficit_m", "double", 1 },
	{ 2,  "infiltration_depth_m",             "double", 1 },
	{ 3,  "gw_reservoir_storage_deficit_m",   "double", 1 },
	{ 4,  "timestep_h",                       "double", 1 },
	//--------------------------------
	// Vars in soil reservoir struct
	// type: conceptual_reservoir
	//--------------------------------
	{ 5,  "is_exponential",                   "int",    1 },
	{ 6,  "storage_max_m",                    "double", 1 },
	{ 7,  "storage_m",                        "double", 1 },
	{ 8,  "coeff_primary",                    "double", 1 },
	{ 9,  "exponent_primary",                 "double", 1 },
	{ 10, "storage_threshold_primary_m",      "double", 1 },
	{ 11, "storage_threshold_secondary_m",    "double", 1 },
	{ 12, "coeff_secondary",                  "double", 1 },
	{ 13, "exponent_secondary",               "double", 1 },
	//------------------------------
	// Vars in gw reservoir struct
	// type: conceptual_reservoir
	//------------------------------
	{ 14, "is_exponential",                   "int",    1 },
	{ 15, "storage_max_m",                    "double", 1 },
	{ 16, "storage_m",                        "double", 1 },
	{ 17, "coeff_primary",                    "double", 1 },
	{ 18, "exponent_primary",                 "double", 1 },
	{ 19, "storage_threshold_primary_m",      "double", 1 },
	{ 20, "storage_threshold_secondary_m",    "double", 1 },
	{ 21, "coeff_secondary",                  "double", 1 },
	{ 22, "exponent_secondary",               "double", 1 },
	//---------------------------------
	// Vars in NWM_soil_params struct
	//---------------------------------
	{ 23, "smcmax",          "double", 1 },
	{ 24, "wltsmc",          "double", 1 },
	{ 25, "satdk",           "double", 1 },
	{ 26, "satpsi",          "double", 1 },
	{ 27, "bb",              "double", 1 },
	{ 28, "slop",            "double", 1 },
	{ 29, "D",               "double", 1 },
	{ 30, "wilting_point_m", "double", 1 },
	//--------------------
	// Vars in et_struct
	//--------------------
	{ 31, "potential_et_m_per_s",        "double", 1 },
	{ 32, "potential_et_m_per_timestep", "double", 1 }, //TODO Doesn't appear to be used, but this is uninitialized
	{ 33, "actual_et_m_per_timestep",    "double", 1 },
	//------------------------------
	// Vars in vol_tracking_struct
	//------------------------------
	{ 34, "vol_sch_runoff",        "double", 1 },
	{ 35, "vol_sch_infilt",        "double", 1 },
	{ 36, "vol_to_soil",           "double", 1 },
	{ 37, "vol_to_gw",             "double", 1 },
	{ 38, "vol_soil_to_gw",        "double", 1 },
	{ 39, "vol_soil_to_lat_flow",  "double", 1 },
	{ 40, "volstart",              "double", 1 },
	{ 41, "volout",                "double", 1 },
	{ 42, "volin",                 "double", 1 },
	{ 43, "vol_from_gw",           "double", 1 },
	{ 44, "vol_out_giuh",          "double", 1 },
	{ 45, "vol_in_nash",           "double", 1 },
	{ 46, "vol_out_nash",          "double", 1 },
	{ 47, "vol_in_gw_start",       "double", 1 },
	{ 48, "vol_soil_start",        "double", 1 },
	//-----------------------------------------         
	// More top-level, static allocation vars
	//-----------------------------------------
	{ 49, "epoch_start_time",      "long", 1 },
	{ 50, "num_timesteps",         "int",  1 },
	{ 51, "current_time_step",     "int",  1 },
	{ 52, "time_step_size",        "int",  1 },
	{ 53, "is_forcing_from_bmi",   "int",  1 },
	{ 54, "forcing_file",                                 "string", 1 },  // strlen
	//{ 56, "Schaake_adjusted_magic_constant_by_soil_type", "double", 1 }, 
    	{ 55, "surface_partitioning_scheme", "int", 1 }, // from direct_runoff_params_struct
	{ 56, "num_lateral_flow_nash_reservoirs",             "int",    1 },
	{ 57, "K_lf",                                         "double", 1 },
	{ 58, "K_nash",                                       "double", 1 },
	{ 59, "num_giuh_ordinates",                           "int",    1 },
	//---------------------------------------
	// Vars in aorc_forcing_data_cfe struct
	//---------------------------------------
	{ 60, "precip_kg_per_m2",               "double", 1 },
	{ 61, "incoming_longwave_W_per_m2",     "double", 1 },
	{ 62, "incoming_shortwave_W_per_m2",    "double", 1 },
	{ 63, "surface_pressure_Pa",            "double", 1 },
	{ 64, "specific_humidity_2m_kg_per_kg", "double", 1 },
	{ 65, "air_temperature_2m_K",           "double", 1 },
	{ 66, "u_wind_speed_10m_m_per_s",       "double", 1 },
	{ 67, "v_wind_speed_10m_m_per_s",       "double", 1 },
	{ 68, "latitude",                       "double", 1 },
	{ 69, "longitude",                      "double", 1 },
	{ 70, "time",                           "long",   1 },
	//------------------------------------------         
	// More top-level, dynamic allocation vars
	// (all pointers except verbosity)
	//------------------------------------------
	{ 71, "forcing_data_precip_kg_per_m2",  "double*", 1 },
	{ 72, "forcing_data_time",              "long*",   1 },
	{ 73, "giuh_ordinates",                 "double*", 1 },  // num_giuh
	{ 74, "nash_storage",                   "double*", 1 },  // num_lat_flow
	{ 75, "runoff_queue_m_per_timestep",    "double*", 1 },  // num_giuh
	{ 76, "flux_Schaake_output_runoff_m",   "double*", 1 },
	{ 77, "flux_giuh_runoff_m",             "double*", 1 },
	{ 78, "flux_nash_lateral_runoff_m",     "double*", 1 },
	{ 79, "flux_from_deep_gw_to_chan_m",    "double*", 1 },
	{ 80, "flux_perc_m",                    "double*", 1 },
	{ 81, "flux_lat_m",                     "double*", 1 },
	{ 82, "flux_Qout_m",                    "double*", 1 },
	{ 83, "verbosity",                      "int",     1 },
    //---------------------------------------
    // direct_runoff_params_struct vars
    // xinanjiang or schaake flag [56]
    //---------------------------------------
    { 84, "Schaake_adjusted_magic_constant_by_soil_type",   "double", 1},
    { 85, "a_Xinanjiang_inflection_point_parameter",        "double", 1},
    { 86, "b_Xinanjiang_shape_parameter",                   "double", 1},
    { 87, "x_Xinanjiang_shape_parameter",                   "double", 1}
    //---------------------------------------
};

int i = 0;
int j = 0;

// Don't forget to update Get_value/Get_value_at_indices (and setter) implementation if these are adjusted
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
        "RAIN_RATE",
        /* xinanjiang_dev
        "SCHAAKE_OUTPUT_RUNOFF",*/
        "DIRECT_RUNOFF",
        "GIUH_RUNOFF",
        "NASH_LATERAL_RUNOFF",
        "DEEP_GW_TO_CHANNEL_FLUX",
        "Q_OUT"
};

static const char *output_var_types[OUTPUT_VAR_NAME_COUNT] = {
        "double",
        "double",
        "double",
        "double",
        "double",
        "double"
};

static const int output_var_item_count[OUTPUT_VAR_NAME_COUNT] = {
        1,
        1,
        1,
        1,
        1,
        1
};

static const char *output_var_units[OUTPUT_VAR_NAME_COUNT] = {
        "m",
        "m",
        "m",
        "m",
        "m",
        "m"
};

static const int output_var_grids[OUTPUT_VAR_NAME_COUNT] = {
        0,
        0,
        0,
        0,
        0,
        0
};

static const char *output_var_locations[OUTPUT_VAR_NAME_COUNT] = {
        "node",
        "node",
        "node",
        "node",
        "node",
        "node"
};

// Don't forget to update Get_value/Get_value_at_indices (and setter) implementation if these are adjusted
static const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
        "atmosphere_water__liquid_equivalent_precipitation_rate",
        "water_potential_evaporation_flux"
};

static const char *input_var_types[INPUT_VAR_NAME_COUNT] = {
        "double",
        "double"
};

static const char *input_var_units[INPUT_VAR_NAME_COUNT] = {
        "mm h-1", //"atmosphere_water__liquid_equivalent_precipitation_rate"
        "m s-1"   //"water_potential_evaporation_flux"
};

static const int input_var_item_count[INPUT_VAR_NAME_COUNT] = {
        1,
        1
};

static const char input_var_grids[INPUT_VAR_NAME_COUNT] = {
        0,
        0
};

static const char *input_var_locations[INPUT_VAR_NAME_COUNT] = {
        "node",
        "node"
};

static int Get_start_time (Bmi *self, double * time)
{
    *time = 0.0;
    return BMI_SUCCESS;
}


static int Get_end_time (Bmi *self, double * time)
{
/*    Get_start_time(self, time);
    *time += (((cfe_state_struct *) self->data)->num_timesteps * ((cfe_state_struct *) self->data)->time_step_size);
    return BMI_SUCCESS;*/

    // JG EDIT 09.22.2021
    cfe_state_struct *cfe;
    cfe = (cfe_state_struct *) self->data;
    Get_start_time(self, time);
    
    // see if forcings read in or via BMI (framework to set)
    //jmframe: the docs say that "If the model doesn’t define an end time, a large number"
    //         so even if it is BMI, we should still use num_timesteps
    //         BMI shouldn't always have a very large end time...
    if (cfe->is_forcing_from_bmi == TRUE){
        // if BMI, set to FLT_MAX macro via float.h
        // See https://bmi.readthedocs.io/en/latest/#get-end-time
        *time += FLT_MAX;
        return BMI_SUCCESS;
    }
    else {
        // otherwise, set via numsteps as usual
        *time += cfe->num_timesteps * cfe->time_step_size;
        return BMI_SUCCESS;
    }

    return BMI_FAILURE;
}

// TODO: document that this will get the size of the current time step (the getter can access the full array)
static int Get_time_step (Bmi *self, double * dt)
{
    *dt = ((cfe_state_struct *) self->data)->time_step_size;
    return BMI_SUCCESS;
}


static int Get_time_units (Bmi *self, char * units)
{
    strncpy (units, "s", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int Get_current_time (Bmi *self, double * time)
{
    Get_start_time(self, time);
#if CFE_DEGUG > 1
    printf("Current model time step: '%d'\n", ((cfe_state_struct *) self->data)->current_time_step);
#endif
    *time += (((cfe_state_struct *) self->data)->current_time_step * ((cfe_state_struct *) self->data)->time_step_size);
    return BMI_SUCCESS;
}

/** Count the number of values in a delimited string representing an array of values. */
static int count_delimited_values(char* string_val, char* delimiter)
{
    char *copy, *copy_to_free, *value;
    int count = 0;

    // Make duplicate to avoid changing original string
    // Then work on copy, but keep 2nd pointer to copy so that memory can be freed
    copy_to_free = copy = strdup(string_val);
    while ((value = strsep(&copy, delimiter)) != NULL)
        count++;
    free(copy_to_free);
    return count;
}

/*int read_init_config_cfe(const char* config_file, cfe_state_struct* model, double* alpha_fc, double* soil_storage,
                     int* is_soil_storage_ratio)
{*/
int read_init_config_cfe(const char* config_file, cfe_state_struct* model, double* alpha_fc, double* soil_storage)
{
    int config_line_count, max_config_line_length;
    // Note that this determines max line length including the ending return character, if present
    int count_result = read_file_line_counts_cfe(config_file, &config_line_count, &max_config_line_length);
    if (count_result == -1) {
        printf("Invalid config file '%s'", config_file);
        return BMI_FAILURE;
    }
#if CFE_DEGUG >= 1
    printf("Config file details - Line Count: %d | Max Line Length %d\n", config_line_count, max_config_line_length);
#endif

    FILE* fp = fopen(config_file, "r");
    if (fp == NULL)
        return BMI_FAILURE;

    // TODO: document config file format (<param_key>=<param_val>, where array values are comma delim strings)


    // TODO: things needed in config file:
    //  - forcing file name
    //  - refkdt (optional, defaults to 3.0)
    //  - soil params:
    //    // - D, or depth
    //    //  - bb, or b
    //    //  - mult, or multiplier
    //    //  - satdk
    //    //  - satpsi
    //    //  - slop, or slope
    //    //  - smcmax, or maxsmc
    //    //  - wltsmc
    //    - additional gw res params
    //    //  - max_gw_storage
    //    //  - Cgw
    //    //  - expon
    //    //  - starting S_gw (may be literal or ratio, control by checking for "%")
    //    - additionally lateral flow res params
    //    //  - alpha_fc
    //    //  - starting S_lf (may be literal or ratio, control by checking for "%")
    //    - number of Nash lf reservoirs (optional, defaults to 2, ignored if storage values present)
    //    - K_nash
    //    - initial Nash storage values (optional, defaults to 0.0 for all reservoirs according to number)
    //    - GIUH ordinates

    char config_line[max_config_line_length + 1];

    // TODO: may need to add other variables to track that everything that was required was properly set

    // Keep track of whether required values were set in config
    // TODO: do something more efficient, maybe using bitwise operations
    int is_forcing_file_set = FALSE;
    int is_soil_params__depth_set = FALSE;
    int is_soil_params__bb_set = FALSE;
    int is_soil_params__satdk_set = FALSE;
    int is_soil_params__satpsi_set = FALSE;
    int is_soil_params__slop_set = FALSE;
    int is_soil_params__smcmax_set = FALSE;
    int is_soil_params__wltsmc_set = FALSE;
    int is_soil_params__expon_set = FALSE;
    int is_soil_params__expon2_set = FALSE;
    int is_Cgw_set = FALSE;
    int is_expon_set = FALSE;
    int is_alpha_fc_set = FALSE;
    int is_soil_storage_set = FALSE;
    int is_K_nash_set = FALSE;
    int is_K_lf_set = FALSE;
    int is_num_timesteps_set = FALSE;
    int is_verbosity_set = FALSE;

    /* xinanjiang_dev*/
    int is_direct_runoff_method_set = FALSE;
    int is_a_Xinanjiang_inflection_point_parameter_set = FALSE;
    int is_b_Xinanjiang_shape_parameter_set = FALSE;
    int is_x_Xinanjiang_shape_parameter_set = FALSE;

    // Keep track these in particular, because the "true" storage value may be a ratio and need both storage and max
    int is_gw_max_set = FALSE;
    int is_gw_storage_set = FALSE;

    int is_giuh_originates_string_val_set = FALSE;

    // Default value
    double refkdt = 3.0;

    //int is_gw_storage_ratio = FALSE;
    //double gw_storage_literal;
    // Also keep track of Nash stuff and properly set at the end of reading the config file
    int num_nash_lf = 2;
    char* nash_storage_string_val;
    int is_nash_storage_string_val_set = FALSE;
    // Similarly as for Nash, track stuff for GIUH ordinates
    char* giuh_originates_string_val;


    // Additionally,

    for (i = 0; i < config_line_count; i++) {
        char *param_key, *param_value, *param_units;
        fgets(config_line, max_config_line_length + 1, fp);
#if CFE_DEGUG >= 3
        printf("Line value: ['%s']\n", config_line);
#endif
        char* config_line_ptr = config_line;
        config_line_ptr = strsep(&config_line_ptr, "\n");
        param_key = strsep(&config_line_ptr, "=");
        //param_value = strsep(&config_line_ptr, "=");
        param_value = strsep(&config_line_ptr, "[");
        param_units = strsep(&config_line_ptr, "]");

#if CFE_DEGUG >= 1
        printf("Config Value - Param: '%s' | Value: '%s' | Units: '%s'\n", param_key, param_value, param_units);
#endif
        //printf(" %s | Value: %s | Units: %s\n", param_key, param_value, param_units);

        if (strcmp(param_key, "forcing_file") == 0) {
            model->forcing_file = strdup(param_value);
            is_forcing_file_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "refkdt") == 0) {
            refkdt = strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "soil_params.D") == 0 || strcmp(param_key, "soil_params.depth") == 0) {
            model->NWM_soil_params.D = strtod(param_value, NULL);
            is_soil_params__depth_set = TRUE;
            // Check if units are present and print warning if missing from config file
            if ((param_units == NULL) || (strlen(param_units) < 1)) {
#if CFE_DEGUG >= 1            
                printf ("WARNING: [units] expected for '%s' in config file \n", param_key);
#endif
            }
            continue;
        }
        if (strcmp(param_key, "soil_params.bb") == 0 || strcmp(param_key, "soil_params.b") == 0) {
            model->NWM_soil_params.bb = strtod(param_value, NULL);
            is_soil_params__bb_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.satdk") == 0) {
            model->NWM_soil_params.satdk = strtod(param_value, NULL);
            is_soil_params__satdk_set = TRUE;
            // Check if units are present and print warning if missing from config file
            if ((param_units == NULL) || (strlen(param_units) < 1)) {
#if CFE_DEGUG >= 1            
                printf ("WARNING: [units] expected for '%s' in config file \n", param_key);
#endif
            }
            continue;
        }
        if (strcmp(param_key, "soil_params.satpsi") == 0) {
            model->NWM_soil_params.satpsi = strtod(param_value, NULL);
            is_soil_params__satpsi_set = TRUE;
            // Check if units are present and print warning if missing from config file
            if ((param_units == NULL) || (strlen(param_units) < 1)) {
#if CFE_DEGUG >= 1            
                printf ("WARNING: [units] expected for '%s' in config file \n", param_key);
#endif
            }
            continue;
        }
        if (strcmp(param_key, "soil_params.slope") == 0 || strcmp(param_key, "soil_params.slop") == 0) {
            model->NWM_soil_params.slop = strtod(param_value, NULL);
            is_soil_params__slop_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.smcmax") == 0 || strcmp(param_key, "soil_params.maxsmc") == 0) {
            model->NWM_soil_params.smcmax = strtod(param_value, NULL);
            is_soil_params__smcmax_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.wltsmc") == 0) {
            model->NWM_soil_params.wltsmc = strtod(param_value, NULL);
            is_soil_params__wltsmc_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.expon") == 0) {
            model->soil_reservoir.exponent_primary = strtod(param_value, NULL);
            is_soil_params__expon_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.expon_secondary") == 0) {
            model->soil_reservoir.exponent_secondary = strtod(param_value, NULL);
            is_soil_params__expon2_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "max_gw_storage") == 0) {
            model->gw_reservoir.storage_max_m = strtod(param_value, NULL);
            is_gw_max_set = TRUE;
            // Also set the true storage if storage was already read and was a ratio, and so we were waiting for this
/*            if (is_gw_storage_set == TRUE && is_gw_storage_ratio == TRUE) {
                model->gw_reservoir.storage_m = (gw_storage_literal / 100.0) * model->gw_reservoir.storage_max_m;
            }*/
            continue;
        }
        if (strcmp(param_key, "Cgw") == 0) {
            model->gw_reservoir.coeff_primary = strtod(param_value, NULL);
            is_Cgw_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "expon") == 0) {
            model->gw_reservoir.exponent_primary = strtod(param_value, NULL);
            is_expon_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "gw_storage") == 0) {
            model->gw_reservoir.storage_m = strtod(param_value, NULL);
            is_gw_storage_set = TRUE;
/*            char* trailing_chars;
            gw_storage_literal = strtod(param_value, &trailing_chars);
            if (strcmp(trailing_chars, "%") == 0) {
                is_gw_storage_ratio = TRUE;
            }
            // Can't set the struct value yet unless storage is non-ratio or max storage was already set
            if (is_gw_storage_ratio == FALSE) {
                model->gw_reservoir.storage_m = gw_storage_literal;
            }
            if (is_gw_storage_ratio == TRUE && is_gw_max_set == TRUE) {
                model->gw_reservoir.storage_m = (gw_storage_literal / 100.0) * model->gw_reservoir.storage_max_m;
            }*/
            continue;
        }
        if (strcmp(param_key, "alpha_fc") == 0) {
            *alpha_fc = strtod(param_value, NULL);
            is_alpha_fc_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_storage") == 0) {
/*            char* trailing_chars;
            double parsed_value = strtod(param_value, &trailing_chars);
            *is_soil_storage_ratio = strcmp(trailing_chars, "%") == 0 ? TRUE : FALSE;
            *soil_storage = *is_soil_storage_ratio == TRUE ? (parsed_value / 100.0) : parsed_value;*/
            *soil_storage = strtod(param_value, NULL);
            is_soil_storage_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "number_nash_reservoirs") == 0 || strcmp(param_key, "N_nash") == 0) {
            num_nash_lf = strtol(param_value, NULL, 10);
            continue;
        }
        if (strcmp(param_key, "K_nash") == 0) {
            model->K_nash = strtod(param_value, NULL);
            is_K_nash_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "K_lf") == 0) {
            model->K_lf = strtod(param_value, NULL);
            is_K_lf_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "nash_storage") == 0) {
            nash_storage_string_val = strdup(param_value);
            is_nash_storage_string_val_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "giuh_ordinates") == 0) {
#if CFE_DEGUG >= 1
            printf("Found configured GIUH ordinate values ('%s')\n", param_value);
#endif
            giuh_originates_string_val = strdup(param_value);
            is_giuh_originates_string_val_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "num_timesteps") == 0) {
            model->num_timesteps = strtol(param_value, NULL, 10);
            is_num_timesteps_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "verbosity") == 0) {
            model->verbosity = strtol(param_value, NULL, 10);
            is_verbosity_set = TRUE;
            continue;
        }
        /* xinanjiang_dev: Need the option to run either runoff method in the config file, 
        *//////////////////////////////////////////////////////////////////////////////
        if (strcmp(param_key, "surface_partitioning_scheme") == 0) {
            if (strcmp(param_value, "Schaake")==0 || strcmp(param_value, "schaake")==0 || strcmp(param_value,"1")==0 )
                model->direct_runoff_params_struct.surface_partitioning_scheme = Schaake;
            if (strcmp(param_value, "Xinanjiang")==0 || strcmp(param_value, "xinanjiang")==0 || strcmp(param_value,"2")==0)
                model->direct_runoff_params_struct.surface_partitioning_scheme = Xinanjiang;
            is_direct_runoff_method_set = TRUE;
        }
        if (model->direct_runoff_params_struct.surface_partitioning_scheme == Xinanjiang) {  //Check that logical statement is correct
            if (strcmp(param_key, "a_Xinanjiang_inflection_point_parameter") == 0){
                model->direct_runoff_params_struct.a_Xinanjiang_inflection_point_parameter = strtod(param_value, NULL);
                is_a_Xinanjiang_inflection_point_parameter_set = TRUE;
            }
            if (strcmp(param_key, "b_Xinanjiang_shape_parameter") == 0) {
                model->direct_runoff_params_struct.b_Xinanjiang_shape_parameter = strtod(param_value, NULL);
                is_b_Xinanjiang_shape_parameter_set = TRUE;
            }
            if (strcmp(param_key, "x_Xinanjiang_shape_parameter") == 0) {
                model->direct_runoff_params_struct.x_Xinanjiang_shape_parameter = strtod(param_value, NULL);
                is_x_Xinanjiang_shape_parameter_set = TRUE;
            }
        }
    }

    if (is_forcing_file_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'forcing_file' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__depth_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.depth' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__bb_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.bb' not found in config file\n");
#endif
        return BMI_FAILURE;
    }

    if (is_soil_params__satdk_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.satdk' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__satpsi_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.satpsi' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__slop_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.slop' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__smcmax_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.smcmax' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__wltsmc_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.wltsmc' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_params__expon_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.expon' not found in config file, defaulting to 1 (linear)\n");
#endif
        model->soil_reservoir.exponent_primary = 1.0;
        //is_soil_params__expon_set == TRUE;
        // Don't return BMI_FAILURE, this is a optional config
        //return BMI_FAILURE;
    }
    if (is_soil_params__expon2_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.expon_secondary' not found in config file, defaulting to 1 (linear)\n");
#endif
        model->soil_reservoir.exponent_secondary = 1.0;
        // Don't return BMI_FAILURE, this is a optional config
        //return BMI_FAILURE;
    }
    if (is_Cgw_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'Cgw' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_expon_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'expon' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_alpha_fc_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'alpha_fc' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_soil_storage_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_storage' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_K_nash_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'K_nash' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_K_lf_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'K_nash' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_gw_max_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'max_gw_storage' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_gw_storage_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'gw_storage' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_num_timesteps_set == FALSE && strcmp(model->forcing_file, "BMI")) {
#if CFE_DEGUG >= 1
        printf("Config param 'num_timesteps' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
    if (is_verbosity_set == FALSE) {
        printf("Config param 'verbosity' not found in config file\n");
        printf("setting verbosity to a high value\n");
        model->verbosity = 10;
        return BMI_FAILURE;
    }
    if (is_direct_runoff_method_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'direct_runoff_method' not found in config file\n");
#endif
        return BMI_FAILURE;
    }
/* xinanjiang_dev*/
    if(model->direct_runoff_params_struct.surface_partitioning_scheme == Xinanjiang){
        if (is_a_Xinanjiang_inflection_point_parameter_set == FALSE) {
#if CFE_DEGUG >= 1
            printf("Config param 'a_Xinanjiang_inflection_point_parameter' not found in config file\n");
#endif
            return BMI_FAILURE;
        }
        if (is_b_Xinanjiang_shape_parameter_set == FALSE) {
#if CFE_DEGUG >= 1
            printf("Config param 'b_Xinanjiang_shape_parameter' not found in config file\n");
#endif
            return BMI_FAILURE;
        }
        if (is_x_Xinanjiang_shape_parameter_set == FALSE) {
#if CFE_DEGUG >= 1
            printf("Config param 'x_Xinanjiang_shape_parameter' not found in config file\n");
#endif
            return BMI_FAILURE;
        }   
    }

    if(model->direct_runoff_params_struct.surface_partitioning_scheme == Schaake){
        model->direct_runoff_params_struct.Schaake_adjusted_magic_constant_by_soil_type = refkdt * model->NWM_soil_params.satdk / 0.000002;   
#if CFE_DEGUG >= 1
    printf("Schaake Magic Constant calculated\n");
#endif
    }

#if CFE_DEGUG >= 1
    printf("All CFE config params present\n");
#endif    
    // Used for parsing strings representing arrays of values below
    char *copy, *value;

    // Handle GIUH ordinates, bailing if they were not provided
    if (is_giuh_originates_string_val_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("GIUH ordinate string not set!\n");
#endif
        return BMI_FAILURE;
    }
#if CFE_DEGUG >= 1
    printf("GIUH ordinates string value found in config ('%s')\n", giuh_originates_string_val);
#endif
    model->num_giuh_ordinates = count_delimited_values(giuh_originates_string_val, ",");
#if CFE_DEGUG >= 1
    printf("Counted number of GIUH ordinates (%d)\n", model->num_giuh_ordinates);
#endif
    if (model->num_giuh_ordinates < 1)
        return BMI_FAILURE;

    model->giuh_ordinates = malloc(sizeof(double) * model->num_giuh_ordinates);
    // Work with copy of the string pointer to make sure the original pointer remains unchanged, so mem can be freed at end
    copy = giuh_originates_string_val;
    // Now iterate back through and get the values (this modifies the string, which is why we needed the full string copy above)
    int i = 0;
    while ((value = strsep(&copy, ",")) != NULL)
        model->giuh_ordinates[i++] = strtod(value, NULL);
    // Finally, free the original string memory
    free(giuh_originates_string_val);

    // Now handle the Nash storage array properly
    if (is_nash_storage_string_val_set == TRUE) {
        // First, when there are values, read how many there are, and have that override any set count value
        int value_count = count_delimited_values(nash_storage_string_val, ",");
        // TODO: consider adding a warning if value_count and N_nash (assuming it was read from the config and not default) disagree
        // Ignore the values if there are not enough, and use whatever was set, or defaults
        if (value_count < 2) {
            model->num_lateral_flow_nash_reservoirs = num_nash_lf;
            model->nash_storage = malloc(sizeof(double) * num_nash_lf);
            for (j = 0; j < num_nash_lf; j++)
                model->nash_storage[j] = 0.0;
        }
        else {
            model->num_lateral_flow_nash_reservoirs = value_count;
            model->nash_storage = malloc(sizeof(double) * value_count);
            // Work with copy the string pointer to make sure the original remains unchanged, so it can be freed at end
            copy = nash_storage_string_val;
            // Now iterate back through and get the values
            int k = 0;
            while ((value = strsep(&copy, ",")) != NULL)
                model->nash_storage[k++] = strtod(value, NULL);
        }
        // Make sure at the end to free this too, since it was a copy
        free(nash_storage_string_val);
    }
    // If Nash storage values weren't set, initialize them to 0.0
    else {
        model->num_lateral_flow_nash_reservoirs = num_nash_lf;
        model->nash_storage = malloc(sizeof(double) * num_nash_lf);
        for (j = 0; j < num_nash_lf; j++)
            model->nash_storage[j] = 0.0;
    }

#if CFE_DEGUG >= 1
    printf("Finished function parsing CFE config\n");
#endif

    return BMI_SUCCESS;
}


static int Initialize (Bmi *self, const char *file)
{
    //FIXME, we can use the input file to help imply "framework" support or "standalone"
    //an empty init file string indicates things will come from set_value???
    //what happens when both occur, that is we have a config file and framewrok
    //using set_value after init???

    //Consider enumeration of failure states and how that might look across multiple languages
    //integrating into a framework

    cfe_state_struct* cfe_bmi_data_ptr;

    if (!self || !file)
        return BMI_FAILURE;
    else
        cfe_bmi_data_ptr = (cfe_state_struct *) self->data;

    cfe_bmi_data_ptr->current_time_step = 0;

    double alpha_fc, /*max_soil_storage,*/ S_soil;
    //int is_S_soil_ratio;

    int config_read_result = read_init_config_cfe(file, cfe_bmi_data_ptr, &alpha_fc, &S_soil);
    if (config_read_result == BMI_FAILURE)
        return BMI_FAILURE;

    // time_step_size is set to 3600 in the "new_bmi_cfe" function.
    cfe_bmi_data_ptr->timestep_h = cfe_bmi_data_ptr->time_step_size / 3600.0;
     
    // JG NOTE: this is done in init_soil_reservoir
    //max_soil_storage = cfe_bmi_data_ptr->NWM_soil_params.D * cfe_bmi_data_ptr->NWM_soil_params.smcmax;


    /***********************************************************************
       JMFRAME: Moved these up before the read forcing line,
                Since we need them even if we don't read forcings from file.
    ************************************************************************/

    /* xinanjiang_dev
        changing the name to the more general "direct runoff"
    cfe_bmi_data_ptr->flux_Schaake_output_runoff_m = malloc(sizeof(double));*/
    cfe_bmi_data_ptr->flux_output_direct_runoff_m  = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_output_direct_runoff_m = 0.0;
    cfe_bmi_data_ptr->flux_Qout_m = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_Qout_m = 0.0;
    cfe_bmi_data_ptr->flux_from_deep_gw_to_chan_m = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_from_deep_gw_to_chan_m = 0.0;
    cfe_bmi_data_ptr->flux_giuh_runoff_m = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_giuh_runoff_m = 0.0;
    cfe_bmi_data_ptr->flux_lat_m = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_lat_m = 0.0;
    cfe_bmi_data_ptr->flux_nash_lateral_runoff_m = malloc(sizeof(double));
    *cfe_bmi_data_ptr->flux_nash_lateral_runoff_m = 0.0;
    cfe_bmi_data_ptr->flux_perc_m = malloc(sizeof(double));
    //This needs an initial value, it is used in computations in CFE and only set towards the end of the model
    //See issue #31
    *cfe_bmi_data_ptr->flux_perc_m = 0.0;

    /*******************************************************
       JMFRAME: Check to see where forcings come from
                The options as of May 2021 are:
                1. Read your own forcings from a file
                2. Get the forcings passed in through BMI
    *******************************************************/
    if (strcmp(cfe_bmi_data_ptr->forcing_file, "BMI") == 0){
        cfe_bmi_data_ptr->is_forcing_from_bmi = TRUE;
        cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2 = malloc(sizeof(double));
        *cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2 = 0.0;
        cfe_bmi_data_ptr->forcing_data_time = malloc(sizeof(long));
        *cfe_bmi_data_ptr->forcing_data_time = 0;
    }
    else
    {
        cfe_bmi_data_ptr->is_forcing_from_bmi = FALSE;
    
        // Figure out the number of lines first (also char count)
        int forcing_line_count, max_forcing_line_length;
        int count_result = read_file_line_counts_cfe(cfe_bmi_data_ptr->forcing_file, &forcing_line_count, &max_forcing_line_length);
        if (count_result == -1) {
            printf("Configured forcing file '%s' could not be opened for reading\n", cfe_bmi_data_ptr->forcing_file);
            return BMI_FAILURE;
        }
        if (forcing_line_count == 1) {
            printf("Invalid header-only forcing file '%s'\n", cfe_bmi_data_ptr->forcing_file);
            return BMI_FAILURE;
        }
        // Infer the number of time steps: assume a header, so equal to the number of lines minus 1
        cfe_bmi_data_ptr->num_timesteps = forcing_line_count - 1;
    
    #if CFE_DEGUG > 0
        printf("Counts - Lines: %d | Max Line: %d | Num Time Steps: %d\n", forcing_line_count, max_forcing_line_length,
               cfe_bmi_data_ptr->num_timesteps);
    #endif

        // Now initialize empty arrays that depend on number of time steps
        cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2 = malloc(sizeof(double) * (cfe_bmi_data_ptr->num_timesteps + 1));
        cfe_bmi_data_ptr->forcing_data_time = malloc(sizeof(long) * (cfe_bmi_data_ptr->num_timesteps + 1));
    
        // Now open it again to read the forcings
        FILE* ffp = fopen(cfe_bmi_data_ptr->forcing_file, "r");
        // Ensure still exists
        if (ffp == NULL) {
            printf("Forcing file '%s' disappeared!", cfe_bmi_data_ptr->forcing_file);
            return BMI_FAILURE;
        }
    
        // Read forcing file and parse forcings
        char line_str[max_forcing_line_length + 1];
        long year, month, day, hour, minute;
        double dsec;
        // First read the header line
        fgets(line_str, max_forcing_line_length + 1, ffp);
        
        aorc_forcing_data_cfe forcings;
        for (i = 0; i < cfe_bmi_data_ptr->num_timesteps; i++) {
            fgets(line_str, max_forcing_line_length + 1, ffp);  // read in a line of AORC data.
            parse_aorc_line_cfe(line_str, &year, &month, &day, &hour, &minute, &dsec, &forcings);
    #if CFE_DEGUG > 0
            printf("Forcing data: [%s]\n", line_str);
            printf("Forcing details - s_time: %ld | precip: %f\n", forcings.time, forcings.precip_kg_per_m2);
    #endif
            cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2[i] = forcings.precip_kg_per_m2; //* ((double)cfe_bmi_data_ptr->time_step_size);
            cfe_bmi_data_ptr->forcing_data_time[i] = forcings.time;

            // TODO: make sure some kind of conversion isn't needed for the rain rate data
            // assumed 1000 kg/m3 density of water.  This result is mm/h;
            //rain_rate[i] = (double) aorc_data.precip_kg_per_m2;
        }
    
        cfe_bmi_data_ptr->epoch_start_time = cfe_bmi_data_ptr->forcing_data_time[0];
    } // end if is_forcing_from_bmi
    
    // divide by 1000 to convert from mm/h to m w/ 1h timestep as per t-shirt_0.99f
    cfe_bmi_data_ptr->timestep_rainfall_input_m = cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2[0] / 1000;

    // Initialize the rest of the groundwater conceptual reservoir (some was done when reading in the config)
    cfe_bmi_data_ptr->gw_reservoir.is_exponential = TRUE;
    cfe_bmi_data_ptr->gw_reservoir.storage_threshold_primary_m = 0.0;    // 0.0 means no threshold applied
    cfe_bmi_data_ptr->gw_reservoir.storage_threshold_secondary_m = 0.0;  // 0.0 means no threshold applied
    cfe_bmi_data_ptr->gw_reservoir.coeff_secondary = 0.0;                // 0.0 means that secondary outlet is not applied
    cfe_bmi_data_ptr->gw_reservoir.exponent_secondary = 1.0;             // linear

    // Initialize soil conceptual reservoirs
    init_soil_reservoir(cfe_bmi_data_ptr, alpha_fc, S_soil);

    // Initialize the runoff queue to empty to start with
    cfe_bmi_data_ptr->runoff_queue_m_per_timestep = malloc(sizeof(double) * cfe_bmi_data_ptr->num_giuh_ordinates + 1);
    for (i = 0; i < cfe_bmi_data_ptr->num_giuh_ordinates + 1; i++)
        cfe_bmi_data_ptr->runoff_queue_m_per_timestep[i] = 0.0;


    // jmframe: initialize zero potential and actual evapotranspiration.
    //          Will subtract this from the soil moisture and infiltration mass.
    cfe_bmi_data_ptr->et_struct.potential_et_m_per_s = 0;
    cfe_bmi_data_ptr->et_struct.actual_et_m_per_timestep = 0;
    cfe_bmi_data_ptr->et_struct.potential_et_m_per_timestep = 0;
    cfe_bmi_data_ptr->et_struct.reduced_potential_et_m_per_timestep = 0;
    cfe_bmi_data_ptr->et_struct.actual_et_from_rain_m_per_timestep = 0;
    cfe_bmi_data_ptr->et_struct.actual_et_from_soil_m_per_timestep = 0;
    // Set all the mass balance trackers to zero.
    initialize_volume_trackers(cfe_bmi_data_ptr);

    return BMI_SUCCESS;
}


static int Update (Bmi *self)
{
    // TODO: look at how the time step size (in seconds) effects 'coeff_primary' for lat flow reservoir, and whether
    //  this dictates that, if time step size is not fixed, an adjustment needs to be made to reservoir on each call
    //  here according to the size of the next time step.

//    double current_time, end_time;
//    self->get_current_time(self, &current_time);
//    self->get_end_time(self, &end_time);
//    printf("end time: %lf\n", end_time);
//    if (current_time >= end_time) {
//        return BMI_FAILURE;
//    }
            
    cfe_state_struct* cfe_ptr = ((cfe_state_struct *) self->data);

    // Two modes to get forcing data... 0/FALSE) read from file, 1/TRUE) pass with bmi    
    if (cfe_ptr->is_forcing_from_bmi == TRUE)
      // BMI sets the precipitation to the aorc structure.
      // divide by 1000 to convert from mm/h to m w/ 1h timestep as per t-shirt_0.99f
      cfe_ptr->timestep_rainfall_input_m = cfe_ptr->aorc.precip_kg_per_m2 /1000;
    else
      // Set the current rainfall input to the right place in the forcing array.
      // divide by 1000 to convert from mm/h to m w/ 1h timestep as per t-shirt_0.99f
      cfe_ptr->timestep_rainfall_input_m = cfe_ptr->forcing_data_precip_kg_per_m2[cfe_ptr->current_time_step] /1000;
    
    //Adjust the rainfall input by a potential fraction of the time step
    cfe_ptr->timestep_rainfall_input_m *= cfe_ptr->time_step_fraction;
    //Accumulate volume for mass balance
    cfe_ptr->vol_struct.volin += cfe_ptr->timestep_rainfall_input_m;
    
    run_cfe(cfe_ptr);
        
    // Advance the model time 
    cfe_ptr->current_time_step += 1;

    return BMI_SUCCESS;
}


static int Update_until (Bmi *self, double t)
{
    // https://bmi.readthedocs.io/en/latest/#update-until
    // "the time argument can be a non-integral multiple of time steps"

    cfe_state_struct* cfe_ptr = ((cfe_state_struct *) self->data);
    
    double dt;
    double now;

    if(self->get_time_step (self, &dt) == BMI_FAILURE)
        return BMI_FAILURE;

    if(self->get_current_time(self, &now) == BMI_FAILURE)
        return BMI_FAILURE;    

    {
    
    int n;
    double frac;
    const double n_steps = (t - now) / dt;
    for (n=0; n<(int)n_steps; n++) {
        Update (self);
    }
    frac = n_steps - (int)n_steps;
    if (frac > 0){
        printf("WARNING: CFE trying to update a fraction of a timestep\n");
        
        // change timestep to remaining fraction & call update()
        cfe_ptr->time_step_size = frac * dt;
        //Record the time step fraction so `Update` can adjust input if needed
        cfe_ptr->time_step_fraction = frac;
        Update (self);
        // set back to original
        cfe_ptr->time_step_fraction = 1.0;
        cfe_ptr->time_step_size = dt;
    }
    
    }

  return BMI_SUCCESS;
}


//    cfe_state_struct* cfe_ptr = ((cfe_state_struct *) self->data);
//    
//    int t_int = (int) t;
//    if ((t - ((double)t_int)) != 0)
//        return BMI_FAILURE;
//
//    for (int j = 0; j < t_int; j++){
//
//        self->update(self); 
//        if (cfe_ptr->verbosity > 1)
//            print_cfe_flux_at_timestep(cfe_ptr);
//
//    }
//    return BMI_SUCCESS;
//}


static int Finalize (Bmi *self)
{
    // Function assumes everything that is needed is retrieved from the model before Finalize is called.
    // mass_balance_check(self);
    if (self) {
        cfe_state_struct* model = (cfe_state_struct *)(self->data);
        if( model->forcing_data_precip_kg_per_m2 != NULL )
            free(model->forcing_data_precip_kg_per_m2);
        if( model->forcing_data_time != NULL )
            free(model->forcing_data_time);
        
        if( model->giuh_ordinates != NULL )
            free(model->giuh_ordinates);
        if( model->nash_storage != NULL )
            free(model->nash_storage);
        if( model->runoff_queue_m_per_timestep != NULL )
            free(model->runoff_queue_m_per_timestep);
        if( model->flux_Qout_m != NULL )
            free(model->flux_Qout_m);
        
        /* xinanjiang_dev: changing name to the more general "direct runoff"
        if( model->flux_Schaake_output_runoff_m != NULL )
            free(model->flux_Schaake_output_runoff_m);*/
        if( model->flux_output_direct_runoff_m != NULL )
            free(model->flux_output_direct_runoff_m);

        if( model->flux_from_deep_gw_to_chan_m != NULL )
            free(model->flux_from_deep_gw_to_chan_m);
        if( model->flux_giuh_runoff_m != NULL )
            free(model->flux_giuh_runoff_m);
        if( model->flux_lat_m != NULL )
            free(model->flux_lat_m);
        if( model->flux_nash_lateral_runoff_m != NULL )
            free(model->flux_nash_lateral_runoff_m);
        if( model->flux_perc_m != NULL )
            free(model->flux_perc_m);
        
        free(self->data);
    }

    return BMI_SUCCESS;
}


static int Get_adjusted_index_for_variable(const char *name)
{
    // Get an "adjusted index" value for the associated variable, where this is its index within the
    // associated names array, plus either:
    //      0 if it is in the output variable array or
    //      OUTPUT_VAR_NAME_COUNT if it is in the input variable array
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++){
        if (strcmp(name, output_var_names[i]) == 0)
            return i;
    }

    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++){
        if (strcmp(name, input_var_names[i]) == 0)
            return i + OUTPUT_VAR_NAME_COUNT;
    }

    for (i = 0; i < PARAM_VAR_NAME_COUNT; i++){
    if (strcmp(name, param_var_names[i]) == 0)
        return i + INPUT_VAR_NAME_COUNT + OUTPUT_VAR_NAME_COUNT;
    }

    return -1;
}


static int Get_var_grid(Bmi *self, const char *name, int *grid)
{
    // Check to see if in output array first
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            *grid = output_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            *grid = input_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    grid[0] = '\0';
    
    return BMI_FAILURE;
}


static int Get_var_type (Bmi *self, const char *name, char * type)
{
    //TODO may need to expose type info for "hidden" parameter values
    //or impose BMI side casting of double to whatever type
    // Check to see if in output array first
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(type, output_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(type, input_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }

    for (i = 0; i < PARAM_VAR_NAME_COUNT; i++) {
        if (strcmp(name, param_var_names[i]) == 0) {
            strncpy(type, param_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }

    // If we get here, it means the variable name wasn't recognized
    type[0] = '\0';
    return BMI_FAILURE;
}


static int Get_var_itemsize (Bmi *self, const char *name, int * size)
{
    //TODO may need to implement for "hidden" parameter variables
    //or impose BMI side casting from double to desired type
    char type[BMI_MAX_TYPE_NAME];
    int type_result = Get_var_type(self, name, type);
    if (type_result != BMI_SUCCESS) {
        return BMI_FAILURE;
    }

    if (strcmp (type, "double") == 0) {
        *size = sizeof(double);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "double") == 0) {
        *size = sizeof(double);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "int") == 0) {
        *size = sizeof(int);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "short") == 0) {
        *size = sizeof(short);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "long") == 0) {
        *size = sizeof(long);
        return BMI_SUCCESS;
    }
    else {
        *size = 0;
        return BMI_FAILURE;
    }
}


static int Get_var_location (Bmi *self, const char *name, char * location)
{
    // Check to see if in output array first
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(location, output_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(location, input_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    location[0] = '\0';
    return BMI_FAILURE;
}



static int Get_var_units (Bmi *self, const char *name, char * units)
{
    // Check to see if in output array first
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(units, output_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(units, input_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    units[0] = '\0';
    return BMI_FAILURE;
}


static int Get_var_nbytes (Bmi *self, const char *name, int * nbytes)
{
    int item_size;
    int item_size_result = Get_var_itemsize(self, name, &item_size);
    if (item_size_result != BMI_SUCCESS) {
        return BMI_FAILURE;
    }
    int item_count = -1;
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            item_count = input_var_item_count[i];
            break;
        }
    }
    if (item_count < 1) {
        for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
            if (strcmp(name, output_var_names[i]) == 0) {
                item_count = output_var_item_count[i];
                break;
            }
        }
    }
    if (item_count < 1)
        item_count = ((cfe_state_struct *) self->data)->num_timesteps;

    *nbytes = item_size * item_count;
    return BMI_SUCCESS;
}


static int Get_value_ptr (Bmi *self, const char *name, void **dest)
{
    /*********** Calibration Params Hacked ************/
    if (strcmp (name, "maxsmc") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->NWM_soil_params.smcmax;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "satdk") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->NWM_soil_params.satdk;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "slope") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->NWM_soil_params.slop;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "b") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->NWM_soil_params.bb;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "multiplier") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->NWM_soil_params.mult;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "Klf") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->soil_reservoir.coeff_secondary;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "Kn") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->K_nash;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "Cgw") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->gw_reservoir.coeff_primary;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "expon") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->gw_reservoir.exponent_primary;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "max_gw_storage") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->gw_reservoir.storage_max_m;
        return BMI_SUCCESS;
    }
    
    
    //NOT MESSING WITH nash_n (number of nash cascades) cause it is a bit of a side
    //effecting value when it changes

    /***********************************************************/
    /***********    OUTPUT   ***********************************/
    /***********************************************************/
    if (strcmp (name, "RAIN_RATE") == 0) {   //jmframe: Seems unnecessary to have rain rate as an output variable.
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->aorc.precip_kg_per_m2;
        return BMI_SUCCESS;
    }

    /* xinanjiang_dev
        changing the name to the more general "direct runoff"
    if (strcmp (name, "SCHAAKE_OUTPUT_RUNOFF") == 0) {
        *dest = (void*) ((cfe_state_struct *)(self->data))->flux_Schaake_output_runoff_m;
        return BMI_SUCCESS;
    }*/
    if (strcmp (name, "DIRECT_RUNOFF") == 0) {
        *dest = (void*) ((cfe_state_struct *)(self->data))->flux_output_direct_runoff_m;
        return BMI_SUCCESS;
    }

    if (strcmp (name, "GIUH_RUNOFF") == 0) {
        *dest = (void *) ((cfe_state_struct *)(self->data))->flux_giuh_runoff_m;
        return BMI_SUCCESS;
    }

    if (strcmp (name, "NASH_LATERAL_RUNOFF") == 0) {
        *dest = (void *) ((cfe_state_struct *)(self->data))->flux_nash_lateral_runoff_m;
        return BMI_SUCCESS;
    }

    if (strcmp (name, "DEEP_GW_TO_CHANNEL_FLUX") == 0) {
        *dest = (void *) ((cfe_state_struct *)(self->data))->flux_from_deep_gw_to_chan_m;
        return BMI_SUCCESS;
    }

    if (strcmp (name, "Q_OUT") == 0) {
        *dest = ((cfe_state_struct *)(self->data))->flux_Qout_m;
        return BMI_SUCCESS;
    }

    /***********************************************************/
    /***********    INPUT    ***********************************/
    /***********************************************************/
    if (strcmp (name, "water_potential_evaporation_flux") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr-> et_struct.potential_et_m_per_s;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "atmosphere_water__liquid_equivalent_precipitation_rate") == 0) {
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->aorc.precip_kg_per_m2;
        return BMI_SUCCESS;
    }

    return BMI_FAILURE;
}


static int Get_value_at_indices (Bmi *self, const char *name, void *dest, int *inds, int len)
{
    if (len < 1)
        return BMI_FAILURE;

    // Start by getting an "adjusted index" value for the associated variable, where this is its index within the
    // associated names array, plus either:
    //      0 if it is in the output variable array or
    //      OUTPUT_VAR_NAME_COUNT if it is in the input variable array
    int adjusted_index = Get_adjusted_index_for_variable(name);
    if (adjusted_index < 0)
        return BMI_FAILURE;

    int var_item_size;
    int status = Get_var_itemsize(self, name, &var_item_size);
    if (status == BMI_FAILURE)
        return BMI_FAILURE;

    // For now, all variables are non-array scalar values, with only 1 item of type double
    
    // Thus, there is only ever one value to return (len must be 1) and it must always be from index 0
    if (len > 1 || inds[0] != 0) 
        return BMI_FAILURE;

    void* ptr;
    status = Get_value_ptr(self, name, &ptr);
    if (status == BMI_FAILURE)
        return BMI_FAILURE;
    memcpy(dest, ptr, var_item_size * len);
    return BMI_SUCCESS;
}


static int Get_value (Bmi *self, const char *name, void *dest)
{
    // Use nested call to "by index" version

    // Here, for now at least, we know all the variables are scalar, so
    int inds[] = {0};

    // Then we can just ...
    return Get_value_at_indices(self, name, dest, inds, 1);
}


static int Set_value_at_indices (Bmi *self, const char *name, int * inds, int len, void *src)
{
    if (len < 1)
        return BMI_FAILURE;

    // Get "adjusted_index" for variable
    int adjusted_index = Get_adjusted_index_for_variable(name);
    if (adjusted_index < 0)
        return BMI_FAILURE;

    int var_item_size;
    int status = Get_var_itemsize(self, name, &var_item_size);
    if (status == BMI_FAILURE)
        return BMI_FAILURE;

    // For now, all variables are non-array scalar values, with only 1 item of type double

    // Thus, there is only ever one value to return (len must be 1) and it must always be from index 0
    if (len > 1 || inds[0] != 0)
        return BMI_FAILURE;

    void* ptr;
    status = Get_value_ptr(self, name, &ptr);
    if (status == BMI_FAILURE)
        return BMI_FAILURE;
    memcpy(ptr, src, var_item_size * len);
    /*
    * If we want to modify the number of nash cascades, we must also side effect other
    * variables.  This "should" work, based on current program structure
    */
    /*
    if (strcmp (name, "number_nash_reservoirs") ) {
    
            //Note that this only allows adjustment to the number of cascades,
            //and each one will have initial storage value of 0.0
            //Reallocate nash_storage
            if( model->nash_storage != NULL ) free(model->nash_storage);
            model->nash_storage = malloc(sizeof(double) * model->num_lateral_flow_nash_reservoirs);
            if( model->nash_storage == NULL ) return BMI_FAILURE;
            for (j = 0; j < model->num_lateral_flow_nash_reservoirs; j++)
                model->nash_storage[j] = 0.0;
        
    }
    */
    return BMI_SUCCESS;
}


static int Set_value (Bmi *self, const char *name, void *array)
{
    // Avoid using set value, call instead set_value_at_index
    // Use nested call to "by index" version

    // Here, for now at least, we know all the variables are scalar, so
    int inds[] = {0};

    // Then we can just ...
    return Set_value_at_indices(self, name, inds, 1, array);

/*  This is the sample code from read the docs
    void * dest = NULL;
    int nbytes = 0;

    if (self->get_value_ptr(self, name, &dest) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_nbytes(self, name, &nbytes) == BMI_FAILURE)
        return BMI_FAILURE;

    memcpy (dest, array, nbytes);

    return BMI_SUCCESS;
*/    
}


static int Get_component_name (Bmi *self, char * name)
{
    strncpy (name, "The CFE Model", BMI_MAX_COMPONENT_NAME);
    return BMI_SUCCESS;
}


static int Get_input_item_count (Bmi *self, int * count)
{
    *count = INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int Get_input_var_names (Bmi *self, char ** names)
{
    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        strncpy (names[i], input_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}


static int Get_output_item_count (Bmi *self, int * count)
{
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int Get_output_var_names (Bmi *self, char ** names)
{
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        strncpy (names[i], output_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}

// ***********************************************************
// *****  NEW BMI: STATE VAR GETTER & SETTER FUNCTIONS   *****
// *****  Proposed extensions to support serialization.  *****
// ***********************************************************
static int Get_state_var_count (Bmi *self, int * count)
{
    if (!self){
        return BMI_FAILURE;
    }

    *count = STATE_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}
//--------------------------------------------------------------------------
static int Get_state_var_names (Bmi *self, char ** names)
{
    //---------------------------------------------------
    // Note: This pulls information from the var_info
    // structure defined at the top, which helps to
    // prevent implementation errors.
    //---------------------------------------------------
    // This is used for model state serialization, and
    // returns a string array of all state variable
    // names, in same order as defined in state struct.
    // These names can simply be internal vs. standard
    // names because they are not used for coupling.
    //---------------------------------------------------
    if (!self){
        return BMI_FAILURE;   
    }

    int n_state_vars = STATE_VAR_NAME_COUNT;
    int MAX_NAME_LEN = 512;
    //int MAX_NAME_LEN = BMI_MAX_VAR_NAME;
    
    for (int i = 0; i < n_state_vars; i++) {
        strncpy(names[i], var_info[i].name, MAX_NAME_LEN);

        //--------------------------------  
        // Option to print all the names
        //--------------------------------
        // if (i==0) printf(" State variable names:");
        // printf(" var name[%d] = %s\n", i, names[i]);
    }
        
    return BMI_SUCCESS;
}

//--------------------------------------------------------------------------
static int Get_state_var_types (Bmi *self, char ** types)
{
    //---------------------------------------------------
    // Note: This pulls information from the var_info
    // structure defined at the top, which helps to 
    // prevent implementation errors.   
    //---------------------------------------------------
    // This is used for model state serialization, and
    // returns a string array of all state variable
    // types, in same order as defined in state struct.
    // Later, bmi.get_var_type() may be extended to
    // get more than input & output variable types.
    //---------------------------------------------------
    if (!self){
        return BMI_FAILURE;   
    }

    int n_state_vars = STATE_VAR_NAME_COUNT;   
    int MAX_NAME_LEN = 512;
    //int MAX_NAME_LEN = BMI_MAX_VAR_NAME;

    for (int i = 0; i < n_state_vars; i++) {
        strncpy(types[i], var_info[i].type, MAX_NAME_LEN);       
        //--------------------------------  
        // Option to print all the types
        //--------------------------------
        // if (i==0) printf(" State var_types:");
        // printf(" var type[%d] = %s\n", i, types[i]);
    }
        
    return BMI_SUCCESS;
}

//--------------------------------------------------------------------------
static int Get_state_var_sizes (Bmi *self, unsigned int size_list[])
{
    //---------------------------------------------------
    // Note: This pulls information from the var_info
    // structure defined at the top, which helps to 
    // prevent implementation errors.   
    //---------------------------------------------------
    // This is used for model state serialization, and
    // returns a string array of all state variable
    // sizes, in same order as defined in state struct.
    // Size is number of array elements (not bytes).
    // Just number of elements, even for n-dim arrays.
    //---------------------------------------------------
    if (!self){
        return BMI_FAILURE;   
    }

    cfe_state_struct *state;
    state = (cfe_state_struct*) self->data;  // typecast self->data
    int n_state_vars = STATE_VAR_NAME_COUNT;

    //---------------------------------------------------
    // In this file, see the functions:
    //    itwo_alloc_cfe(), dtwo_alloc_cfe(),
    //    i_alloc_cfe() and d_alloc_cfe().
    // They add 1 "for safety" & also Fortran indexing.
    //---------------------------------------------------
    unsigned int ff_len = strlen( state->forcing_file );  //##############
    unsigned int num_giuh = state->num_giuh_ordinates + 1;
    unsigned int num_lat_flow = state->num_lateral_flow_nash_reservoirs + 1;
    // unsigned int num_giuh = state->num_giuh_ordinates;
    // unsigned int num_lat_flow = state->num_lateral_flow_nash_reservoirs;

    //-------------------------------------------------
    // Overwrite the sizes that are not 1 (now known)
    //-------------------------------------------------
    var_info[55].size = ff_len;
    var_info[74].size = num_giuh;
    var_info[75].size = num_lat_flow;
    var_info[76].size = num_giuh;

    for (int i = 0; i < n_state_vars; i++) {
        size_list[i] = var_info[i].size;
    }
  
    return BMI_SUCCESS;
}

//-----------------------------------------------------------------------
static int Get_state_var_ptrs (Bmi *self, void *ptr_list[])
{
    //----------------------------------------------
    // Return array of pointers to state variables
    // in same order as defined in state struct.
    //----------------------------------------------
    if (!self){
        return BMI_FAILURE;   
    }

    cfe_state_struct *state;
    state = (cfe_state_struct*) self->data;  // typecast self->data
    //--------------------------------------------------
    // From bmi_cfe.h:
    // struct conceptual_reservoir soil_reservoir;
    // struct conceptual_reservoir gw_reservoir;
    // struct NWM_soil_parameters NWM_soil_params;
    // struct evapotranspiration_structure et_struct;
    // struct vol_tracking_struct vol_struct;
    // struct aorc_forcing_data_cfe aorc;
    //--------------------------------------------------

    ptr_list[0]  = &(state->timestep_rainfall_input_m);
    ptr_list[1]  = &(state->soil_reservoir_storage_deficit_m);    
    ptr_list[2]  = &(state->infiltration_depth_m);
    ptr_list[3]  = &(state->gw_reservoir_storage_deficit_m);
    ptr_list[4]  = &(state->timestep_h);
    //--------------------------------
    // Vars in soil reservoir struct
    //--------------------------------
    ptr_list[5]  = &(state->soil_reservoir.is_exponential );    
    ptr_list[6]  = &(state->soil_reservoir.storage_max_m );  
    ptr_list[7]  = &(state->soil_reservoir.storage_m );  
    ptr_list[8]  = &(state->soil_reservoir.coeff_primary );  
    ptr_list[9]  = &(state->soil_reservoir.exponent_primary );  
    ptr_list[10] = &(state->soil_reservoir.storage_threshold_primary_m );  
    ptr_list[11] = &(state->soil_reservoir.storage_threshold_secondary_m );  
    ptr_list[12] = &(state->soil_reservoir.coeff_secondary );      
    ptr_list[13] = &(state->soil_reservoir.exponent_secondary );
    //------------------------------
    // Vars in gw reservoir struct
    //------------------------------ 
    ptr_list[14] = &(state->gw_reservoir.is_exponential );
    ptr_list[15] = &(state->gw_reservoir.storage_max_m );
    ptr_list[16] = &(state->gw_reservoir.storage_m );
    ptr_list[17] = &(state->gw_reservoir.coeff_primary );
    ptr_list[18] = &(state->gw_reservoir.exponent_primary );
    ptr_list[19] = &(state->gw_reservoir.storage_threshold_primary_m );    
    ptr_list[20] = &(state->gw_reservoir.storage_threshold_secondary_m );
    ptr_list[21] = &(state->gw_reservoir.coeff_secondary );
    ptr_list[22] = &(state->gw_reservoir.exponent_secondary );
    //---------------------------------
    // Vars in NWM_soil_params struct
    //---------------------------------
    ptr_list[23] = &(state->NWM_soil_params.smcmax );
    ptr_list[24] = &(state->NWM_soil_params.wltsmc);
    ptr_list[25] = &(state->NWM_soil_params.satdk);
    ptr_list[26] = &(state->NWM_soil_params.satpsi); 
    ptr_list[27] = &(state->NWM_soil_params.bb);
    ptr_list[28] = &(state->NWM_soil_params.slop);
    ptr_list[29] = &(state->NWM_soil_params.D);
    ptr_list[30] = &(state->NWM_soil_params.wilting_point_m);
    //--------------------
    // Vars in et_struct
    //--------------------       
    ptr_list[31] = &(state->et_struct.potential_et_m_per_s );
    ptr_list[32] = &(state->et_struct.potential_et_m_per_timestep );
    ptr_list[33] = &(state->et_struct.actual_et_m_per_timestep );
    //------------------------------
    // Vars in vol_tracking_struct
    //------------------------------
    ptr_list[34] = &(state->vol_struct.vol_runoff );
    ptr_list[35] = &(state->vol_struct.vol_infilt );
    ptr_list[36] = &(state->vol_struct.vol_to_soil );
    ptr_list[37] = &(state->vol_struct.vol_to_gw );
    ptr_list[38] = &(state->vol_struct.vol_soil_to_gw );
    ptr_list[39] = &(state->vol_struct.vol_soil_to_lat_flow );
    ptr_list[40] = &(state->vol_struct.volstart );
    ptr_list[41] = &(state->vol_struct.volout );
    ptr_list[42] = &(state->vol_struct.volin );
    ptr_list[43] = &(state->vol_struct.vol_from_gw ); 
    ptr_list[44] = &(state->vol_struct.vol_out_giuh );
    ptr_list[45] = &(state->vol_struct.vol_in_nash );
    ptr_list[46] = &(state->vol_struct.vol_out_nash );
    ptr_list[47] = &(state->vol_struct.vol_in_gw_start );
    ptr_list[48] = &(state->vol_struct.vol_soil_start );
    //-----------------------------------------         
    // More top-level, static allocation vars
    //-----------------------------------------
    ptr_list[49] = &(state->epoch_start_time ); 
    ptr_list[50] = &(state->num_timesteps );
    ptr_list[51] = &(state->current_time_step );
    ptr_list[52] = &(state->time_step_size );
    ptr_list[53] = &(state->is_forcing_from_bmi );
    ptr_list[54] = state->forcing_file;
    // ####### ptr_list[55] = &(state->forcing_file ); 
    ptr_list[55] = &(state->direct_runoff_params_struct.surface_partitioning_scheme );
    // ptr_list[56] = &(state->Schaake_adjusted_magic_constant_by_soil_type );
    ptr_list[56] = &(state->num_lateral_flow_nash_reservoirs);
    ptr_list[57] = &(state->K_lf);
    ptr_list[58] = &(state->K_nash);
    ptr_list[59] = &(state->num_giuh_ordinates);
    //---------------------------------------
    // Vars in aorc_forcing_data_cfe struct
    //---------------------------------------
    ptr_list[60] = &(state->aorc.precip_kg_per_m2 );
    ptr_list[61] = &(state->aorc.incoming_longwave_W_per_m2 );
    ptr_list[62] = &(state->aorc.incoming_shortwave_W_per_m2 );
    ptr_list[63] = &(state->aorc.surface_pressure_Pa );
    ptr_list[64] = &(state->aorc.specific_humidity_2m_kg_per_kg );
    ptr_list[65] = &(state->aorc.air_temperature_2m_K );
    ptr_list[66] = &(state->aorc.u_wind_speed_10m_m_per_s );
    ptr_list[67] = &(state->aorc.v_wind_speed_10m_m_per_s );
    ptr_list[68] = &(state->aorc.latitude );
    ptr_list[69] = &(state->aorc.longitude );
    ptr_list[70] = &(state->aorc.time );
    //------------------------------------------         
    // More top-level, dynamic allocation vars
    //----------------------------------------------------
    // These vars ARE pointers to different-sized arrays
    // Do not append "&".
    // Last one, verbosity, is not a pointer.
    //------------------------------------------
    ptr_list[71] = state->forcing_data_precip_kg_per_m2;
    ptr_list[72] = state->forcing_data_time;
    ptr_list[73] = state->giuh_ordinates;
    ptr_list[74] = state->nash_storage;
    ptr_list[75] = state->runoff_queue_m_per_timestep;
    ptr_list[76] = state->flux_output_direct_runoff_m;
    ptr_list[77] = state->flux_giuh_runoff_m;
    ptr_list[78] = state->flux_nash_lateral_runoff_m;
    ptr_list[79] = state->flux_from_deep_gw_to_chan_m;
    ptr_list[80] = state->flux_perc_m;
    ptr_list[81] = state->flux_lat_m;
    ptr_list[82] = state->flux_Qout_m;
    ptr_list[83] = &(state->verbosity); 
    //---------------------------------------
    // direct_runoff_params_struct vars
    // xinanjiang or schaake flag [56]
    //---------------------------------------
    ptr_list[84] = &(state->direct_runoff_params_struct.Schaake_adjusted_magic_constant_by_soil_type );
    ptr_list[85] = &(state->direct_runoff_params_struct.a_Xinanjiang_inflection_point_parameter );
    ptr_list[86] = &(state->direct_runoff_params_struct.b_Xinanjiang_shape_parameter );
    ptr_list[87] = &(state->direct_runoff_params_struct.x_Xinanjiang_shape_parameter );
    //-------------------------------------------------------------                
    return BMI_SUCCESS;
}

//-----------------------------------------------------------------------
// IDEA:  A function something like this would make it possible
//        to use a variable's data type string to typecast void
//        pointer to pointer of required type.  This would help
//        prevent implementation errors.  But this can't work
//        as shown since function's return type is not fixed.
//-----------------------------------------------------------------------
//SOME_TYPE typecast_ptr( void *src, char* type )
//
//    if (strcmp(type, 'int') == 0){
//        return *(int *) src; }
//    else if (strcmp(type, 'long') == 0){
//        return *(long *) src; }
//    else if (strcmp(type, 'float') == 0){
//        return *(float *) src; }
//    else if (strcmp(type, 'double') == 0){
//        return *(double *) src; }
//    else if (strcmp(type, 'string') == 0){
//        // return *(double *) src; }
//        }
//
//-----------------------------------------------------------------------
/*static int Set_state_var (Bmi *self, void *src, int index)
{
    // NOTE: 11.05.2021: this and other serialization functions (*_state_var_*)
    // are in current working development outside of master branch

    //----------------------------------------------------
    // Set the value (or values) for a state variable
    // using its position index within the state struct.
    //----------------------------------------------------
    if (!self){
        return BMI_FAILURE;
    }

    int n_state_vars, i;
    self->get_state_var_count(self, &n_state_vars);  
    unsigned int sizes[ n_state_vars ];
    self->get_state_var_sizes(self, sizes);
    unsigned int size = sizes[ index ];

    cfe_state_struct *state;
    state = (cfe_state_struct*) self->data;  // typecast self->data

    //----------------------------------
    // For strings, this seems to work
    //----------------------------------
//        memcpy(&state->title, src, size); }
//    else if (index == 7){
//        memcpy(&state->subcat, src, size); }
    //------------------------------
    // And this seems to work also
    //------------------------------
//     else if (index == 6){
//        for (i=0; i<size; i++) {
//            state->title[i] = *( (char *)src + i); } }
//     else if (index == 7){     
//        for (i=0; i<size; i++) {
//            state->subcat[i] = *( (char *)src + i); } }

    //---------------------------------------------
    // Set value of state variable given by index
    //--------------------------------------------------
    if (index == 0){
        state->timestep_rainfall_input_m = *(double *)src; }
    else if (index == 1){
        state->soil_reservoir_storage_deficit_m = *(double *)src; }
    else if (index == 2){
        state->infiltration_depth_m = *(double *)src; }  
    else if (index == 3){
        state->gw_reservoir_storage_deficit_m = *(double *)src; }      
    else if (index == 4){
        state->timestep_h = *(double *)src; }
    //----------------------------------------------------------------
    // soil_reservoir vars
    //----------------------------------------------------------------
    else if (index == 5){
        state->soil_reservoir.is_exponential = *(int *)src; }
    else if (index == 6){
        state->soil_reservoir.storage_max_m = *(double *)src; }
    else if (index == 7){
        state->soil_reservoir.storage_m = *(double *)src; }
    else if (index == 8){
        state->soil_reservoir.coeff_primary = *(double *)src; }  
    else if (index == 9){
        state->soil_reservoir.exponent_primary = *(double *)src; }
    else if (index == 10){
        state->soil_reservoir.storage_threshold_primary_m = *(double *)src; } 
    else if (index == 11){
        state->soil_reservoir.storage_threshold_secondary_m = *(double *)src; }
    else if (index == 12){
        state->soil_reservoir.coeff_secondary = *(double *)src; } 
    else if (index == 13){
        state->soil_reservoir.exponent_secondary = *(double *)src; }
    //----------------------------------------------------------------
    // gw_reservoir vars
    //---------------------------------------------------------------- 
    else if (index == 14){
        state->gw_reservoir.is_exponential = *(int *)src; }
    else if (index == 15){
        state->gw_reservoir.storage_max_m = *(double *)src; }                               
    else if (index == 16){
        state->gw_reservoir.storage_m = *(double *)src; } 
    else if (index == 17){
        state->gw_reservoir.coeff_primary = *(double *)src; } 
    else if (index == 18){
        state->gw_reservoir.exponent_primary = *(double *)src; }                 
    else if (index == 19){
        state->gw_reservoir.storage_threshold_primary_m = *(double *)src; }
    else if (index == 20){
        state->gw_reservoir.storage_threshold_secondary_m = *(double *)src; }         
    else if (index == 21){
        state->gw_reservoir.coeff_secondary = *(double *)src; } 
    else if (index == 22){
        state->gw_reservoir.exponent_secondary = *(double *)src; }
    //---------------------------------------------------------------- 
    // NWM_soil_params vars
    //----------------------------------------------------------------  
    else if (index == 23){
        state->NWM_soil_params.smcmax = *(double *)src; }  
    else if (index == 24){
        state->NWM_soil_params.wltsmc = *(double *)src; }                  
    else if (index == 25){
        state->NWM_soil_params.satdk = *(double *)src; }  
    else if (index == 26){
        state->NWM_soil_params.satpsi= *(double *)src; }  
    else if (index == 27){
        state->NWM_soil_params.bb = *(double *)src; }                          
    else if (index == 28){
        state->NWM_soil_params.slop = *(double *)src; }  
    else if (index == 29){
        state->NWM_soil_params.D = *(double *)src; }          
    else if (index == 30){
        state->NWM_soil_params.wilting_point_m = *(double *)src; }
    //----------------------------------------------------------------
    // et_struct vars
    //----------------------------------------------------------------
    else if (index == 31){
        state->et_struct.potential_et_m_per_s = *(double *)src; }  
    else if (index == 32){
        state->et_struct.potential_et_m_per_timestep = *(double *)src; }          
    else if (index == 33){
        state->et_struct.actual_et_m_per_timestep = *(double *)src; }
    //----------------------------------------------------------------
    // vol_struct vars
    //----------------------------------------------------------------
    else if (index == 34){
        state->vol_struct.vol_runoff = *(double *)src; }  
    else if (index == 35){
        state->vol_struct.vol_infilt = *(double *)src; }          
    else if (index == 36){
        state->vol_struct.vol_to_soil = *(double *)src; }
    else if (index == 37){
        state->vol_struct.vol_to_gw = *(double *)src; }  
    else if (index == 38){
        state->vol_struct.vol_soil_to_gw = *(double *)src; }          
    else if (index == 39){
        state->vol_struct.vol_soil_to_lat_flow = *(double *)src; }
    else if (index == 40){
        state->vol_struct.volstart = *(double *)src; }  
    else if (index == 41){
        state->vol_struct.volout = *(double *)src; }          
    else if (index == 42){
        state->vol_struct.volin = *(double *)src; }
    else if (index == 43){
        state->vol_struct.vol_from_gw = *(double *)src; }  
    else if (index == 44){
        state->vol_struct.vol_out_giuh = *(double *)src; }          
    else if (index == 45){
        state->vol_struct.vol_in_nash = *(double *)src; }
    else if (index == 46){
        state->vol_struct.vol_out_nash = *(double *)src; }  
    else if (index == 47){
        state->vol_struct.vol_in_gw_start = *(double *)src; }          
    else if (index == 48){
        state->vol_struct.vol_soil_start = *(double *)src; }
    //----------------------------------------------------------------
    // More top-level vars
    //----------------------------------------------------------------
    else if (index == 49){
        state->epoch_start_time = *(long *)src; }  
    else if (index == 50){
        state->num_timesteps = *(int *)src; }    //########  LONG?      
    else if (index == 51){
        state->current_time_step = *(int *)src; }
    else if (index == 52){
        state->time_step_size = *(int *)src; }  
    else if (index == 53){
        state->is_forcing_from_bmi= *(int *)src; }          
    else if (index == 54){
        // forcing_file is a string
        memcpy(state->forcing_file, src, size); }
        // state->forcing_file = (char *)src; }    // Doesn't work
    else if (index == 55){ 
        state->direct_runoff_params_struct.surface_partitioning_scheme = *(int *)src; }  
    else if (index == 56){
        state->num_lateral_flow_nash_reservoirs = *(int *)src; }          
    else if (index == 57){
        state->K_lf = *(double *)src; }
    else if (index == 58){
        state->K_nash = *(double *)src; }  
    else if (index == 59){
        state->num_giuh_ordinates = *(int *)src; }    
    //----------------------------------------------------------------
    // aorc forcing vars
    //----------------------------------------------------------------
    else if (index == 60){
        state->aorc.precip_kg_per_m2 = *(double *)src; }  
    else if (index == 61){
        state->aorc.incoming_longwave_W_per_m2 = *(double *)src; }     
    else if (index == 62){
        state->aorc.incoming_shortwave_W_per_m2 = *(double *)src; }
    else if (index == 63){
        state->aorc.surface_pressure_Pa = *(double *)src; }  
    else if (index == 64){
        state->aorc.specific_humidity_2m_kg_per_kg= *(double *)src; }          
    else if (index == 65){
        state->aorc.air_temperature_2m_K = *(double *)src; }
    else if (index == 66){
        state->aorc.u_wind_speed_10m_m_per_s = *(double *)src; }  
    else if (index == 67){
        state->aorc.v_wind_speed_10m_m_per_s = *(double *)src; }          
    else if (index == 68){
        state->aorc.latitude = *(double *)src; }
    else if (index == 69){
        state->aorc.longitude= *(double *)src; }  
    else if (index == 70){
        state->aorc.time = *(long *)src; } 
    //----------------------------------------------------------------
    // More top-level dynamically-allocated vars
    //  (pointers to scalars or arrays)
    //----------------------------------------------------------------
    // NOTE! Typecast ptr first, then add offset,
    //       into array, then dereference the ptr
    //       CORRECT:    *( ((double *)src) + i)  ??
    //       INCORRECT:  *( (double *)(src + i))
    //       INCORRECT:  *( (double *)src + i)  ??
    //       INCORRECT:  *(double *)src + i 
    //---------------------------------------------
    // Note: state->X is a pointer to an array
    //       We don't need to change that pointer,
    //       just the values in the array.
    //---------------------------------------------     
    else if (index == 71){
        for (i=0; i<size; i++) {
            state->forcing_data_precip_kg_per_m2[i] = *( ((double *)src) + i); } }
    else if (index == 72){
        for (i=0; i<size; i++) {
            state->forcing_data_time[i] = *( ((long *)src) + i); } }
    else if (index == 73){
        for (i=0; i<size; i++) {
            state->giuh_ordinates[i] = *( ((double *)src) + i); } }
    else if (index == 74){
        for (i=0; i<size; i++) {
            state->nash_storage[i] = *( ((double *)src) + i); } }                                    
    else if (index == 75){
        for (i=0; i<size; i++) {
            state->runoff_queue_m_per_timestep[i] = *( ((double *)src) + i); } }   
    else if (index == 76){
        for (i=0; i<size; i++) {
            state->flux_output_direct_runoff_m[i] = *( ((double *)src) + i); } } 
    else if (index == 77){
        for (i=0; i<size; i++) {
            state->flux_giuh_runoff_m[i] = *( ((double *)src) + i); } } 
    else if (index == 78){
        for (i=0; i<size; i++) {
            state->flux_nash_lateral_runoff_m[i] = *( ((double *)src) + i); } }             
    else if (index == 79){
        for (i=0; i<size; i++) {
            state->flux_from_deep_gw_to_chan_m[i] = *( ((double *)src) + i); } }                                     
    else if (index == 80){
        for (i=0; i<size; i++) {
            state->flux_perc_m[i] = *( ((double *)src) + i); } } 
    else if (index == 81){
        for (i=0; i<size; i++) {
            state->flux_lat_m[i] = *( ((double *)src) + i); } } 
    else if (index == 82){
        for (i=0; i<size; i++) {
            state->flux_Qout_m[i] = *( ((double *)src) + i); } } 
    else if (index == 83){
        // verbosity is not a pointer
        state->verbosity = *(int *)src; }
    //--------------------------------------------------------------------------
    // direct_runoff_params_struc vars (includes xinanjiang AND schaake)
    //--------------------------------------------------------------------------
    else if (index == 84){ 
        state->direct_runoff_params_struct.Schaake_adjusted_magic_constant_by_soil_type = *(double *)src; }
    else if (index == 85){ 
        state->direct_runoff_params_struct.a_Xinanjiang_inflection_point_parameter = *(double *)src; }
    else if (index == 86){ 
        state->direct_runoff_params_struct.b_Xinanjiang_shape_parameter = *(double *)src; }
    else if (index == 87){ 
        state->direct_runoff_params_struct.x_Xinanjiang_shape_parameter = *(double *)src; }       
    
    return BMI_SUCCESS;
}*/

/* Grid information */
static int Get_grid_rank (Bmi *self, int grid, int * rank)
{
    if (grid == 0) {
        *rank = 1;
        return BMI_SUCCESS;
    }
    else {
        *rank = -1;
        return BMI_FAILURE;
    }
}


static int Get_grid_size(Bmi *self, int grid, int * size)
{
    if (grid == 0) {
        *size = 1;
        return BMI_SUCCESS;
    }
    else {
        *size = -1;
        return BMI_FAILURE;
    }
}


static int Get_grid_type (Bmi *self, int grid, char * type)
{
    int status = BMI_FAILURE;

    if (grid == 0) {
        strncpy(type, "scalar", BMI_MAX_TYPE_NAME);
        status = BMI_SUCCESS;
    }
    else {
        type[0] = '\0';
        status = BMI_FAILURE;
    }
    return status;
}


/* Uniform rectilinear (grid type) */
static int Get_grid_shape(Bmi *self, int grid, int *shape)
{
    return BMI_FAILURE;
}


static int Get_grid_spacing(Bmi *self, int grid, double *spacing)
{
    return BMI_FAILURE;
}


static int Get_grid_origin(Bmi *self, int grid, double *origin)
{
    return BMI_FAILURE;
}


/* Non-uniform rectilinear, curvilinear (grid type)*/
static int Get_grid_x(Bmi *self, int grid, double *x)
{
    return BMI_FAILURE;
}


static int Get_grid_y(Bmi *self, int grid, double *y)
{
    return BMI_FAILURE;
}


static int Get_grid_z(Bmi *self, int grid, double *z)
{
    return BMI_FAILURE;
}


/*Unstructured (grid type)*/
static int Get_grid_node_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}


static int Get_grid_edge_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}


static int Get_grid_face_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}


static int Get_grid_edge_nodes(Bmi *self, int grid, int *edge_nodes)
{
    return BMI_FAILURE;
}


static int Get_grid_face_edges(Bmi *self, int grid, int *face_edges)
{
    return BMI_FAILURE;
}


static int Get_grid_face_nodes(Bmi *self, int grid, int *face_nodes)
{
    return BMI_FAILURE;
}


static int Get_grid_nodes_per_face(Bmi *self, int grid, int *nodes_per_face)
{
    return BMI_FAILURE;
}


int read_file_line_counts_cfe(const char* file_name, int* line_count, int* max_line_length)
{

    *line_count = 0;
    *max_line_length = 0;
    int current_line_length = 0;
    FILE* fp = fopen(file_name, "r");
    // Ensure exists
    if (fp == NULL) {
        return -1;
    }
    int seen_non_whitespace = 0;
    char c;
    for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
        // keep track if this line has seen any char other than space or tab
        if (c != ' ' && c != '\t' && c != '\n')
            seen_non_whitespace++;
        // Update line count, reset non-whitespace count, adjust max_line_length (if needed), and reset current line count
        if (c == '\n') {
            *line_count += 1;
            seen_non_whitespace = 0;
            if (current_line_length > *max_line_length)
                *max_line_length = current_line_length;
            current_line_length = 0;
        }
        else {
            current_line_length += 1;
        }
    }
    fclose(fp);

    // If we saw some non-whitespace char on last line, assume last line didn't have its own \n, so count needs to be
    // incremented by 1.
    if (seen_non_whitespace > 0) {
        *line_count += 1;
    }

    // Before returning, increment the max line length by 1, since the \n will be on the line also.
    *max_line_length += 1;

    return 0;
}


cfe_state_struct *new_bmi_cfe(void)
{
    cfe_state_struct *data;
    data = (cfe_state_struct *) malloc(sizeof(cfe_state_struct));
    data->time_step_size = 3600;
    data->time_step_fraction = 1.0;
    data->forcing_data_precip_kg_per_m2 = NULL;
    data->forcing_data_time = NULL;
    data->giuh_ordinates = NULL;
    data->nash_storage = NULL;
    data->runoff_queue_m_per_timestep = NULL;
    data->flux_Qout_m = NULL;

    /* xinanjiang_dev
        changing the name to the more general "direct runoff"
    data->flux_Schaake_output_runoff_m = NULL;*/
    data->flux_output_direct_runoff_m = NULL;

    data->flux_from_deep_gw_to_chan_m = NULL;
    data->flux_giuh_runoff_m = NULL;
    data->flux_lat_m = NULL;
    data->flux_nash_lateral_runoff_m = NULL;
    data->flux_perc_m = NULL;

    return data;
}

Bmi* register_bmi_cfe(Bmi *model) {
    if (model) {
        model->data = (void*)new_bmi_cfe();

        model->initialize = Initialize;
        model->update = Update;
        model->update_until = Update_until;
        model->finalize = Finalize;

        model->get_component_name = Get_component_name;
        model->get_input_item_count = Get_input_item_count;
        model->get_output_item_count = Get_output_item_count;
        model->get_input_var_names = Get_input_var_names;
        model->get_output_var_names = Get_output_var_names;

        model->get_var_grid = Get_var_grid;
        model->get_var_type = Get_var_type;
        model->get_var_itemsize = Get_var_itemsize;
        model->get_var_units = Get_var_units;
        model->get_var_nbytes = Get_var_nbytes;
        model->get_var_location = Get_var_location;

        model->get_current_time = Get_current_time;
        model->get_start_time = Get_start_time;
        model->get_end_time = Get_end_time;
        model->get_time_units = Get_time_units;
        model->get_time_step = Get_time_step;

        model->get_value = Get_value;
        model->get_value_ptr = Get_value_ptr;   // TODO: needs finished implementation
        model->get_value_at_indices = Get_value_at_indices;

        model->set_value = Set_value;
        model->set_value_at_indices = Set_value_at_indices;

        model->get_grid_size = Get_grid_size;    
        model->get_grid_rank = Get_grid_rank;    
        model->get_grid_type = Get_grid_type;    

        model->get_grid_shape = Get_grid_shape;    // N/a for grid type scalar
        model->get_grid_spacing = Get_grid_spacing;    // N/a for grid type scalar
        model->get_grid_origin = Get_grid_origin;    // N/a for grid type scalar

        model->get_grid_x = Get_grid_x;    // N/a for grid type scalar
        model->get_grid_y = Get_grid_y;    // N/a for grid type scalar
        model->get_grid_z = Get_grid_z;    // N/a for grid type scalar

        model->get_grid_node_count = Get_grid_node_count;    // N/a for grid type scalar
        model->get_grid_edge_count = Get_grid_edge_count;    // N/a for grid type scalar
        model->get_grid_face_count = Get_grid_face_count;    // N/a for grid type scalar
        model->get_grid_edge_nodes = Get_grid_edge_nodes;    // N/a for grid type scalar
        model->get_grid_face_edges = Get_grid_face_edges;    // N/a for grid type scalar
        model->get_grid_face_nodes = Get_grid_face_nodes;    // N/a for grid type scalar
        model->get_grid_nodes_per_face = Get_grid_nodes_per_face;    // N/a for grid type scalar

    }

    return model;
}

extern void run_cfe(cfe_state_struct* cfe_ptr){
    cfe(
        &cfe_ptr->soil_reservoir_storage_deficit_m,               // Set in cfe function
        cfe_ptr->NWM_soil_params,     // Set by config file
        &cfe_ptr->soil_reservoir,          // Set in "init_soil_reservoir" function 
        cfe_ptr->timestep_h,                                     // Set in initialize

    /* xinanjiang_dev
        changing the name to the more general "direct runoff"
        cfe_ptr->Schaake_adjusted_magic_constant_by_soil_type,   // Set by config file*/
        cfe_ptr->direct_runoff_params_struct,   // Set by config file, includes parameters for Schaake and/or XinanJiang*/

        cfe_ptr->timestep_rainfall_input_m,                      // Set by bmi (set value) or read from file.

     /* xinanjiang_dev
        cfe_ptr->flux_Schaake_output_runoff_m,                  // Set by cfe function*/
        cfe_ptr->flux_output_direct_runoff_m,

        &cfe_ptr->infiltration_depth_m,                          // Set by Schaake partitioning scheme
        cfe_ptr->flux_perc_m,                                   // Set to zero in definition.
        cfe_ptr->flux_lat_m,                                    // Set by CFE function after soil_resevroir calc
        &cfe_ptr->gw_reservoir_storage_deficit_m,                // Set by CFE function after soil_resevroir calc
        &cfe_ptr->gw_reservoir,      // Set in initialize and from config file
        cfe_ptr->flux_from_deep_gw_to_chan_m,                   // Set by CFE function after gw_reservoir calc
        cfe_ptr->flux_giuh_runoff_m,                            // Set in CFE by convolution_integral
        cfe_ptr->num_giuh_ordinates,                             // Set by config file with func. count_delimited_values
        cfe_ptr->giuh_ordinates,                            // Set by configuration file.
        cfe_ptr->runoff_queue_m_per_timestep,               // Set in initialize
        cfe_ptr->flux_nash_lateral_runoff_m,                    // Set in CFE from nash_cascade function
        cfe_ptr->num_lateral_flow_nash_reservoirs,               // Set from config file
        cfe_ptr->K_nash,                                         // Set from config file
        cfe_ptr->nash_storage,                              // Set from config file
        &cfe_ptr->et_struct,                                    // Set to zero with initalize. Set by BMI (set_value) during run
        cfe_ptr->flux_Qout_m,                                    // Set by CFE function
        &cfe_ptr->vol_struct,
        cfe_ptr->time_step_size
    );
}

// Functions for setting up CFE data, i.e., initializing...
/*extern void init_soil_reservoir(cfe_state_struct* cfe_ptr, double alpha_fc, double max_storage, double storage,
                                int is_storage_ratios)
{*/
extern void init_soil_reservoir(cfe_state_struct* cfe_ptr, double alpha_fc, double storage)
{
    // calculate the activation storage for the secondary lateral flow outlet in the soil nonlinear reservoir.
    // following the method in the NWM/t-shirt parameter equivalence document, assuming field capacity soil
    // suction pressure = 1/3 atm= field_capacity_atm_press_fraction * atm_press_Pa.

    // solve the integral given by Eqn. 5 in the parameter equivalence document.
    // this equation calculates the amount of water stored in the 2 m thick soil column when the water content
    // at the center of the bottom discretization (trigger_z_m, below 0.5) is at field capacity
    // Initial parentheses calc equation 3 from param equiv. doc
#define STANDARD_ATMOSPHERIC_PRESSURE_PASCALS 101325
    // This may need to be changed as follows later, but for now, use the constant value
    //double Omega = (alpha_fc * cfe->forcing_data_surface_pressure_Pa[0] / WATER_SPECIFIC_WEIGHT) - 0.5;
    double Omega = (alpha_fc * STANDARD_ATMOSPHERIC_PRESSURE_PASCALS / WATER_SPECIFIC_WEIGHT) - 0.5;
    double lower_lim = pow(Omega, (1.0 - 1.0 / cfe_ptr->NWM_soil_params.bb)) / 
                                  (1.0 - 1.0 / cfe_ptr->NWM_soil_params.bb);
    double upper_lim = pow(Omega + cfe_ptr->NWM_soil_params.D, 
                                       (1.0 - 1.0 / cfe_ptr->NWM_soil_params.bb)) /
                       (1.0 - 1.0 / cfe_ptr->NWM_soil_params.bb);

    // JMFRAME adding this, not sure if is correct. Ask FRED OGDEN
    cfe_ptr->NWM_soil_params.wilting_point_m = cfe_ptr->NWM_soil_params.wltsmc * cfe_ptr->NWM_soil_params.D;

    // soil conceptual reservoir first, two outlets, two thresholds, linear (exponent=1.0).
    cfe_ptr->soil_reservoir.is_exponential = FALSE;  // set this TRUE to use the exponential form of the discharge equation
    // this should NEVER be set to true in the soil reservoir.
    cfe_ptr->soil_reservoir.storage_max_m = cfe_ptr->NWM_soil_params.smcmax * cfe_ptr->NWM_soil_params.D;
    //  vertical percolation parameters ------------------------------------------------
    // Units of primary coefficient are m per time step
    cfe_ptr->soil_reservoir.coeff_primary = cfe_ptr->NWM_soil_params.satdk * cfe_ptr->NWM_soil_params.slop * cfe_ptr->time_step_size;
    // 1.0=linear
    //cfe_ptr->soil_reservoir.exponent_primary = 1.0;
    // i.e., field_capacity_storage_threshold_m
    cfe_ptr->soil_reservoir.storage_threshold_primary_m =
            cfe_ptr->NWM_soil_params.smcmax * pow(1.0 / cfe_ptr->NWM_soil_params.satpsi, (-1.0 / cfe_ptr->NWM_soil_params.bb)) *
            (upper_lim - lower_lim);
    // lateral flow parameters --------------------------------------------------------
    // 0.0 to deactiv. else =lateral_flow_linear_reservoir_constant;   // m per ts
    // TODO: look at whether K_lf needs to be a derived (i.e., via get_K_lf_for_time_step()) or explicit parameter
    cfe_ptr->soil_reservoir.coeff_secondary = cfe_ptr->K_lf;
    // 1.0=linear
    //cfe_ptr->soil_reservoir.exponent_secondary = 1.0;
    // making them the same, but they don't have 2B
    cfe_ptr->soil_reservoir.storage_threshold_secondary_m = cfe_ptr->soil_reservoir.storage_threshold_primary_m;

    // Negative amounts are always ignored and just considered emtpy
    if (storage < 0.0) storage = 0.0;
    cfe_ptr->soil_reservoir.storage_m = storage;

    //cfe_ptr->soil_reservoir.storage_m = init_reservoir_storage(is_storage_ratios, storage, max_storage);
}

/*extern double init_reservoir_storage(int is_ratio, double amount, double max_amount) {
    // Negative amounts are always ignored and just considered emtpy
    if (amount < 0.0) {
        return 0.0;
    }
    // When not a ratio (and positive), just return the literal amount
    if (is_ratio == FALSE) {
        return amount;
    }
    // When between 0 and 1, return the simple ratio computation
    if (amount <= 1.0) {
        return max_amount * amount;
    }
    // Otherwise, just return the literal amount, and assume the is_ratio value was invalid
    // TODO: is this the best way to handle this?
    else {
        return amount;
    }
}*/

extern void initialize_volume_trackers(cfe_state_struct* cfe_ptr){
    cfe_ptr->vol_struct.volin = 0;
    cfe_ptr->vol_struct.vol_runoff = 0;
    cfe_ptr->vol_struct.vol_infilt = 0;
    cfe_ptr->vol_struct.vol_to_soil = 0;
    cfe_ptr->vol_struct.vol_to_gw = 0;
    cfe_ptr->vol_struct.vol_soil_to_gw = 0;
    cfe_ptr->vol_struct.vol_soil_to_lat_flow = 0;
    cfe_ptr->vol_struct.volout = 0;
    cfe_ptr->vol_struct.vol_from_gw = 0;
    cfe_ptr->vol_struct.vol_out_giuh = 0;
    cfe_ptr->vol_struct.vol_in_nash = 0;
    cfe_ptr->vol_struct.vol_out_nash = 0;
    cfe_ptr->vol_struct.volstart       += cfe_ptr->gw_reservoir.storage_m;    // initial mass balance checks in g.w. reservoir
    cfe_ptr->vol_struct.vol_in_gw_start = cfe_ptr->gw_reservoir.storage_m;  
    cfe_ptr->vol_struct.volstart          += cfe_ptr->soil_reservoir.storage_m;    // initial mass balance checks in soil reservoir
    cfe_ptr->vol_struct.vol_soil_start     = cfe_ptr->soil_reservoir.storage_m;
    cfe_ptr->vol_struct.vol_et_from_soil = 0;
    cfe_ptr->vol_struct.vol_et_from_rain = 0;
    cfe_ptr->vol_struct.vol_et_to_atm = 0;
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS VERBOSITY AND TROUBLESHOOTING */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
extern void print_cfe_flux_header(){
    printf("#    ,            hourly ,  direct,   giuh ,lateral,  base,   total\n");
    printf("#Time,           rainfall,  runoff,  runoff, flow  ,  flow,  discharge\n");
    printf("# (h),             (mm)   ,  (mm) ,   (mm) , (mm)  ,  (mm),    (mm)\n");
}
extern void print_cfe_flux_at_timestep(cfe_state_struct* cfe_ptr){
    printf("%d %lf %lf %lf %lf %lf %lf\n",
                           cfe_ptr->current_time_step,
                           cfe_ptr->timestep_rainfall_input_m*1000.0,
                           
                           /* xinanjiang_dev
                           *cfe_ptr->flux_Schaake_output_runoff_m*1000.0,*/
                           *cfe_ptr->flux_output_direct_runoff_m*1000.0,

                           *cfe_ptr->flux_giuh_runoff_m*1000.0,
                           *cfe_ptr->flux_nash_lateral_runoff_m*1000.0, 
                           *cfe_ptr->flux_from_deep_gw_to_chan_m*1000.0,
                           *cfe_ptr->flux_Qout_m*1000.0 );
}

extern void mass_balance_check(cfe_state_struct* cfe_ptr){
    //-----------------------------------------------------------
    // PERFORM MASS BALANCE CHECKS AND WRITE RESULTS TO stderr.
    //-----------------------------------------------------------
    
    double volend= cfe_ptr->soil_reservoir.storage_m+cfe_ptr->gw_reservoir.storage_m;
    double vol_in_gw_end = cfe_ptr->gw_reservoir.storage_m;
    double vol_end_giuh = 0.0;
    double vol_in_nash_end = 0.0;
    double vol_soil_end;
    
    // the GIUH queue might have water in it at the end of the simulation, so sum it up.
    for(i=0;i<cfe_ptr->num_giuh_ordinates;i++) vol_end_giuh+=cfe_ptr->runoff_queue_m_per_timestep[i];
    
    for(i=0;i<cfe_ptr->num_lateral_flow_nash_reservoirs;i++)  vol_in_nash_end+=cfe_ptr->nash_storage[i];
    
    vol_soil_end=cfe_ptr->soil_reservoir.storage_m;
    
    double global_residual;

    /* xinanjiang_dev
    double schaake_residual;*/
    double direct_residual;

    double giuh_residual;
    double soil_residual;
    double nash_residual;
    double gw_residual;
    
    global_residual = cfe_ptr->vol_struct.volstart + cfe_ptr->vol_struct.volin - 
                      cfe_ptr->vol_struct.volout - volend - vol_end_giuh;
    printf("GLOBAL MASS BALANCE\n");
    printf("  initial volume: %8.4lf m\n",cfe_ptr->vol_struct.volstart);
    printf("    volume input: %8.4lf m\n",cfe_ptr->vol_struct.volin);
    printf("   volume output: %8.4lf m\n",cfe_ptr->vol_struct.volout);
    printf("    final volume: %8.4lf m\n",volend);
    printf("        residual: %6.4e m\n",global_residual);
    if(cfe_ptr->vol_struct.volin>0.0) printf("global pct. err: %6.4e percent of inputs\n",global_residual/cfe_ptr->vol_struct.volin*100.0);
    else          printf("global pct. err: %6.4e percent of initial\n",global_residual/cfe_ptr->vol_struct.volstart*100.0);
    if(!is_fabs_less_than_epsilon(global_residual,1.0e-12)) 
                  printf("WARNING: GLOBAL MASS BALANCE CHECK FAILED\n");
    
    /* xinanjiang_dev
    schaake_residual = cfe_ptr->vol_struct.volin - cfe_ptr->vol_struct.vol_sch_runoff - cfe_ptr->vol_struct.vol_sch_infilt;
    printf(" SCHAAKE MASS BALANCE\n");
    printf("  surface runoff: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_runoff);
    printf("    infiltration: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_infilt);
    printf("schaake residual: %6.4e m\n",schaake_residual);  // should equal 0.0
    if(!is_fabs_less_than_epsilon(schaake_residual,1.0e-12))
                  printf("WARNING: SCHAAKE PARTITIONING MASS BALANCE CHECK FAILED\n");*/
    direct_residual = cfe_ptr->vol_struct.volin - cfe_ptr->vol_struct.vol_runoff - cfe_ptr->vol_struct.vol_infilt-cfe_ptr->vol_struct.vol_et_from_rain;
    printf(" DIRECT RUNOFF MASS BALANCE\n");
    printf("  surface runoff: %8.4lf m\n",cfe_ptr->vol_struct.vol_runoff);
    printf("    infiltration: %8.4lf m\n",cfe_ptr->vol_struct.vol_infilt);
    printf("    vol_et_from_rain: %8.4lf m\n",cfe_ptr->vol_struct.vol_et_from_rain);
    printf("direct residual: %6.4e m\n",direct_residual);  // should equal 0.0
    if(!is_fabs_less_than_epsilon(direct_residual,1.0e-12))
                  printf("WARNING: DIRECT RUNOFF PARTITIONING MASS BALANCE CHECK FAILED\n");
    
    /* xinanjiang_dev
    giuh_residual = cfe_ptr->vol_struct.vol_out_giuh - cfe_ptr->vol_struct.vol_sch_runoff - vol_end_giuh;   */
    giuh_residual = cfe_ptr->vol_struct.vol_runoff - cfe_ptr->vol_struct.vol_out_giuh - vol_end_giuh;
    printf(" GIUH MASS BALANCE\n");

    /* xinanjiang_dev
    printf("  vol. into giuh: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_runoff);    */
    printf("  vol. into giuh: %8.4lf m\n",cfe_ptr->vol_struct.vol_runoff);
    fprintf(stderr,"   vol. out giuh: %8.4lf m\n",cfe_ptr->vol_struct.vol_out_giuh);
    fprintf(stderr," vol. end giuh q: %8.4lf m\n",cfe_ptr->vol_struct.vol_end_giuh);
    printf("   giuh residual: %6.4e m\n",giuh_residual);  // should equal zero
    if(!is_fabs_less_than_epsilon(giuh_residual,1.0e-12))
                  printf("WARNING: GIUH MASS BALANCE CHECK FAILED\n");

    /* xinanjiang_dev 
    soil_residual=cfe_ptr->vol_struct.vol_soil_start + cfe_ptr->vol_struct.vol_sch_infilt -      */
    soil_residual=cfe_ptr->vol_struct.vol_soil_start + cfe_ptr->vol_struct.vol_infilt -
                  cfe_ptr->vol_struct.vol_soil_to_lat_flow - vol_soil_end - cfe_ptr->vol_struct.vol_to_gw - cfe_ptr->vol_struct.vol_et_from_soil;
                  
    printf(" SOIL WATER CONCEPTUAL RESERVOIR MASS BALANCE\n");
    printf("   init soil vol: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_start);     

    /* xinanjiang_dev
    printf("  vol. into soil: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_infilt);    */
    printf("  vol. into soil: %8.4lf m\n",cfe_ptr->vol_struct.vol_infilt);
    printf("vol.soil2latflow: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_to_lat_flow);
    printf(" vol. soil to gw: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_to_gw);
    printf(" final vol. soil: %8.4lf m\n",vol_soil_end);  
    printf(" vol. et from soil: %8.4lf m\n",cfe_ptr->vol_struct.vol_et_from_soil);  
    printf("vol. soil resid.: %6.4e m\n",soil_residual);
    if(!is_fabs_less_than_epsilon(soil_residual,1.0e-12))
                   printf("WARNING: SOIL CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");
    
    nash_residual=cfe_ptr->vol_struct.vol_in_nash - cfe_ptr->vol_struct.vol_out_nash - vol_in_nash_end;
    printf(" NASH CASCADE CONCEPTUAL RESERVOIR MASS BALANCE\n");
    printf("    vol. to nash: %8.4lf m\n",cfe_ptr->vol_struct.vol_in_nash);
    printf("  vol. from nash: %8.4lf m\n",cfe_ptr->vol_struct.vol_out_nash);
    printf(" final vol. nash: %8.4lf m\n",vol_in_nash_end);
    printf("nash casc resid.: %6.4e m\n",nash_residual);
    if(!is_fabs_less_than_epsilon(nash_residual,1.0e-12))
                   printf("WARNING: NASH CASCADE CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");
    
    
    gw_residual = cfe_ptr->vol_struct.vol_in_gw_start + cfe_ptr->vol_struct.vol_to_gw - cfe_ptr->vol_struct.vol_from_gw - vol_in_gw_end;
    printf(" GROUNDWATER CONCEPTUAL RESERVOIR MASS BALANCE\n");
    printf("init gw. storage: %8.4lf m\n",cfe_ptr->vol_struct.vol_in_gw_start);
    printf("       vol to gw: %8.4lf m\n",cfe_ptr->vol_struct.vol_to_gw);
    printf("     vol from gw: %8.4lf m\n",cfe_ptr->vol_struct.vol_from_gw);
    printf("final gw.storage: %8.4lf m\n",vol_in_gw_end);
    printf("    gw. residual: %6.4e m\n",gw_residual);
    if(!is_fabs_less_than_epsilon(gw_residual,1.0e-12))
                   fprintf(stderr,"WARNING: GROUNDWATER CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS JUST UTILITY MEMORY AND TIME FUNCTION CODE */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


/*####################################################################*/
/*########################### PARSE LINE #############################*/
/*####################################################################*/
void parse_aorc_line_cfe(char *theString,long *year,long *month, long *day,long *hour,long *minute, double *second,
                struct aorc_forcing_data_cfe *aorc)
{
    int start,end;
    int wordlen;
    char theWord[150];
    
    //len=strlen(theString);
    
    start=0; /* begin at the beginning of theString */
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *year=atol(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *month=atol(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *day=atol(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *hour=atol(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *minute=atol(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    *second=(double)atof(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->precip_kg_per_m2=atof(theWord);
    //printf("%s, %s, %lf, %lf\n", theString, theWord, (double)atof(theWord), aorc->precip_kg_per_m2);
                  
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->incoming_longwave_W_per_m2=(double)atof(theWord);   
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->incoming_shortwave_W_per_m2=(double)atof(theWord);   
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->surface_pressure_Pa=(double)atof(theWord);           
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->specific_humidity_2m_kg_per_kg=(double)atof(theWord);
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->air_temperature_2m_K=(double)atof(theWord);          
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->u_wind_speed_10m_m_per_s=(double)atof(theWord);      
    
    get_word_cfe(theString,&start,&end,theWord,&wordlen);
    aorc->v_wind_speed_10m_m_per_s=(double)atof(theWord);      

    //----------------------------------------------    
    // Is this needed?  It has not been set. (SDP)
    //----------------------------------------------
    aorc->time = -9999.0;
    //######################################################
     
    return;
    }

/*####################################################################*/
/*############################## GET WORD ############################*/
/*####################################################################*/
void get_word_cfe(char *theString, int *start, int *end, char *theWord, int *wordlen) {
    int i, lenny, j;
    lenny = strlen(theString);

    while (theString[*start] == ' ' || theString[*start] == '\t') {
        (*start)++;
    };

    j = 0;
    for (i = (*start); i < lenny; i++) {
        if (theString[i] != ' ' && theString[i] != '\t' && theString[i] != ',' && theString[i] != ':' &&
            theString[i] != '/') {
            theWord[j] = theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]==',') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]=='e') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else {
            break;
        }
    }
    theWord[j] = '\0';
    *start = i + 1;
    *wordlen = strlen(theWord);
    return;
}

/****************************************/
void itwo_alloc_cfe(int ***array, int rows, int cols) {
    int i, frows, fcols;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (int **) malloc(frows * sizeof(int *));
    if (*array) {
        memset((*array), 0, frows * sizeof(int *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (int *) malloc(fcols * sizeof(int));
            if ((*array)[i] == NULL) {
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(int));
        }
    }
    return;
}


void dtwo_alloc_cfe(double ***array, int rows, int cols) {
    int i, frows, fcols;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (double **) malloc(frows * sizeof(double *));
    if (*array) {
        memset((*array), 0, frows * sizeof(double *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (double *) malloc(fcols * sizeof(double));
            if ((*array)[i] == NULL) {
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(double));
        }
    }
    return;
}


void d_alloc_cfe(double **var, int size) {
    size++;  /* just for safety */

    *var = (double *) malloc(size * sizeof(double));
    if (*var == NULL) {
        printf("Problem allocating memory for array in d_alloc.\n");
        return;
    } else
        memset(*var, 0, size * sizeof(double));
    return;
}

void i_alloc_cfe(int **var, int size) {
    size++;  /* just for safety */

    *var = (int *) malloc(size * sizeof(int));
    if (*var == NULL) {
        printf("Problem allocating memory in i_alloc\n");
        return;
    } else
        memset(*var, 0, size * sizeof(int));
    return;
}

/*
 * convert Gregorian days to Julian date
 *
 * Modify as needed for your application.
 *
 * The Julian day starts at noon of the Gregorian day and extends
 * to noon the next Gregorian day.
 *
 */
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
*/


double greg_2_jul_cfe(
        long year,
        long mon,
        long day,
        long h,
        long mi,
        double se) {
    long m = mon, d = day, y = year;
    long c, ya, j;
    double seconds = h * 3600.0 + mi * 60 + se;

    if (m > 2)
        m -= 3;
    else {
        m += 9;
        --y;
    }
    c = y / 100L;
    ya = y - (100L * c);
    j = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
    if (seconds < 12 * 3600.0) {
        j--;
        seconds += 12.0 * 3600.0;
    } else {
        seconds = seconds - 12.0 * 3600.0;
    }
    return (j + (seconds / 3600.0) / 24.0);
}

/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
**
** Modified by FLO 4/99 to account for nagging round off error 
**
*/
void calc_date_cfe(double jd, long *y, long *m, long *d, long *h, long *mi,
               double *sec) {

    long j;
    double tmp;
    double frac;

    j = (long) jd;
    frac = jd - j;

    if (frac >= 0.5) {
        frac = frac - 0.5;
        j++;
    } else {
        frac = frac + 0.5;
    }

    j -= 1721119L;
    *y = (4L * j - 1L) / 146097L;
    j = 4L * j - 1L - 146097L * *y;
    *d = j / 4L;
    j = (4L * *d + 3L) / 1461L;
    *d = 4L * *d + 3L - 1461L * j;
    *d = (*d + 4L) / 4L;
    *m = (5L * *d - 3L) / 153L;
    *d = 5L * *d - 3 - 153L * *m;
    *d = (*d + 5L) / 5L;
    *y = 100L * *y + j;
    if (*m < 10)
        *m += 3;
    else {
        *m -= 9;
        *y = *y + 1; /* Invalid use: *y++. Modified by Tony */
    }

    /* if (*m < 3) *y++; */
    /* incorrectly repeated the above if-else statement. Deleted by Tony.*/

    tmp = 3600.0 * (frac * 24.0);
    *h = (long) (tmp / 3600.0);
    tmp = tmp - *h * 3600.0;
    *mi = (long) (tmp / 60.0);
    *sec = tmp - *mi * 60.0;
}

int dayofweek(double j) {
    j += 0.5;
    return (int) (j + 1) % 7;
}

