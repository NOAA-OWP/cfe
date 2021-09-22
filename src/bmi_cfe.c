#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"
#include <time.h>
#ifndef WATER_SPECIFIC_WEIGHT
#define WATER_SPECIFIC_WEIGHT 9810
#endif

#define CFE_DEGUG 0

#define INPUT_VAR_NAME_COUNT 2
#define OUTPUT_VAR_NAME_COUNT 6

int i = 0;
int j = 0;

// Don't forget to update Get_value/Get_value_at_indices (and setter) implementation if these are adjusted
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
        "RAIN_RATE",
        "SCHAAKE_OUTPUT_RUNOFF",
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
        "kg m-2", //"atmosphere_water__liquid_equivalent_precipitation_rate"
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

int read_init_config_cfe(const char* config_file, cfe_state_struct* model, double* alpha_fc, double* soil_storage,
                     int* is_soil_storage_ratio)
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
    int is_soil_params__mult_set = FALSE;
    int is_soil_params__satdk_set = FALSE;
    int is_soil_params__satpsi_set = FALSE;
    int is_soil_params__slop_set = FALSE;
    int is_soil_params__smcmax_set = FALSE;
    int is_soil_params__wltsmc_set = FALSE;
    int is_Cgw_set = FALSE;
    int is_expon_set = FALSE;
    int is_alpha_fc_set = FALSE;
    int is_soil_storage_set = FALSE;
    int is_K_nash_set = FALSE;
    int is_K_lf_set = FALSE;
    int is_num_timesteps_set = FALSE;
    int is_verbosity_set = FALSE;

    // Keep track these in particular, because the "true" storage value may be a ratio and need both storage and max
    int is_gw_max_set = FALSE;
    int is_gw_storage_set = FALSE;

    int is_giuh_originates_string_val_set = FALSE;

    // Default value
    double refkdt = 3.0;

    int is_gw_storage_ratio = FALSE;
    double gw_storage_literal;
    // Also keep track of Nash stuff and properly set at the end of reading the config file
    int num_nash_lf = 2;
    char* nash_storage_string_val;
    int is_nash_storage_string_val_set = FALSE;
    // Similarly as for Nash, track stuff for GIUH ordinates
    int num_giuh_ordinates = 1;
    char* giuh_originates_string_val;


    // Additionally,

    for (i = 0; i < config_line_count; i++) {
        char *param_key, *param_value;
        fgets(config_line, max_config_line_length + 1, fp);
#if CFE_DEGUG >= 3
        printf("Line value: ['%s']\n", config_line);
#endif
        char* config_line_ptr = config_line;
        config_line_ptr = strsep(&config_line_ptr, "\n");
        param_key = strsep(&config_line_ptr, "=");
        param_value = strsep(&config_line_ptr, "=");

#if CFE_DEGUG >= 2
        printf("Config Value - Param: '%s' | Value: '%s'\n", param_key, param_value);
#endif

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
            continue;
        }
        if (strcmp(param_key, "soil_params.bb") == 0 || strcmp(param_key, "soil_params.b") == 0) {
            model->NWM_soil_params.bb = strtod(param_value, NULL);
            is_soil_params__bb_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.multiplier") == 0 || strcmp(param_key, "soil_params.mult") == 0) {
            model->NWM_soil_params.mult = strtod(param_value, NULL);
            is_soil_params__mult_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.satdk") == 0) {
            model->NWM_soil_params.satdk = strtod(param_value, NULL);
            is_soil_params__satdk_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_params.satpsi") == 0) {
            model->NWM_soil_params.satpsi = strtod(param_value, NULL);
            is_soil_params__satpsi_set = TRUE;
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
        if (strcmp(param_key, "max_gw_storage") == 0) {
            model->gw_reservoir.storage_max_m = strtod(param_value, NULL);
            is_gw_max_set = TRUE;
            // Also set the true storage if storage was already read and was a ratio, and so we were waiting for this
            if (is_gw_storage_set == TRUE && is_gw_storage_ratio == TRUE) {
                model->gw_reservoir.storage_m = (gw_storage_literal / 100.0) * model->gw_reservoir.storage_max_m;
            }
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
            is_gw_storage_set = TRUE;
            char* trailing_chars;
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
            }
            continue;
        }
        if (strcmp(param_key, "alpha_fc") == 0) {
            *alpha_fc = strtod(param_value, NULL);
            is_alpha_fc_set = TRUE;
            continue;
        }
        if (strcmp(param_key, "soil_storage") == 0) {
            char* trailing_chars;
            double parsed_value = strtod(param_value, &trailing_chars);
            *is_soil_storage_ratio = strcmp(trailing_chars, "%") == 0 ? TRUE : FALSE;
            *soil_storage = *is_soil_storage_ratio == TRUE ? (parsed_value / 100.0) : parsed_value;
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
    if (is_soil_params__mult_set == FALSE) {
#if CFE_DEGUG >= 1
        printf("Config param 'soil_params.mult' not found in config file\n");
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
    if (is_num_timesteps_set == FALSE && model->forcing_file == "BMI") {
#if CFE_DEGUG >= 1
        printf("Config param 'num_timesteps' not found in config file\n");
#endif
        return BMI_FAILURE;
    }

#if CFE_DEGUG >= 1
    printf("All CFE config params present\n");
#endif

    model->Schaake_adjusted_magic_constant_by_soil_type = refkdt * model->NWM_soil_params.satdk / 0.000002;

#if CFE_DEGUG >= 1
    printf("Schaake Magic Constant calculated\n");
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

    cfe_state_struct* cfe_bmi_data_ptr;

    if (!self || !file)
        return BMI_FAILURE;
    else
        cfe_bmi_data_ptr = (cfe_state_struct *) self->data;

    cfe_bmi_data_ptr->current_time_step = 0;

    double alpha_fc, max_soil_storage, S_soil;
    int is_S_soil_ratio;

    int config_read_result = read_init_config_cfe(file, cfe_bmi_data_ptr, &alpha_fc, &S_soil, &is_S_soil_ratio);
    if (config_read_result == BMI_FAILURE)
        return BMI_FAILURE;

    // time_step_size is set to 3600 in the "new_bmi_cfe" function.
    cfe_bmi_data_ptr->timestep_h = cfe_bmi_data_ptr->time_step_size / 3600.0;
     
    max_soil_storage = cfe_bmi_data_ptr->NWM_soil_params.D * cfe_bmi_data_ptr->NWM_soil_params.smcmax;


    /***********************************************************************
       JMFRAME: Moved these up before the read forcing line,
                Since we need them even if we don't read forcings from file.
    ************************************************************************/
//    cfe_bmi_data_ptr->flux_overland_m = malloc(sizeof(double));    //NOT NEEDED, redundant with flux_Schaake_output_runoff_m
    cfe_bmi_data_ptr->flux_Schaake_output_runoff_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_Qout_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_from_deep_gw_to_chan_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_giuh_runoff_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_lat_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_nash_lateral_runoff_m = malloc(sizeof(double));
    cfe_bmi_data_ptr->flux_perc_m = malloc(sizeof(double));

    /*******************************************************
       JMFRAME: Check to see where forcings come from
                The options as of May 2021 are:
                1. Read your own forcings from a file
                2. Get the forcings passed in through BMI
    *******************************************************/
    if (strcmp(cfe_bmi_data_ptr->forcing_file, "BMI") == 0){
        cfe_bmi_data_ptr->is_forcing_from_bmi = 1;
        cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2 = malloc(sizeof(double));
        cfe_bmi_data_ptr->forcing_data_time = malloc(sizeof(long));
    }
    else
    {
        cfe_bmi_data_ptr->is_forcing_from_bmi = 0;
    
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
            cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2[i] = forcings.precip_kg_per_m2 * ((double)cfe_bmi_data_ptr->time_step_size);
            cfe_bmi_data_ptr->forcing_data_time[i] = forcings.time;

            // TODO: make sure some kind of conversion isn't needed for the rain rate data
            // assumed 1000 kg/m3 density of water.  This result is mm/h;
            //rain_rate[i] = (double) aorc_data.precip_kg_per_m2;
        }
    
        cfe_bmi_data_ptr->epoch_start_time = cfe_bmi_data_ptr->forcing_data_time[0];
    } // end if is_forcing_from_bmi
    
    cfe_bmi_data_ptr->timestep_rainfall_input_m = cfe_bmi_data_ptr->forcing_data_precip_kg_per_m2[0];

    // Initialize the rest of the groundwater conceptual reservoir (some was done when reading in the config)
    cfe_bmi_data_ptr->gw_reservoir.is_exponential = TRUE;
    cfe_bmi_data_ptr->gw_reservoir.storage_threshold_primary_m = 0.0;    // 0.0 means no threshold applied
    cfe_bmi_data_ptr->gw_reservoir.storage_threshold_secondary_m = 0.0;  // 0.0 means no threshold applied
    cfe_bmi_data_ptr->gw_reservoir.coeff_secondary = 0.0;                // 0.0 means that secondary outlet is not applied
    cfe_bmi_data_ptr->gw_reservoir.exponent_secondary = 1.0;             // linear

    // Initialize soil conceptual reservoirs
    init_soil_reservoir(cfe_bmi_data_ptr, alpha_fc, max_soil_storage, S_soil, is_S_soil_ratio);

    // Initialize the runoff queue to empty to start with
    cfe_bmi_data_ptr->runoff_queue_m_per_timestep = malloc(sizeof(double) * cfe_bmi_data_ptr->num_giuh_ordinates + 1);
    for (i = 0; i < cfe_bmi_data_ptr->num_giuh_ordinates + 1; i++)
        cfe_bmi_data_ptr->runoff_queue_m_per_timestep[i] = 0.0;


    // jmframe: initialize zero potential and actual evapotranspiration.
    //          Will subtract this from the soil moisture and infiltration mass.
    cfe_bmi_data_ptr->et_struct.potential_et_m_per_s = 0;
    cfe_bmi_data_ptr->et_struct.actual_et_m_per_timestep = 0;

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

    // Two modes to get forcing data... 1) read from file, 2) pass with bmi    
    if (cfe_ptr->is_forcing_from_bmi)
      // BMI sets the precipitation to the aorc structure.
      cfe_ptr->timestep_rainfall_input_m = cfe_ptr->aorc.precip_kg_per_m2;
    else
      // Set the current rainfall input to the right place in the forcing array.
      cfe_ptr->timestep_rainfall_input_m = cfe_ptr->forcing_data_precip_kg_per_m2[cfe_ptr->current_time_step];
    
    cfe_ptr->vol_struct.volin += cfe_ptr->timestep_rainfall_input_m;
    // Delete me...    printf("_______PRECIP IN  CFE UPDATE FUNCTION________: %lf\n", cfe_ptr->aorc.precip_kg_per_m2);
    run_cfe(cfe_ptr);

    // Advance the model time 
    cfe_ptr->current_time_step += 1;

    return BMI_SUCCESS;
}


static int Update_until (Bmi *self, double t)
{
    // Since this model's time units are seconds, it is assumed that the param is either a valid time in seconds, a
    // relative number of time steps into the future, or invalid

    // Don't support negative parameter values
    if (t < 0.0)
        return BMI_FAILURE;

    // Don't continue if current time is at or beyond end time (or we can't determine this)
    double current_time, end_time;
    int current_time_result = self->get_current_time(self, &current_time);
    if (current_time_result == BMI_FAILURE)
        return BMI_FAILURE;
    int end_time_result = self->get_end_time(self, &end_time);
    if (end_time_result == BMI_FAILURE || current_time >= end_time)
        return BMI_FAILURE;

    // Handle easy case of t == current_time by just returning success
    if (t == current_time)
        return BMI_SUCCESS;

    cfe_state_struct* cfe_ptr = ((cfe_state_struct *) self->data);

    // First, determine if t is some future time that will be arrived at exactly after some number of future time steps
    int is_exact_future_time = (t == end_time) ? TRUE : FALSE;
    // Compare to time step endings unless obvious that t lines up (i.e., t == end_time) or doesn't (t <= current_time)
    if (is_exact_future_time == FALSE && t > current_time) {
        int future_time_step = cfe_ptr->current_time_step;
        double future_time_step_time = current_time;
        while (future_time_step < cfe_ptr->num_timesteps && future_time_step_time < end_time) {
            future_time_step_time += cfe_ptr->time_step_size;
            if (future_time_step_time == t) {
                is_exact_future_time = TRUE;
                break;
            }
        }
    }
    // If it is an exact time, advance to that time step
    if (is_exact_future_time == TRUE) {
        while (current_time < t) {
            
            run_cfe(cfe_ptr);
            
            // Advance the model time 
            cfe_ptr->current_time_step += 1;
        
            // Set the current rainfall input to the right place in the forcing.
            cfe_ptr->timestep_rainfall_input_m = cfe_ptr->forcing_data_precip_kg_per_m2[cfe_ptr->current_time_step];

            self->get_current_time(self, &current_time);

        }
        return BMI_SUCCESS;
    }

    // If t is not an exact time, it could be a number of time step forward to proceed

    // The model doesn't support partial time step value args (i.e., fractions)
    int t_int = (int) t;
    if ((t - ((double)t_int)) != 0)
        return BMI_FAILURE;

    // Keep in mind the current_time_step hasn't been processed yet (hence, using <= for this test)
    // E.g., if (unprocessed) current_time_step = 0, t = 2, num_timesteps = 2, this is valid a valid t (run 0, run 1)
    if ((cfe_ptr->current_time_step + t_int) <= cfe_ptr->num_timesteps) {
        for (i = 0; i < t_int; i++)
            
            // Set the current rainfall input to the right place in the forcing.
            cfe_ptr->timestep_rainfall_input_m = cfe_ptr->forcing_data_precip_kg_per_m2[i];

            run_cfe(cfe_ptr);

        return BMI_SUCCESS;
    }

    // If we arrive here, t wasn't an exact time at end of a time step or a valid relative time step jump, so invalid.
    return BMI_FAILURE;
}


static int Finalize (Bmi *self)
{
    // Function assumes everything that is needed is retrieved from the model before Finalize is called.
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

 //       if( model->flux_overland_m != NULL )    //NOT NEEDED redundant with flux_Schaake_runoff_m
 //           free(model->flux_overland_m);
        if( model->flux_Qout_m != NULL )
            free(model->flux_Qout_m);
        if( model->flux_Schaake_output_runoff_m != NULL )
            free(model->flux_Schaake_output_runoff_m);
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
    for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, output_var_names[i]) == 0)
            return i;

    for (i = 0; i < INPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, input_var_names[i]) == 0)
            return i + OUTPUT_VAR_NAME_COUNT;

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
    // If we get here, it means the variable name wasn't recognized
    type[0] = '\0';
    return BMI_FAILURE;
}


static int Get_var_itemsize (Bmi *self, const char *name, int * size)
{
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
    /***********************************************************/
    /***********    OUTPUT   ***********************************/
    /***********************************************************/
    if (strcmp (name, "RAIN_RATE") == 0) {   //jmframe: Seems unnecessary to have rain rate as an output variable.
        cfe_state_struct *cfe_ptr;
        cfe_ptr = (cfe_state_struct *) self->data;
        *dest = (void*)&cfe_ptr->aorc.precip_kg_per_m2;
        return BMI_SUCCESS;
    }

    if (strcmp (name, "SCHAAKE_OUTPUT_RUNOFF") == 0) {
        *dest = (void*) ((cfe_state_struct *)(self->data))->flux_Schaake_output_runoff_m;
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
    
    data->forcing_data_precip_kg_per_m2 = NULL;
    data->forcing_data_time = NULL;
    data->giuh_ordinates = NULL;
    data->nash_storage = NULL;
    data->runoff_queue_m_per_timestep = NULL;
    data->flux_Qout_m = NULL;
    data->flux_Schaake_output_runoff_m = NULL;
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
        cfe_ptr->Schaake_adjusted_magic_constant_by_soil_type,   // Set by config file
        cfe_ptr->timestep_rainfall_input_m,                      // Set by bmi (set value) or read from file.
        cfe_ptr->flux_Schaake_output_runoff_m,                  // Set by cfe function
        &cfe_ptr->infiltration_depth_m,                          // Set by Schaake partitioning scheme
//        cfe_ptr->flux_overland_m,                               // Set by CFE function, after Schaake  not needed, redundant with flux_Schaake_runoff_m
        &cfe_ptr->vol_struct.vol_sch_runoff,                     // Set by set_volume_trackers_to_zero
        &cfe_ptr->vol_struct.vol_sch_infilt,                     // Set by set_volume_trackers_to_zero
        cfe_ptr->flux_perc_m,                                   // Set to zero in definition.
        &cfe_ptr->vol_struct.vol_to_soil,                        // Set by set_volume_trackers_to_zero
        cfe_ptr->flux_lat_m,                                    // Set by CFE function after soil_resevroir calc
        &cfe_ptr->gw_reservoir_storage_deficit_m,                // Set by CFE function after soil_resevroir calc
        &cfe_ptr->gw_reservoir,      // Set in initialize and from config file
        &cfe_ptr->vol_struct.vol_to_gw,                          // Set by set_volume_trackers_to_zero
        &cfe_ptr->vol_struct.vol_soil_to_gw,                     // Set by set_volume_trackers_to_zero
        &cfe_ptr->vol_struct.vol_soil_to_lat_flow,               // Set by set_volume_trackers_to_zero
        &cfe_ptr->vol_struct.volout,                             // Set by set_volume_trackers_to_zero
        cfe_ptr->flux_from_deep_gw_to_chan_m,                   // Set by CFE function after gw_reservoir calc
        &cfe_ptr->vol_struct.vol_from_gw,                        // Set by set_volume_trackers_to_zero
        cfe_ptr->flux_giuh_runoff_m,                            // Set in CFE by convolution_integral
        cfe_ptr->num_giuh_ordinates,                             // Set by config file with func. count_delimited_values
        cfe_ptr->giuh_ordinates,                            // Set by configuration file.
        cfe_ptr->runoff_queue_m_per_timestep,               // Set in initialize
        &cfe_ptr->vol_struct.vol_out_giuh,                       // Set by set_volume_trackers_to_zero
        cfe_ptr->flux_nash_lateral_runoff_m,                    // Set in CFE from nash_cascade function
        cfe_ptr->num_lateral_flow_nash_reservoirs,               // Set from config file
        cfe_ptr->K_nash,                                         // Set from config file
        cfe_ptr->nash_storage,                              // Set from config file
        &cfe_ptr->vol_struct.vol_in_nash,                        // Set by set_volume_trackers_to_zero
        &cfe_ptr->vol_struct.vol_out_nash,                       // Set by set_volume_trackers_to_zero
        &cfe_ptr->et_struct,                                    // Set to zero with initalize. Set by BMI (set_value) during run
        cfe_ptr->flux_Qout_m                                    // Set by CFE function
    );
}

// Functions for setting up CFE data, i.e., initializing...
extern void init_soil_reservoir(cfe_state_struct* cfe_ptr, double alpha_fc, double max_storage, double storage,
                                int is_storage_ratios)
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
    cfe_ptr->soil_reservoir.exponent_primary = 1.0;
    // i.e., field_capacity_storage_threshold_m
    cfe_ptr->soil_reservoir.storage_threshold_primary_m =
            cfe_ptr->NWM_soil_params.smcmax * pow(1.0 / cfe_ptr->NWM_soil_params.satpsi, (-1.0 / cfe_ptr->NWM_soil_params.bb)) *
            (upper_lim - lower_lim);
    // lateral flow parameters --------------------------------------------------------
    // 0.0 to deactiv. else =lateral_flow_linear_reservoir_constant;   // m per ts
    // TODO: look at whether K_lf needs to be a derived (i.e., via get_K_lf_for_time_step()) or explicit parameter
    cfe_ptr->soil_reservoir.coeff_secondary = cfe_ptr->K_lf;
    // 1.0=linear
    cfe_ptr->soil_reservoir.exponent_secondary = 1.0;
    // making them the same, but they don't have 2B
    cfe_ptr->soil_reservoir.storage_threshold_secondary_m = cfe_ptr->soil_reservoir.storage_threshold_primary_m;
    cfe_ptr->soil_reservoir.storage_m = init_reservoir_storage(is_storage_ratios, storage, max_storage);
}

extern double init_reservoir_storage(int is_ratio, double amount, double max_amount) {
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
}

extern void initialize_volume_trackers(cfe_state_struct* cfe_ptr){
    cfe_ptr->vol_struct.vol_sch_runoff = 0;
    cfe_ptr->vol_struct.vol_sch_infilt = 0;
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
    printf("%16.8lf %lf %lf %lf %lf %lf %lf\n",
                           cfe_ptr->current_time_step,
                           cfe_ptr->timestep_rainfall_input_m*1000.0,
                           *cfe_ptr->flux_Schaake_output_runoff_m*1000.0,
                           *cfe_ptr->flux_giuh_runoff_m*1000.0,
                           *cfe_ptr->flux_nash_lateral_runoff_m*1000.0, 
                           *cfe_ptr->flux_from_deep_gw_to_chan_m*1000.0,
                           *cfe_ptr->flux_Qout_m*1000.0 );
}

extern void mass_balance_check(cfe_state_struct* cfe_ptr){
    //
    // PERFORM MASS BALANCE CHECKS AND WRITE RESULTS TO stderr.
    //----------------------------------------------------------
    
    double volend= cfe_ptr->soil_reservoir.storage_m+cfe_ptr->gw_reservoir.storage_m;
    double vol_in_gw_end = cfe_ptr->gw_reservoir.storage_m;
    double vol_end_giuh;
    double vol_in_nash_end;
    double vol_soil_end;
    
    // the GIUH queue might have water in it at the end of the simulation, so sum it up.
    for(i=0;i<cfe_ptr->num_giuh_ordinates;i++) vol_end_giuh+=cfe_ptr->runoff_queue_m_per_timestep[i];
    
    for(i=0;i<cfe_ptr->num_lateral_flow_nash_reservoirs;i++)  vol_in_nash_end+=cfe_ptr->nash_storage[i];
    
    
    vol_soil_end=cfe_ptr->soil_reservoir.storage_m;
    
    double global_residual;
    double schaake_residual;
    double giuh_residual;
    double soil_residual;
    double nash_residual;
    double gw_residual;
    
    global_residual = cfe_ptr->vol_struct.volstart + cfe_ptr->vol_struct.volin - cfe_ptr->vol_struct.volout - volend;
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
    
    schaake_residual = cfe_ptr->vol_struct.volin - cfe_ptr->vol_struct.vol_sch_runoff - cfe_ptr->vol_struct.vol_sch_infilt;
    printf(" SCHAAKE MASS BALANCE\n");
    printf("  surface runoff: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_runoff);
    printf("    infiltration: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_infilt);
    printf("schaake residual: %6.4e m\n",schaake_residual);  // should equal 0.0
    if(!is_fabs_less_than_epsilon(schaake_residual,1.0e-12))
                  printf("WARNING: SCHAAKE PARTITIONING MASS BALANCE CHECK FAILED\n");
    
    giuh_residual = cfe_ptr->vol_struct.vol_out_giuh - cfe_ptr->vol_struct.vol_sch_runoff - vol_end_giuh;
    printf(" GIUH MASS BALANCE\n");
    printf("  vol. into giuh: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_runoff);
    printf("   vol. out giuh: %8.4lf m\n",cfe_ptr->vol_struct.vol_out_giuh);
    printf(" vol. end giuh q: %8.4lf m\n",vol_end_giuh);
    printf("   giuh residual: %6.4e m\n",giuh_residual);  // should equal zero
    if(!is_fabs_less_than_epsilon(giuh_residual,1.0e-12))
                  printf("WARNING: GIUH MASS BALANCE CHECK FAILED\n");
    
    soil_residual=cfe_ptr->vol_struct.vol_soil_start + cfe_ptr->vol_struct.vol_sch_infilt -
                  cfe_ptr->vol_struct.vol_soil_to_lat_flow - vol_soil_end - cfe_ptr->vol_struct.vol_to_gw;
    printf(" SOIL WATER CONCEPTUAL RESERVOIR MASS BALANCE\n");
    printf("   init soil vol: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_start);     
    printf("  vol. into soil: %8.4lf m\n",cfe_ptr->vol_struct.vol_sch_infilt);
    printf("vol.soil2latflow: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_to_lat_flow);
    printf(" vol. soil to gw: %8.4lf m\n",cfe_ptr->vol_struct.vol_soil_to_gw);
    printf(" final vol. soil: %8.4lf m\n",vol_soil_end);   
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

    /*   JMFRAME: Something was wrong with the AORC reading code from the frameork, 
                  because I didn't have the correct libraries, that I assume are part of the framework.
                  so I paseted in the AORC reading code from the ET function.
    */

    char str[20];
    long yr, mo, da, hr, mi;
    double mm, julian, se;
    double val;
    int i, start, end, len;
    int yes_pm, wordlen;
    char theWord[150];

    len = strlen(theString);

    char *copy, *copy_to_free, *value;
    copy_to_free = copy = strdup(theString);

    // time
    value = strsep(&copy, ",");
    // TODO: handle this
    struct tm{
      int tm_year;
      int tm_mon; 
      int tm_mday; 
      int tm_hour; 
      int tm_min; 
      int tm_sec; 
      int tm_isdst;
    };
    struct tm t;
    time_t t_of_day;

    t.tm_year = (int)strtol(strsep(&value, "-"), NULL, 10) - 1900;
    t.tm_mon = (int)strtol(strsep(&value, "-"), NULL, 10);
    t.tm_mday = (int)strtol(strsep(&value, " "), NULL, 10);
    t.tm_hour = (int)strtol(strsep(&value, ":"), NULL, 10);
    t.tm_min = (int)strtol(strsep(&value, ":"), NULL, 10);
    t.tm_sec = (int)strtol(value, NULL, 10);
    t.tm_isdst = -1;        // Is DST on? 1 = yes, 0 = no, -1 = unknown
//    aorc->time = mktime(t);

    // APCP_surface
    value = strsep(&copy, ",");
    // Not sure what this is

    // DLWRF_surface
    value = strsep(&copy, ",");
    aorc->incoming_longwave_W_per_m2 = strtof(value, NULL);
    // DSWRF_surface
    value = strsep(&copy, ",");
    aorc->incoming_shortwave_W_per_m2 = strtof(value, NULL);
    // PRES_surface
    value = strsep(&copy, ",");
    aorc->surface_pressure_Pa = strtof(value, NULL);
    // SPFH_2maboveground
    value = strsep(&copy, ",");
    aorc->specific_humidity_2m_kg_per_kg = strtof(value, NULL);;
    // TMP_2maboveground
    value = strsep(&copy, ",");
    aorc->air_temperature_2m_K = strtof(value, NULL);
    // UGRD_10maboveground
    value = strsep(&copy, ",");
    aorc->u_wind_speed_10m_m_per_s = strtof(value, NULL);
    // VGRD_10maboveground
    value = strsep(&copy, ",");
    aorc->v_wind_speed_10m_m_per_s = strtof(value, NULL);
    // precip_rate
    value = strsep(&copy, ",");
    aorc->precip_kg_per_m2 = strtof(value, NULL);

    // Go ahead and free the duplicate copy now
    free(copy_to_free);

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
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit;
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (int **) malloc(frows * sizeof(int *));
    if (*array) {
        memset((*array), 0, frows * sizeof(int *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (int *) malloc(fcols * sizeof(int));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(int));
        }
    }
    return;
}


void dtwo_alloc_cfe(double ***array, int rows, int cols) {
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit;
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (double **) malloc(frows * sizeof(double *));
    if (*array) {
        memset((*array), 0, frows * sizeof(double *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (double *) malloc(fcols * sizeof(double));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
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
    static int ret[4];

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

    ret[3] = (j + 1L) % 7L;
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

