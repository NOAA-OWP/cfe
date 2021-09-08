#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../../include/bmi.h"
#include "../include/aorc.h"
#include "../include/bmi_aorc.h"

#define INPUT_VAR_NAME_COUNT 0 // This BMI reads from files, then passes this information to the Modules.
#define OUTPUT_VAR_NAME_COUNT 8 // All the forcings

static int 
Initialize (Bmi *self, const char *cfg_file)
{
    // don't print here. print in the framework.
    //printf("BMI Initialization AORC ... with config file %s \n", cfg_file);
        
    aorc_model *aorc;
    aorc = (aorc_model *) self->data;

    int config_read_result = read_init_config_aorc(aorc, cfg_file);
    if (config_read_result == BMI_FAILURE)
        return BMI_FAILURE;

    aorc_setup(aorc);

    if (aorc->bmi.verbose >1)
        printf("BMI Initialization AORC ... setup just finished \n");
    
    // Figure out the number of lines first (also char count)
    int forcing_line_count, max_forcing_line_length;

    int count_result = read_file_line_counts_aorc(aorc->forcing_file, &forcing_line_count, &max_forcing_line_length);
    if (count_result == -1) {
        printf("Configured forcing file '%s' could not be opened for reading\n", aorc->forcing_file);
        return BMI_FAILURE;
    }
    if (forcing_line_count == 1) {
        printf("Invalid header-only forcing file '%s'\n", aorc->forcing_file);
        return BMI_FAILURE;
    }

    // Now initialize empty arrays that depend on number of time steps
    aorc->forcing_data_precip_kg_per_m2 = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_surface_pressure_Pa = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_time = malloc(sizeof(long) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_incoming_shortwave_W_per_m2 = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_incoming_longwave_W_per_m2 = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_specific_humidity_2m_kg_per_kg = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_air_temperature_2m_K = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_u_wind_speed_10m_m_per_s = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));
    aorc->forcing_data_v_wind_speed_10m_m_per_s = malloc(sizeof(double) * (aorc->bmi.num_timesteps + 1));

    // Now open it again to read the forcings
    FILE* ffp = fopen(aorc->forcing_file, "r");
    // Ensure still exists
    if (ffp == NULL) {
        printf("Forcing file '%s' disappeared!", aorc->forcing_file);
        return BMI_FAILURE;
    }

    // Read forcing file and parse forcings
    char line_str[max_forcing_line_length + 1];
    long year, month, day, hour, minute;
    double dsec;
    // First read the header line
    fgets(line_str, max_forcing_line_length + 1, ffp);

    if (aorc->bmi.verbose > 2) 
        printf("the number of time steps from the forcing file is: %8.6e \n", aorc->bmi.num_timesteps);

    aorc_forcing_data forcings;
    for (int i = 0; i < aorc->bmi.num_timesteps; i++) {
        fgets(line_str, max_forcing_line_length + 1, ffp);  // read in a line of AORC data.
        parse_aorc_line_aorc(line_str, &year, &month, &day, &hour, &minute, &dsec, &forcings);
        aorc->forcing_data_precip_kg_per_m2[i] = forcings.precip_kg_per_m2 * ((double)aorc->bmi.time_step_size);
        if (aorc->bmi.verbose >4)
            printf("precip %f \n", aorc->forcing_data_precip_kg_per_m2[i]);
        aorc->forcing_data_surface_pressure_Pa[i] = forcings.surface_pressure_Pa;
        if (aorc->bmi.verbose >4)
            printf("surface pressure %f \n", aorc->forcing_data_surface_pressure_Pa[i]);
        aorc->forcing_data_incoming_longwave_W_per_m2[i] = forcings.incoming_longwave_W_per_m2;
        if (aorc->bmi.verbose >4)
            printf("longwave %f \n", aorc->forcing_data_incoming_longwave_W_per_m2[i]);
        aorc->forcing_data_incoming_shortwave_W_per_m2[i] = forcings.incoming_shortwave_W_per_m2;
        if (aorc->bmi.verbose >4)
            printf("shortwave %f \n", aorc->forcing_data_incoming_shortwave_W_per_m2[i]);
        aorc->forcing_data_specific_humidity_2m_kg_per_kg[i] = forcings.specific_humidity_2m_kg_per_kg;
        if (aorc->bmi.verbose >4)
            printf("humidity %f \n", aorc->forcing_data_specific_humidity_2m_kg_per_kg[i]);
        aorc->forcing_data_air_temperature_2m_K[i] = forcings.air_temperature_2m_K;
        if (aorc->bmi.verbose >4)
            printf("air temperature %f \n", aorc->forcing_data_air_temperature_2m_K[i]);
        aorc->forcing_data_u_wind_speed_10m_m_per_s[i] = forcings.u_wind_speed_10m_m_per_s;
        if (aorc->bmi.verbose >4)
            printf("u wind speed %f \n", aorc->forcing_data_u_wind_speed_10m_m_per_s[i]);
        aorc->forcing_data_v_wind_speed_10m_m_per_s[i] = forcings.v_wind_speed_10m_m_per_s;
        if (aorc->bmi.verbose >4)
            printf("v wind speed %f \n", aorc->forcing_data_v_wind_speed_10m_m_per_s[i]);


        aorc->forcing_data_time[i] = forcings.time;
        if (i == 0)
            aorc->bmi.current_time =forcings.time;
    }

    // Set the current time step to the first idem in the forcing time series.
    // But should this be an option? Would we ever initialize to a point in the
    //     middle of a forcing file?
    aorc->bmi.current_step = 0;

    return BMI_SUCCESS;
}

static int 
Update (Bmi *self)
{
    aorc_model *aorc = (aorc_model *) self->data;
  
    if (aorc->bmi.verbose >1)
      printf("BMI Update AORC ...\n");
  
    run_aorc(aorc);
    
    aorc->bmi.current_time_step += aorc->bmi.time_step_size; // Seconds since start of run
    aorc->bmi.current_step +=1;                            // time steps since start of run
    aorc->bmi.current_time += aorc->bmi.time_step_size;   // Seconds since 1970

    return BMI_SUCCESS;
}

static int
Update_until(Bmi *self, double t)
{

    aorc_model* aorc = ((aorc_model *) self->data);
    
    // Since this model's time units are seconds, it is assumed that the param is either a valid time in seconds, a
    // relative number of time steps into the future, or invalid

    if (aorc->bmi.verbose >= 1)
        printf("running update until time %lld \n", t);

    // Don't support negative parameter values
    if (t < 0.0)
        return BMI_FAILURE;

    // jmframe troubleshooting. delete this printout
    if (aorc->bmi.verbose > 3)
        printf("t is greater than zero, please delete this printout\n");

    // Don't continue if current time is at or beyond end time (or we can't determine this)
    double current_time = aorc->bmi.current_time_step;
    double end_time = 0;
    int current_time_result = self->get_current_time(self, &current_time);
    if (current_time_result == BMI_FAILURE)
        return BMI_FAILURE;

    // jmframe troubleshooting. delete this printout
    if (aorc->bmi.verbose > 3)
        printf("passed current time result, please delete this printout\n");

    int end_time_result = self->get_end_time(self, &end_time);

    // jmframe troubleshooting. Please delete this printout
    if (aorc->bmi.verbose > 3){
        printf("end time %8.6e, please delete this printout\n", &end_time);
        printf("the current time is %8.6e \n", current_time);
        printf("the end time result is %d \n", end_time_result);
    }

    if (end_time_result == BMI_FAILURE || current_time >= end_time){
        printf("BMI Failued due to end_time_result OR current time greater than end time \n");
        if (end_time_result == BMI_FAILURE)
            printf("end_time_result failed %s \n");
        if (current_time >= end_time)
            printf("current_time %8.6e is greater than end time %8.6e \n", current_time, end_time);
        return BMI_FAILURE;
    }

    // jmframe troubleshooting. Please delete this printout
    if (aorc->bmi.verbose > 3)
        printf("passed end time result, please delete this printout\n");

    // Handle easy case of t == current_time by just returning success
    if (t == current_time)
        return BMI_SUCCESS;

    // jmframe troubleshooting. Please delete this printout
    if (aorc->bmi.verbose > 3)
        printf("About to extract future time, please delete this printout\n");

    // First, determine if t is some future time that will be arrived at exactly after some number of future time steps
    int is_exact_future_time = (t == end_time) ? TRUE : FALSE;
    // Compare to time step endings unless obvious that t lines up (i.e., t == end_time) or doesn't (t <= current_time)
    if (is_exact_future_time == FALSE && t > current_time) {
        int future_time_step = aorc->bmi.current_time_step;
        double future_time_step_time = current_time;
        while (future_time_step < aorc->bmi.num_timesteps && future_time_step_time < end_time) {
            future_time_step_time += aorc->bmi.time_step_size;
            if (future_time_step_time == t) {
                is_exact_future_time = TRUE;
                break;
            }
        }
    }

    // jmframe troubleshooting. Please delete
    if (aorc->bmi.verbose > 3)
        printf("Extracted future time, please delete this printout\n");

    // If it is an exact time, advance to that time step
    if (is_exact_future_time == TRUE) {
        while (current_time < t) {
            
            self->update(self);
            
            // jmframe troubleshooting. Please delete
            if (aorc->bmi.verbose > 3){
                printf("current time step %f \n", aorc->bmi.current_time_step);
                printf("current step %f \n", aorc->bmi.current_step);
                printf("current time %f \n", aorc->bmi.current_time);
            }

        }
        return BMI_SUCCESS;
    }

    // jmframe troubleshooting. Please delete
    if (aorc->bmi.verbose > 3)
        printf("Checked if exact time, please delete this printout \n");

    // If t is not an exact time, it could be a number of time step forward to proceed

    // The model doesn't support partial time step value args (i.e., fractions)
    int t_int = (int) t;
    if ((t - ((double)t_int)) != 0)
        return BMI_FAILURE;

    // jmframe troubleshooting. Please delete
    if (aorc->bmi.verbose > 3)
        printf("Checked setting t_int to t, please delete this line \n");

    // Keep in mind the current_time_step hasn't been processed yet (hence, using <= for this test)
    // E.g., if (unprocessed) current_time_step = 0, t = 2, num_timesteps = 2, this is valid a valid t (run 0, run 1)
    if ((aorc->bmi.current_step + t_int) <= aorc->bmi.num_timesteps) {
        for (int i = 0; i < t_int; i++){

            // jmframe troubleshooting. Please delete
            if (aorc->bmi.verbose > 3)
                printf("Updating AORC forcing data with time step\n");
            
            // Call the update function, instead of run
            // because update also advances the bmi time
            self->update(self);

            // jmframe troubleshooting. Please delete
            if (aorc->bmi.verbose > 3){
                printf("current time step %f \n", aorc->bmi.current_time_step);
                printf("current step %ld \n", aorc->bmi.current_step);
                printf("current time %f \n", aorc->bmi.current_time);
            }
        } // end time loop
        return BMI_SUCCESS;
    }

    // If we arrive here, t wasn't an exact time at end of a time step or a valid relative time step jump, so invalid.
    return BMI_FAILURE;
    
}

aorc_model *
new_bmi_aorc()
{
    aorc_model *data;
    data = (aorc_model*) malloc(sizeof(aorc_model));

    return data;
}

static int 
Finalize (Bmi *self)
{

  // Perform Unit Tests, if set in config file.
  aorc_model *aorc = (aorc_model *) self->data;
  if (aorc->bmi.run_unit_tests == 1)
    aorc_unit_tests(aorc);

  if (self){
    aorc_model* model = (aorc_model *)(self->data);
    self->data = (void*)new_bmi_aorc();
  }
  return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------
static int 
Get_start_time (Bmi *self, double * time)
{
    // TODO: Get the real starting time, not zero. 
    // jmframe: I think this should be the time of the first line in the forcing file. In seconds since 1970
    *time = 0.0;
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {

    // Eight basic forcings. 
    "land_surface_radiation~incoming~longwave__energy_flux",
    "land_surface_air__pressure",
    "atmosphere_air_water~vapor__relative_saturation",
    "atmosphere_water__liquid_equivalent_precipitation_rate",
    "land_surface_radiation~incoming~shortwave__energy_flux",
    "land_surface_air__temperature",
    "land_surface_wind__x_component_of_velocity",
    "land_surface_wind__y_component_of_velocity"

    // Forcings and forcing derived values specifically used by Evapotranspiration
        

};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_types[OUTPUT_VAR_NAME_COUNT] = {
  "double",
  "double",
  "double",
  "double",
  "double",
  "double",
  "double",
  "double"
};

//---------------------------------------------------------------------------------------------------------------------
static const int output_var_item_count[OUTPUT_VAR_NAME_COUNT] = {
  1,
  1,
  1,
  1,
  1,
  1,
  1,
  1
};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_units[OUTPUT_VAR_NAME_COUNT] = {
  "W m-2",
  "Pa",
  "kg kg-1",
  "kg m-2",
  "W m-2",
  "K",
  "m s-1",
  "m s-1"
};

//---------------------------------------------------------------------------------------------------------------------
static const int output_var_grids[OUTPUT_VAR_NAME_COUNT] = {
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0
};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_locations[OUTPUT_VAR_NAME_COUNT] = {
        "node",
        "node",
        "node",
        "node",
        "node",
        "node",
        "node",
        "node"
};

//---------------------------------------------------------------------------------------------------------------------
static const int input_var_item_count[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------
static const char input_var_grids[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_locations[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------
// Don't forget to update Get_value/Get_value_at_indices 
// (and setter) implementation if these are adjusted
static const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_types[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_units[INPUT_VAR_NAME_COUNT] = {
};

//---------------------------------------------------------------------------------------------------------------------

static int Get_end_time (Bmi *self, double * time)
{
  
  Get_start_time(self, time);

  *time += (((aorc_model *) self->data)->bmi.num_timesteps * 
            ((aorc_model *) self->data)->bmi.time_step_size);

  return BMI_SUCCESS;
}


//---------------------------------------------------------------------------------------------------------------------

static int Get_time_step (Bmi *self, double * dt)
{
    *dt = ((aorc_model *) self->data)->bmi.time_step_size;
    return BMI_SUCCESS;
}

static int Get_time_units (Bmi *self, char * units)
{
    strncpy (units, "s", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------

static int Get_current_time (Bmi *self, double * time)
{
    Get_start_time(self, time);
    if (((aorc_model *) self->data)->bmi.verbose > 2){
        printf("Current model time step: '%ld'\n", ((aorc_model *) self->data)->bmi.current_time_step);
    }
    *time += (((aorc_model *) self->data)->bmi.current_step * 
              ((aorc_model *) self->data)->bmi.time_step_size);
    return BMI_SUCCESS;
} // end Get_current_time


//---------------------------------------------------------------------------------------------------------------------
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
} //end count_delimited_values

//---------------------------------------------------------------------------------------------------------------------
int read_file_line_counts_aorc(const char* file_name, int* line_count, int* max_line_length)
{

    *line_count = 0;
    *max_line_length = 0;
    int current_line_length = 0;
    FILE* fp = fopen(file_name, "r");
    // Ensure exists
    if (fp == NULL) {
        printf("File does not exist.\n Failed in function read_file_line_counts\n");
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
}  // end: read_file_line_counts

//---------------------------------------------------------------------------------------------------------------------
int read_init_config_aorc(aorc_model* model, const char* config_file)//,
{
    int config_line_count, max_config_line_length;
    // Note that this determines max line length including the ending return character, if present
    int count_result = read_file_line_counts_aorc(config_file, &config_line_count, &max_config_line_length);
    if (count_result == -1) {
        return BMI_FAILURE;
    }

    FILE* fp = fopen(config_file, "r");
    if (fp == NULL)
        return BMI_FAILURE;

    // TODO: document config file format (<param_key>=<param_val>, where array values are comma delim strings)

    char config_line[max_config_line_length + 1];

    for (int i = 0; i < config_line_count; i++) {
        char *param_key, *param_value;
        fgets(config_line, max_config_line_length + 1, fp);
        char* config_line_ptr = config_line;
        config_line_ptr = strsep(&config_line_ptr, "\n");
        param_key = strsep(&config_line_ptr, "=");
        param_value = strsep(&config_line_ptr, "=");

        if (strcmp(param_key, "verbose") == 0){
            model->bmi.verbose = strtod(param_value, NULL);
            if(model->bmi.verbose > 1){
                printf("printing some stuff (level > 1) for unit tests and troubleshooting \n");
            }
            if(model->bmi.verbose > 2){
                printf("printing a lot of stuff (level > 2) for unit tests and troubleshooting \n");
            }
        }
        if (strcmp(param_key, "forcing_file") == 0) {
            model->forcing_file = strdup(param_value);
            if(model->bmi.verbose >=2){
                printf("set forcing file from config file \n");
                printf("%s\n", model->forcing_file);
            }
            continue;
        }
        if (strcmp(param_key, "wind_speed_measurement_height_m") == 0) {
            model->aorc_params.wind_speed_measurement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("set wind speed measurement height from config file \n");
                printf("%lf\n", model->aorc_params.wind_speed_measurement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "humidity_measurement_height_m") == 0) {
            model->aorc_params.humidity_measurement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("set humidity measurement height from config file \n");
                printf("%lf\n", model->aorc_params.humidity_measurement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "vegetation_height_m") == 0) {
            model->aorc_params.vegetation_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("vegetation height from config file \n");
                printf("%lf\n", model->aorc_params.vegetation_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "zero_plane_displacement_height_m") == 0) {
            model->aorc_params.zero_plane_displacement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("zero_plane_displacement height from config file \n");
                printf("%lf\n", model->aorc_params.zero_plane_displacement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "momentum_transfer_roughness_length_m") == 0) {
            model->aorc_params.momentum_transfer_roughness_length_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("momentum_transfer_roughness_length_m from config file \n");
                printf("%lf\n", model->aorc_params.momentum_transfer_roughness_length_m);
            }
            continue;
        }
        if (strcmp(param_key, "surface_longwave_emissivity") == 0) {
            model->surf_rad_params.surface_longwave_emissivity = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("surface_longwave_emissivity from config file \n");
                printf("%lf\n", model->surf_rad_params.surface_longwave_emissivity);
            }
            continue;
        }
        if (strcmp(param_key, "surface_shortwave_albedo") == 0) {
            model->surf_rad_params.surface_shortwave_albedo = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("surface_shortwave_albedo from config file \n");
                printf("%lf\n", model->surf_rad_params.surface_shortwave_albedo);
            }
            continue;
        }
        if (strcmp(param_key, "surface_shortwave_albedo") == 0) {
            model->surf_rad_params.surface_shortwave_albedo = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("surface_shortwave_albedo from config file \n");
                printf("%lf\n", model->surf_rad_params.surface_shortwave_albedo);
            }
            continue;
        }
        if (strcmp(param_key, "latitude_degrees") == 0) {
            model->solar_params.latitude_degrees = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("latitude_degrees from config file \n");
                printf("%lf\n", model->solar_params.latitude_degrees);
            }
            continue;
        }
        if (strcmp(param_key, "longitude_degrees") == 0) {
            model->solar_params.longitude_degrees = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("longitude_degrees from config file \n");
                printf("%lf\n", model->solar_params.longitude_degrees);
            }
            continue;
        }
        if (strcmp(param_key, "site_elevation_m") == 0) {
            model->solar_params.site_elevation_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("site_elevation_m from config file \n");
                printf("%lf\n", model->solar_params.site_elevation_m);
            }
            continue;
        }
        if (strcmp(param_key, "time_step_size") == 0) {
            model->bmi.time_step_size = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("time_step_size from config file \n");
                printf("%d\n", model->bmi.time_step_size);
            }
            continue;
        }
        if (strcmp(param_key, "num_timesteps") == 0) {
            model->bmi.num_timesteps = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("num_timesteps from config file \n");
                printf("%d\n", model->bmi.num_timesteps);
            }
            continue;
        }
        if (strcmp(param_key, "run_unit_tests") == 0) {
            model->bmi.run_unit_tests = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("Running unit tests \n");
            }
            continue;
        }

    } // end loop through config

    return BMI_SUCCESS;
} // end: read_init_config

static int Get_var_type (Bmi *self, const char *name, char * type)
{
    // Check to see if in output array first
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(type, output_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(type, input_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    printf("Error: The item was not recognized when trying to get type\n");
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
    else if (strcmp (type, "float") == 0) {
        *size = sizeof(float);
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



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
static int Get_var_location (Bmi *self, const char *name, char * location)
{
    // Check to see if in output array first
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(location, output_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(location, input_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    location[0] = '\0';
    return BMI_FAILURE;
}

//------------------------------------------------------------------------------
static int Get_var_grid(Bmi *self, const char *name, int *grid)
{

    // Check to see if in output array first
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            *grid = output_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            *grid = input_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    grid[0] = '\0';
    return BMI_FAILURE;
}

// ***********************************************************
// ********* BMI: VARIABLE GETTER & SETTER FUNCTIONS *********
// ***********************************************************
static int Get_value_ptr (Bmi *self, const char *name, void **dest)
{
        aorc_model *aorc;
        aorc = (aorc_model *) self->data;
        void *src = NULL;

    if (strcmp (name, "land_surface_radiation~incoming~longwave__energy_flux") == 0) {
        *dest = (void*)&aorc->aorc.incoming_longwave_W_per_m2;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_radiation~incoming~shortwave__energy_flux") == 0) {
        *dest = (void*)&aorc->aorc.incoming_shortwave_W_per_m2;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_air__pressure") == 0) {
        *dest = (void*)&aorc->aorc.surface_pressure_Pa;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "atmosphere_air_water~vapor__relative_saturation") == 0) {
        *dest = (void*)&aorc->aorc.specific_humidity_2m_kg_per_kg;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "atmosphere_water__liquid_equivalent_precipitation_rate") == 0) {
        src = (void*)&aorc->aorc.precip_kg_per_m2;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_air__temperature") == 0) {
        *dest = (void*)&aorc->aorc.air_temperature_2m_K;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_wind__x_component_of_velocity") == 0) {
        *dest = (void*)&aorc->aorc.u_wind_speed_10m_m_per_s;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_wind__y_component_of_velocity") == 0) {
        *dest = (void*)&aorc->aorc.v_wind_speed_10m_m_per_s;
        return BMI_SUCCESS;
    }

    return BMI_FAILURE;
}

//------------------------------------------------------------------------------
static int Get_value_at_indices (Bmi *self, const char *name, void *dest, int * inds, int len)
{
    void *src = NULL;
    int itemsize = 0;

    if (self->get_value_ptr(self, name, &src) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_itemsize(self, name, &itemsize) == BMI_FAILURE)
        return BMI_FAILURE;

    { /* Copy the data */
        size_t i;
        size_t offset;
        char * ptr;
        for (i=0, ptr=(char*)dest; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy (ptr, (char*)src + offset, itemsize);
        }
    }

    return BMI_SUCCESS;
}

//------------------------------------------------------------------------------
static int Get_value(Bmi * self, const char * name, void *dest)
{
    void *src = NULL;
    int nbytes = 0;
 
    if (self->get_value_ptr (self, name, &src) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_nbytes (self, name, &nbytes) == BMI_FAILURE)
        return BMI_FAILURE;
 
    memcpy(dest, src, nbytes);

    return BMI_SUCCESS;
}

static int Set_value (Bmi *self, const char *name, void *array)
{
    void * dest = NULL;
    int nbytes = 0;

    if (self->get_value_ptr(self, name, &dest) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_nbytes(self, name, &nbytes) == BMI_FAILURE)
        return BMI_FAILURE;

    memcpy (dest, array, nbytes);

    return BMI_SUCCESS;
}


//------------------------------------------------------------------------------
static int Set_value_at_indices (Bmi *self, const char *name, int * inds, int len, void *src)
{
    void * to = NULL;
    int itemsize = 0;

    if (self->get_value_ptr (self, name, &to) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_itemsize(self, name, &itemsize) == BMI_FAILURE)
        return BMI_FAILURE;

    { /* Copy the data */
        size_t i;
        size_t offset;
        char * ptr;
        for (i=0, ptr=(char*)src; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy ((char*)to + offset, ptr, itemsize);
        }
    }
    return BMI_SUCCESS;
}

// ***********************************************************
// ************ BMI: MODEL INFORMATION FUNCTIONS *************
// ***********************************************************
static int Get_component_name (Bmi *self, char * name)
{
    strncpy (name, "AORC Forcing Data", BMI_MAX_COMPONENT_NAME);
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_input_item_count (Bmi *self, int * count)
{
    *count = INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_input_var_names (Bmi *self, char ** names)
{
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        strncpy (names[i], input_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_output_item_count (Bmi *self, int * count)
{
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_output_var_names (Bmi *self, char ** names)
{
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        strncpy (names[i], output_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_var_units (Bmi *self, const char *name, char * units)
{
    // Check to see if in output array first
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(units, output_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(units, input_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    units[0] = '\0';
    return BMI_FAILURE;
}

//----------------------------------------------------------------------
static int Get_var_nbytes (Bmi *self, const char *name, int * nbytes)
{
    int item_size;
    int item_size_result = Get_var_itemsize(self, name, &item_size);
    if (item_size_result != BMI_SUCCESS) {
        return BMI_FAILURE;
    }
    int item_count = -1;
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            item_count = input_var_item_count[i];
            break;
        }
    }
    if (item_count < 1) {
        for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
            if (strcmp(name, output_var_names[i]) == 0) {
                item_count = output_var_item_count[i];
                break;
            }
        }
    }
    if (item_count < 1)
        item_count = ((aorc_model *) self->data)->bmi.num_timesteps;

    *nbytes = item_size * item_count;
    return BMI_SUCCESS;
}


// ***********************************************************
// **************** BMI: MODEL GRID FUNCTIONS ****************
// ***********************************************************
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

//----------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------------------------------
/* Uniform rectilinear (grid type) */
static int Get_grid_shape(Bmi *self, int grid, int *shape)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_spacing(Bmi *self, int grid, double *spacing)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_origin(Bmi *self, int grid, double *origin)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
/* Non-uniform rectilinear, curvilinear (grid type)*/
static int Get_grid_x(Bmi *self, int grid, double *x)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_y(Bmi *self, int grid, double *y)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_z(Bmi *self, int grid, double *z)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
/*Unstructured (grid type)*/
static int Get_grid_node_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_edge_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_edge_nodes(Bmi *self, int grid, int *edge_nodes)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_edges(Bmi *self, int grid, int *face_edges)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_nodes(Bmi *self, int grid, int *face_nodes)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_nodes_per_face(Bmi *self, int grid, int *nodes_per_face)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
Bmi* 
register_bmi_aorc(Bmi *model)
{
    if (model) {
        model->data = (void*)new_bmi_aorc();
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
        model->get_value_ptr = Get_value_ptr;   
        model->get_value_at_indices = Get_value_at_indices;

        model->set_value = Set_value;
        model->set_value_at_indices = Set_value_at_indices;

        model->get_grid_size = Get_grid_size;    
        model->get_grid_rank = Get_grid_rank;    
        model->get_grid_type = Get_grid_type;    
        model->get_grid_shape = Get_grid_shape;
        model->get_grid_spacing = Get_grid_spacing;
        model->get_grid_origin = Get_grid_origin;
        model->get_grid_x = Get_grid_x;
        model->get_grid_y = Get_grid_y;
        model->get_grid_z = Get_grid_z;
        model->get_grid_node_count = Get_grid_node_count;
        model->get_grid_edge_count = Get_grid_edge_count;
        model->get_grid_face_count = Get_grid_face_count;
        model->get_grid_edge_nodes = Get_grid_edge_nodes;
        model->get_grid_face_edges = Get_grid_face_edges;
        model->get_grid_face_nodes = Get_grid_face_nodes;
        model->get_grid_nodes_per_face = Get_grid_nodes_per_face;

    }

    return model;
}  // end: register_bmi_aorc
