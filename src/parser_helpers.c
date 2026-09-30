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
 
 // CFE v3 config file parser (keyword-based format).
 // Legacy v2 configs must be converted using the cfe_migrate_config utility.
 // FLO 9/2025

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser_helpers.h"
#include "cfe_config.h"
#include "cfe.h"  /* physics constants: GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2, etc. */

/* Safe string copy with guaranteed null termination.
 * Used instead of snprintf(dst, sz, "%s", src) because GCC -Wformat-truncation
 * warns when the source buffer is larger than the destination, even though
 * snprintf handles truncation safely. strncpy avoids the diagnostic. */
static inline void safe_strcpy(char *dst, size_t dst_size, const char *src) {
    strncpy(dst, src, dst_size - 1);
    dst[dst_size - 1] = '\0';
}

// Function to trim whitespace from a string
//###################
char* trim_whitespace(char* str) {
    if (!str) return NULL;

    char* end;

    while (isspace((unsigned char)*str)) str++;

    if (*str == 0) return str;

    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;

    end[1] = '\0';
    return str;
}


// Function to extract units from bracketed section
//################
void extract_units(const char* line, char* units, size_t units_size) {
    const char* start = strchr(line, '[');
    const char* end = strchr(line, ']');
    
    if (start && end && end > start) {
        start++; // Move past the '['
        size_t len = end - start;
        if (len < units_size - 1) {
            strncpy(units, start, len);
            units[len] = '\0';
            trim_whitespace(units);
        }
    } else {
        units[0] = '\0'; // No units found
    }
}

// Function to parse array values (comma-separated)
//####################
int parse_double_array(const char* value_str, double* array, int max_elements) {
    size_t len = strlen(value_str);
    char* str_copy = (char*)malloc(len + 1);
    if (!str_copy) return 0;

    memcpy(str_copy, value_str, len);
    str_copy[len] = '\0';

    char* token = strtok(str_copy, ",");
    int count = 0;

    while (token && count < max_elements) {
        char* trimmed = trim_whitespace(token);
        char* endptr;
        double value = strtod(trimmed, &endptr);
        if (endptr == trimmed || *endptr != '\0') {
            fprintf(stderr, "ERROR: Invalid number format: '%s'\n", trimmed);
            free(str_copy);
            return -1; // Return error code
        }
        array[count] = value; 
        count++;
        token = strtok(NULL, ",");
    }

    if (token != NULL) {
        fprintf(stderr, "Warning: parse_double_array truncated input after %d elements.\n", max_elements);
    }

    free(str_copy);
    return count;
}

// Portable case-insensitive string comparison
//############################
int string_compare_ignore_case(const char* str1, const char* str2) {
    while (*str1 && *str2) {
        char c1 = (*str1 >= 'A' && *str1 <= 'Z') ? *str1 + 32 : *str1;
        char c2 = (*str2 >= 'A' && *str2 <= 'Z') ? *str2 + 32 : *str2;
        if (c1 != c2) return c1 - c2;
        str1++;
        str2++;
    }
    return *str1 - *str2;
}

// Function to parse boolean values
//###############
int parse_boolean(const char* value_str) {
    char* trimmed = trim_whitespace((char*)value_str);
    if (string_compare_ignore_case(trimmed, "TRUE") == 0 || 
        string_compare_ignore_case(trimmed, "1") == 0) {
        return 1;
    }
    return 0;
}

// Function to clean quoted strings (remove quotes iff present)
//######################
void clean_quoted_string(char* str) {
    size_t len = strlen(str);
    if (len >= 2 && str[0] == '"' && str[len-1] == '"') {
        // Remove quotes
        memmove(str, str + 1, len - 2);
        str[len - 2] = '\0';
    }
}

// Helper function to make sure that N values provided when N required.
//######################
int validate_array_count(const char* keyword, int expected_count, int actual_count, 
                        double* array, int allow_single_zero) {
    
    // Special case: if only one value provided and it's 0.0, fill all with 0.0
    if (allow_single_zero && actual_count == 1 && array[0] == 0.0) {
        for (int i = 1; i < expected_count; i++) {
            array[i] = 0.0;
        }
        return expected_count; // Return corrected count
    }
    
    // Otherwise, require exact match
    if (actual_count != expected_count) {
        fprintf(stderr, "ERROR: %s requires exactly %d values, got %d\n", 
                keyword, expected_count, actual_count);
        return -1;
    }
    
    return actual_count;
}

// Returns cfe version iff CFE_CONFIG_VERSION= keyword string is found, ellse returns 0.0
//############################
double read_cfe_config_version(const char* cfg_path)
{
    FILE* f = fopen(cfg_path, "r");
    if (!f) {
        perror(cfg_path);
        return 0.0;
    }

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "cfe_config_version")) {
            double version;
            if (sscanf(line, "cfe_config_version = %lf", &version) == 1) {
                fclose(f);
                return version;
            }
        }
    }

    fclose(f);
    return 0.0;  // No version key found
}

// CFE v3 config file parser
//##########################
int parse_cfe_config(const char* filename, CFE_CONFIG* config) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return -1;
    }

    PARSER_ARRAY_COUNTS array_counts = {0};
    
    memset(config, 0, sizeof(CFE_CONFIG));
    config->catchment_vegetated_fraction = 1.0;
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {  // or whatever max for subsurface
        config->subsurface_routing_nash_cascade_init_storage_m[i] = 0.0;
    }
    char line[1024];
    char keyword[256];
    char value_part[768];
    char units[64];

    snprintf(config->surface_routing_scheme_name, sizeof(config->surface_routing_scheme_name), "%s", "giuh");  // default and only option
    
    while (fgets(line, sizeof(line), file)) {
        // Skip empty lines and comments
        if (line[0] == '\n' || line[0] == '#') continue;
        char* line_start = line;
        while (*line_start && isspace(*line_start)) line_start++;   // skip leading white spaces   
        char* comment_pos = strchr(line_start, '#');  // strip comments that might exist at the end of a keyword=val line.
        char* cpp_comment_pos = strstr(line_start, "//");
        if (comment_pos && cpp_comment_pos) {
            // Use whichever comes first
            if (cpp_comment_pos < comment_pos) comment_pos = cpp_comment_pos;
        } else if (cpp_comment_pos) {
            comment_pos = cpp_comment_pos;
        }
        if (comment_pos) *comment_pos = '\0';
        // Find the equals sign
        char* equals_pos = strchr(line_start, '=');
        if (!equals_pos) continue;
        // Extract keyword: start from beginning, stop at equals, trim whitespace
        char* keyword_start = line_start;
        char* keyword_end = equals_pos;
        // Skip leading whitespace in keyword
        while (keyword_start < keyword_end && isspace(*keyword_start)) {
            keyword_start++;
        }
        // Skip trailing whitespace in keyword
        while (keyword_end > keyword_start && isspace(*(keyword_end - 1))) {
            keyword_end--;
        }

        // Copy keyword
        size_t keyword_len = keyword_end - keyword_start;
        if (keyword_len >= sizeof(keyword)) keyword_len = sizeof(keyword) - 1;
        strncpy(keyword, keyword_start, keyword_len);
        keyword[keyword_len] = '\0';

        // Extract value part (everything after '=')
        char* value_start = equals_pos + 1;
        while (*value_start && isspace(*value_start)) value_start++;  

        // Handle units in brackets [unit]
        char* bracket_pos = strchr(value_start, '[');
        if (bracket_pos) {
            size_t value_len = bracket_pos - value_start;
            strncpy(value_part, value_start, value_len);
            value_part[value_len] = '\0';
        } else {
            strncpy(value_part, value_start, sizeof(value_part) - 1);
            value_part[sizeof(value_part) - 1] = '\0';
            // Remove newline if present
            char* newline = strchr(value_part, '\n');
            if (newline) *newline = '\0';
        }

        // Trim whitespace from value
        trim_whitespace(value_part);

        // Extract units
        extract_units(line, units, sizeof(units));

        // Process keywords (same as before)
        if (string_compare_ignore_case(keyword, "cfe_config_version") == 0) {
            config->version = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_model_timestep_h") == 0) {
            config->timestep_h = atof(value_part);
            snprintf(config->timestep_units, sizeof(config->timestep_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_id") == 0) {
            safe_strcpy(config->cat_id, sizeof(config->cat_id), value_part);
        }
        // Catchment metadata — debug/diagnostic output only, not used by the model.
        else if (string_compare_ignore_case(keyword, "catchment_latitude_decimal_degree") == 0) {
            config->cat_latitude = atof(value_part);
            snprintf(config->cat_latitude_units, sizeof(config->cat_latitude_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_longitude_decimal_degree") == 0) {
            config->cat_longitude = atof(value_part);
            snprintf(config->cat_longitude_units, sizeof(config->cat_longitude_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_elevation") == 0) {
            config->cat_elev = atof(value_part);
            snprintf(config->cat_elev_units, sizeof(config->cat_elev_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_area_km2") == 0) {
            config->cat_area_km2 = atof(value_part);
            snprintf(config->cat_area_units, sizeof(config->cat_area_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_impervious_fraction_0-1") == 0) {
            config->cat_impervious_fraction = atof(value_part);
            snprintf(config->cat_impervious_units, sizeof(config->cat_impervious_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "control_input_forcing_filename") == 0) {
            safe_strcpy(config->control_input_forcing_filename, sizeof(config->control_input_forcing_filename), value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_total_num_simulation_timesteps") == 0) {
            config->total_timesteps = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_verbosity") == 0) {
            config->verbosity = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_simulation_start_date") == 0) {
            safe_strcpy(config->simulation_start_date, sizeof(config->simulation_start_date), value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_depth_m") == 0) {
            config->soil_depth_m = atof(value_part);
            snprintf(config->soil_depth_units, sizeof(config->soil_depth_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_Clapp_Hornberger_exponent_b") == 0) {
            config->soil_Clapp_Hornberger_exponent_b = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_sat_hydraulic_conductivity_cm_per_h") == 0) {
            config->soil_sat_hydraulic_conductivity_cm_per_h = atof(value_part);
            snprintf(config->soil_sat_hydraulic_conductivity_units, sizeof(config->soil_sat_hydraulic_conductivity_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_sat_capillary_head_cm") == 0) {
            config->soil_sat_capillary_head_cm = atof(value_part);
            snprintf(config->soil_sat_capillary_head_units, sizeof(config->soil_sat_capillary_head_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_to_gw_percolation_rate_limiter_0_to_1") == 0) {
            config->soil_to_gw_percolation_rate_limiter_0_to_1 = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_effective_porosity") == 0) {
            config->soil_effective_porosity = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_field_capacity_Pcap_over_Patm_0_1") == 0) {
            config->soil_field_capacity_Pcap_over_Patm_0_1 = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_soil_reservoir_init_storage_m") == 0) {
            config->soil_reservoir_init_storage_m = atof(value_part);
            snprintf(config->soil_reservoir_init_storage_units, sizeof(config->soil_reservoir_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_soil_reservoir_init_discrete_storage_theta") == 0) {
            int parsed_count = parse_double_array(value_part, config->soil_reservoir_init_discrete_storage_theta, NDISCS);

            if (parsed_count == -1) {
                fprintf(stderr, "ERROR: Failed to parse soil_reservoir_init_discrete_storage_theta\n");
                return -1; // parsing error
            }

            if (parsed_count != NDISCS) {
                fprintf(stderr, "ERROR: soil_reservoir_init_discrete_storage_theta requires exactly %d values, got %d\n", 
                        NDISCS, parsed_count);
                return -1; // wrong number of elements
            }

            // Optional: validate theta values are reasonable
            for (int i = 0; i < NDISCS; i++) {
                if (config->soil_reservoir_init_discrete_storage_theta[i] < 0.0 || 
                    config->soil_reservoir_init_discrete_storage_theta[i] > 1.0) {
                    fprintf(stderr, "WARNING: soil_reservoir_init_discrete_storage_theta[%d] = %.3f is outside typical range [0.0, 1.0]\n", 
                            i, config->soil_reservoir_init_discrete_storage_theta[i]);
                }
            }

            snprintf(config->soil_reservoir_init_discrete_storage_theta_units, 
                     sizeof(config->soil_reservoir_init_discrete_storage_theta_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_reservoir_rate_const_to_subsurface_lateral_flow") == 0) {
            config->soil_reservoir_rate_const_to_subsurface_lateral_flow = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_ice_content_impervious_threshold") == 0) {
            config->soil_ice_content_impervious_threshold = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "gw_reservoir_max_storage_m") == 0) {
            config->gw_reservoir_max_storage_m = atof(value_part);
            snprintf(config->gw_reservoir_max_storage_units, sizeof(config->gw_reservoir_max_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_gw_reservoir_init_storage_m") == 0) {
            config->gw_reservoir_init_storage_m = atof(value_part);
            snprintf(config->gw_reservoir_init_storage_units, sizeof(config->gw_reservoir_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "gw_discharge_coeff_m_per_timestep") == 0) {
            config->gw_discharge_coeff_m_per_timestep = atof(value_part);
            snprintf(config->gw_discharge_coeff_m_per_timestep_units, sizeof(config->gw_discharge_coeff_m_per_timestep_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "gw_discharge_exponent") == 0) {
            config->gw_discharge_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "subsurface_routing_nash_reservoir_time_constant_k") == 0) {
            config->subsurface_routing_nash_K = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_subsurface_routing_init_nash_cascade_storage_m") == 0) {
            array_counts.num_subsurf_nash_storages_read = parse_double_array(value_part, 
                                                                             config->subsurface_routing_nash_cascade_init_storage_m, 
                                                                             MAX_NUM_SUBSURFACE_NASH_CASCADE);
            snprintf(config->subsurface_routing_nash_cascade_init_storage_units, sizeof(config->subsurface_routing_nash_cascade_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_simulate_freeze_thaw_true_false") == 0) {
            config->control_soil_simulate_freeze_thaw_true_false = parse_boolean(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_simulate_discrete_soil_moisture_true_false") == 0) {
            config->control_soil_simulate_discrete_soil_moisture_true_false = parse_boolean(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_use_lookup_table_num_points") == 0) {
            config->control_soil_use_lookup_table_num_points = atoi(value_part);
            if(config->verbosity > 1) printf("DEBUG: Parsed control_soil_use_lookup_table_num_points = %s -> %d\n", 
                                             value_part, config->control_soil_use_lookup_table_num_points);
        }
        else if (string_compare_ignore_case(keyword, "control_ET_simulate_Priestley_Taylor") == 0) {
            config->et_alpha_pt = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_ET_deepest_root_zone_discretization") == 0) {
            config->control_ET_deepest_root_zone_discretization = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_simulate_soil_evaporation") == 0) {
            config->control_soil_simulate_soil_evaporation = parse_boolean(value_part);
        }
        else if (string_compare_ignore_case(keyword, "catchment_vegetated_fraction_0-1") == 0 ||
                 string_compare_ignore_case(keyword, "catchment_forested_fraction_0-1") == 0) {
            config->catchment_vegetated_fraction = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_scheme_name") == 0) {
            // Convert to lowercase for consistent comparison
            char temp_value[sizeof(config->partitioning_scheme_name)];
            safe_strcpy(temp_value, sizeof(temp_value), value_part);
            for (int i = 0; temp_value[i]; i++) {
                temp_value[i] = tolower(temp_value[i]);
            }
            snprintf(config->partitioning_scheme_name, sizeof(config->partitioning_scheme_name), "%s", temp_value);
        }
//
//        else if (string_compare_ignore_case(keyword, "surface_routing_scheme_name") == 0) {
//            // Convert to lowercase for consistent comparison
//            char temp_value[sizeof(config->surface_routing_scheme_name)];
//            safe_strcpy(temp_value, sizeof(temp_value), value_part);
//            for (int i = 0; temp_value[i]; i++) {
//                temp_value[i] = tolower(temp_value[i]);
//            }
//            snprintf(config->surface_routing_scheme_name, sizeof(config->surface_routing_scheme_name), "%s", temp_value);
//        }
        else if (string_compare_ignore_case(keyword, "surface_routing_num_giuh_ordinates") == 0) {
            config->surface_routing_num_giuh_ordinates = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "surface_routing_giuh_ordinates") == 0) {
            array_counts.num_giuh_ordinates_read = parse_double_array(value_part,
                                                                      config->surface_routing_giuh_ordinates,
                                                                      MAX_NUM_GIUH_ORDINATES);
            snprintf(config->surface_routing_giuh_ordinates_units,
                     sizeof(config->surface_routing_giuh_ordinates_units),
                     "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_surface_routing_init_giuh_convolution_queue_m") == 0) {
            array_counts.num_giuh_convolution_read = parse_double_array(value_part,
                               config->surface_routing_init_giuh_convolution_queue_m,
                               MAX_NUM_GIUH_ORDINATES);
            snprintf(config->surface_routing_init_giuh_convolution_queue_units,
                     sizeof(config->surface_routing_init_giuh_convolution_queue_units),
                     "%s", units);
        }

        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_tension_water_inflection_point") == 0) {
            config->soil_Xinanjiang_tension_water_inflection_point = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent") == 0) {
            config->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent") == 0) {
            config->soil_Xinanjiang_free_water_soil_moist_distrib_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_skin_temperature_k") == 0) {
            config->state_skin_temperature_k = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_upper_soil_temperature_k") == 0) {
            config->state_upper_soil_temperature_k = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_estimated_annual_air_temperature_k") == 0) {
            config->state_estimated_annual_air_temperature_k = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_air_temperature_time_integral_k_s") == 0) {
            config->state_air_temperature_time_integral_k_s = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_accumulated_time_s") == 0) {
            config->state_accumulated_time_s = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_pet_initialized") == 0) {
            config->state_pet_initialized = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_status_warnings_filename") == 0) {
            safe_strcpy(config->output_status_warnings_filename, sizeof(config->output_status_warnings_filename), value_part);
            clean_quoted_string(config->output_status_warnings_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_internal_fluxes_m_per_timestep_filename") == 0) {
            safe_strcpy(config->output_internal_fluxes_filename, sizeof(config->output_internal_fluxes_filename), value_part);
            clean_quoted_string(config->output_internal_fluxes_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_internal_storages_m_per_timestep_filename") == 0) {
            safe_strcpy(config->output_internal_storages_filename, sizeof(config->output_internal_storages_filename), value_part);
            clean_quoted_string(config->output_internal_storages_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_volume_balance_filename") == 0) {
            safe_strcpy(config->output_volume_balance_filename, sizeof(config->output_volume_balance_filename), value_part);
            clean_quoted_string(config->output_volume_balance_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_soil_moisture_theta_filename") == 0) {
            safe_strcpy(config->output_soil_moisture_theta_filename, sizeof(config->output_soil_moisture_theta_filename), value_part);
            clean_quoted_string(config->output_soil_moisture_theta_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_discharge_m_per_timestep_filename") == 0) {
            safe_strcpy(config->output_discharge_filename, sizeof(config->output_discharge_filename), value_part);
            clean_quoted_string(config->output_discharge_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_total_discharge_m3_per_sec_filename") == 0) {
            snprintf(config->output_total_discharge_m3_per_sec_filename, 
                     sizeof(config->output_total_discharge_m3_per_sec_filename), "%s", value_part);
            clean_quoted_string(config->output_total_discharge_m3_per_sec_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_path_name") == 0) {
            safe_strcpy(config->output_path_name, sizeof(config->output_path_name), value_part);
            clean_quoted_string(config->output_path_name);
        }
        else if (string_compare_ignore_case(keyword, "output_time_standard_format") == 0) {
            safe_strcpy(config->output_time_standard_format, sizeof(config->output_time_standard_format), value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_file_delimiter") == 0) {
            safe_strcpy(config->output_file_delimiter, sizeof(config->output_file_delimiter), value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_value_format") == 0) {
            safe_strcpy(config->output_value_format, sizeof(config->output_value_format), value_part);
            clean_quoted_string(config->output_value_format);

            // ensure that the specified output value format is valid- "%.Ne" or "%.Nf" where 0<N<16
            // iff not, set to defaultt "%.8f"
            if (!validate_and_fix_output_format(config->output_value_format, sizeof(config->output_value_format))) {
                if (config->verbosity > 0) {
                    fprintf(stderr, "WARNING: Invalid output_value_format, using default %%.8f\n");
                }
            }
        }
        else if (string_compare_ignore_case(keyword, "output_new_config_filename_prefix") == 0) {
            safe_strcpy(config->output_new_config_filename, sizeof(config->output_new_config_filename), value_part);
            clean_quoted_string(config->output_new_config_filename);
            
        }

    }

    fclose(file);

    /* Default version if not specified in config */
    if (config->version < 1.0e-04)
        config->version = 3.0;

    // Apply default output_value_format if the config file did not specify one.
    // validate_and_fix_output_format() only runs when the "output_value_format"
    // keyword is present but malformed; if the keyword is absent entirely, this
    // field is left as an empty string by the memset() at the top of the parser.
    // An empty format string passed to fprintf() silently prints nothing for
    // every numeric field in q.out/Q.out/fluxes.out/storage.out/thetas.out, so
    // fill in the same default ("%.8f") used by validate_and_fix_output_format().
    if (strlen(config->output_value_format) == 0) {
        snprintf(config->output_value_format, sizeof(config->output_value_format), "%%.8f");
    }
    return 0;
}

// This function concatenates the output filename to the output path
void build_full_output_path(const char* base_path, const char* filename, 
                              char* full_path, size_t full_path_size) {
      if (base_path && strlen(base_path) > 0) {
          snprintf(full_path, full_path_size, "%s%s", base_path, filename);
      } else {
          snprintf(full_path, full_path_size, "%s", filename);
      }
  }
  
// Note: output_path_name string must end in "/".  If not, add the "/".
//###############################
void validate_and_fix_output_path(char* path_str, size_t buffer_size) {
      if (!path_str || strlen(path_str) == 0) {
          snprintf(path_str, buffer_size, "./");
          return;
      }
      
      size_t len = strlen(path_str);
      if (path_str[len-1] != '/') {
          if (len < buffer_size - 1) {
              strcat(path_str, "/");
          }
      }
  }
  
// Validate and fix output format string iff necessary
// Returns 1 iff format was valid, 0 iff it was corrected to defaultt
//################################
int validate_and_fix_output_format(char* format_str, size_t buffer_size) {
    if (!format_str || strlen(format_str) == 0) {
        snprintf(format_str, buffer_size, "%%.8f");
        return 0;
    }
    
    // Check if format matches %.Nf or %.Ne pattern
    int precision;
    char format_type;
    
    // Try to parse the format string
    int matches = sscanf(format_str, "%%.%d%c", &precision, &format_type);
    
    if (matches == 2 && 
        (format_type == 'f' || format_type == 'e') && 
        precision >= 1 && precision <= 15) {
        // Valid format
        return 1;
    }
    
    // Invalid format - set to default
    snprintf(format_str, buffer_size, "%%.8f");
    return 0;
}

// For each of the initial storages read ensure that the correct number of values was read and applied
//###########################
int validate_giuh_arrays(CFE_CONFIG* config, const PARSER_ARRAY_COUNTS* counts) {

    // GIUH validation
    int expected_num = config->surface_routing_num_giuh_ordinates;
    int ordinates_read = counts->num_giuh_ordinates_read;
    int queue_read = counts->num_giuh_convolution_read;
    
    if (ordinates_read != expected_num) {
        fprintf(stderr, "ERROR: Expected %d GIUH ordinates, read: %d\n", 
                expected_num, ordinates_read);
        return -1;
    }
    
    if (queue_read != expected_num) {
        fprintf(stderr, "ERROR: Expected %d GIUH queue values, read: %d\n", 
                expected_num, queue_read);
        return -1;
    }
    
    // Validate GIUH ordinates sum to 1.0 and are non-negative
    double sum = 0.0;
    int is_negative = FALSE;
    for (int i = 0; i < expected_num; i++) {
        sum += config->surface_routing_giuh_ordinates[i];
        if(config->surface_routing_giuh_ordinates[i] < 0.0) is_negative = TRUE;
    }
    if (fabs(sum - 1.0) > 1e-6) {
        fprintf(stderr, "ERROR: GIUH ordinates sum to %.6f, must sum to 1.0\n", sum);
        return -1;
    }
    if(is_negative) {
        fprintf(stderr,"ERROR: negative GIUH ordinate encountered. All ordinates must be positive and sum to 1.0.\n");
        return -1;
    }
    
    // Check for negative GIUH queue values
    for (int i = 0; i < expected_num; i++) {
        if (config->surface_routing_init_giuh_convolution_queue_m[i] < 0.0) {
            fprintf(stderr, "ERROR: GIUH queue[%d] cannot be negative: %.6f\n", 
                    i, config->surface_routing_init_giuh_convolution_queue_m[i]);
            return -1;
        }
    }
    
    
    // Subsurface Nash validation (always expects 2)-----------------
    int subsurf_read = counts->num_subsurf_nash_storages_read;
    if (subsurf_read != 2) {
        fprintf(stderr, "ERROR: Expected 2 subsurface Nash storage values, read: %d\n", 
                subsurf_read);
        return -1;
    }
    
    // Check for negative subsurface Nash storage values
    for (int i = 0; i < 2; i++) {
        if (config->subsurface_routing_nash_cascade_init_storage_m[i] < 0.0) {
            fprintf(stderr, "ERROR: Subsurface Nash storage[%d] cannot be negative: %.6f\n", 
                    i, config->subsurface_routing_nash_cascade_init_storage_m[i]);
            return -1;
        }
    }
    
    return 0;
}


// Function to print the configuration (ffor testing) - Updated to include new fields
//###############
void print_config(const CFE_CONFIG* config) {
    printf("CFE Configuration:\n");
    printf("Version: %f\n", config->version);
    printf("Timestep: %.1f %s\n", config->timestep_h, config->timestep_units);
    printf("Catchment ID: %s\n", config->cat_id);
    printf("Latitude: %.6f %s\n", config->cat_latitude, config->cat_latitude_units);
    printf("Longitude: %.6f %s\n", config->cat_longitude, config->cat_longitude_units);
    printf("Elevation: %.2f %s\n", config->cat_elev, config->cat_elev_units);
    printf("Area: %.2f %s\n", config->cat_area_km2, config->cat_area_units);
    printf("Impervious fraction: %.3f %s\n", config->cat_impervious_fraction, config->cat_impervious_units);
    printf("Forcing file: %s\n", config->control_input_forcing_filename);
    printf("Total timesteps: %d\n", config->total_timesteps);
    printf("Verbosity: %d\n", config->verbosity);
    printf("Soil depth: %.2f %s\n", config->soil_depth_m, config->soil_depth_units);
    printf("Soil Clapp-Hornberger exponent: %.3f\n", config->soil_Clapp_Hornberger_exponent_b);
    printf("Soil sat hydraulic conductivity: %.4f %s\n", 
           config->soil_sat_hydraulic_conductivity_cm_per_h, 
           config->soil_sat_hydraulic_conductivity_units);
    printf("Soil discrete storage theta: %.2f, %.2f, %.2f, %.2f %s\n",
           config->soil_reservoir_init_discrete_storage_theta[0],
           config->soil_reservoir_init_discrete_storage_theta[1],
           config->soil_reservoir_init_discrete_storage_theta[2],
           config->soil_reservoir_init_discrete_storage_theta[3],
           config->soil_reservoir_init_discrete_storage_theta_units);
    printf("GW max storage: %.3f %s\n", config->gw_reservoir_max_storage_m, config->gw_reservoir_max_storage_units);
    printf("Rainfall partitioning scheme: %s\n", config->partitioning_scheme_name);
    printf("Surface routing scheme: %s\n", config->surface_routing_scheme_name);
    printf("Control discrete soil moisture: %s\n", 
           config->control_soil_simulate_discrete_soil_moisture_true_false ? "TRUE" : "FALSE");
    printf("Control use soil lookup table: %s\n", 
           config->control_soil_use_lookup_table_num_points ? "TRUE" : "FALSE");
    printf("Number of look-up table points to calculate: %d\n",config->control_soil_use_lookup_table_num_points);
    printf("\nOutput Configuration:\n");
    printf("Status warnings: %s\n", config->output_status_warnings_filename);
    printf("Internal fluxes: %s\n", config->output_internal_fluxes_filename);
    printf("Internal storages: %s\n", config->output_internal_storages_filename);
    printf("Volume balance: %s\n", config->output_volume_balance_filename);
    printf("Soil moisture theta: %s\n", config->output_soil_moisture_theta_filename);
    printf("Discharge: %s\n", config->output_discharge_filename);
    printf("Time format: %s\n", config->output_time_standard_format);
    printf("File delimiter: %s\n", config->output_file_delimiter);
}

// Add these functions to parser_helpers.c:

// Convert delimiter name to actual delimiter character
//##############################
const char* get_delimiter_string(const char* delimiter_name) {
    if (!delimiter_name) return ",";  // Default
    
    if (string_compare_ignore_case(delimiter_name, "comma") == 0) {
        return ",";
    } else if (string_compare_ignore_case(delimiter_name, "space") == 0) {
        return " ";
    } else if (string_compare_ignore_case(delimiter_name, "tab") == 0) {
        return "\t";
    } else if (string_compare_ignore_case(delimiter_name, "pipe") == 0) {
        return "|";
    } else {
        return ",";  // Default fallback
    }
}


// Validate time format
//######################
int validate_time_format(const char* time_format) {
    if (!time_format) return 0;
    
    if (string_compare_ignore_case(time_format, "timestep") == 0 ||
        string_compare_ignore_case(time_format, "datetime") == 0 ||
        string_compare_ignore_case(time_format, "datetime_dash") == 0 ||
        string_compare_ignore_case(time_format, "juliandate") == 0) {
        return 1;  // Valid
    }
    
    
    return 0;  // Invalid
}
