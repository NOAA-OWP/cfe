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

// These functions were developed inside cfe_main_driver.c, and moved here to get them out of that file.  They
// are all related to the CFE main() standalone (non-BMI) driver.  FLO 9/2025

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cfe_driver_utils.h"
#include "cfe_helpers.h"
#include "cfe.h"

#define TIME_STRING_LENGTH 64

//#####################
void print_usage(const char* prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  %s -c <config_filename> -f <forcing_filename> [OPTIONS]\n"
        "\n"
        "Required:\n"
        "  -c <file>     Configuration file\n"
        "\n"
        "Optional:\n"
        "  -f <file>     Forcing data file (overrides config)\n"
        "  -q <file>     Discharge per sq. m per timestep output file (overrides config)\n"
        "  -b <file>     Volume balance summary file (overrides config)\n"
        "  -x <file>     Internal fluxes output file (overrides config)\n"
        "  -s <file>     Internal storages output file (overrides config)\n"
        "  -t <file>     Soil moisture theta output file (overrides config, requires discrete soil moisture)\n"
        "  -v <level>    Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)\n"
        "  -dryrun       Run 120 time steps without forcing data (for testing).\n"
        "\n"
        "Examples:\n"
        "  Note: Config file extension is arbitrary (.cf3, .cfn, etc.); version is read from cfe_config_version field.\n"
        "  %s -c config.cfN -f forcings/forcing.dat\n"
        "  %s -c config.cfN -f forcings/forcing.dat -q output/q.out -b output/volbal.out\n"
        "  %s -c config.cfN -f forcings/forcing.dat -x fluxes.out -s stores.out -t thetas.out\n",
        prog, prog, prog, prog);
}

//###########################
int parse_command_line(int argc, char* argv[], cfe_cmdline_args_struct* args)
{
    // Initialize all fields
    memset(args, 0, sizeof(*args));
    args->command_line_verbosity = -1;  // -1 means not set
    args->run_without_forcing = FALSE;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            args->forcing_path = argv[++i];
        }
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            args->volbal_path = argv[++i];
        }
        else if (strcmp(argv[i], "-x") == 0 && i + 1 < argc) {
            args->fluxes_path = argv[++i];
        }
        else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            args->stores_path = argv[++i];
        }
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            args->thetas_path = argv[++i];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            args->qout_path = argv[++i];
        }
        else if (strcmp(argv[i], "-v") == 0 && i + 1 < argc) {
            args->command_line_verbosity = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            args->cfg_path = argv[++i];
        }
        else if (strcmp(argv[i], "-dryrun") == 0 || strcmp(argv[i], "--dryrun") == 0) {
            args->run_without_forcing = TRUE;
            args->forcing_path = NULL;  // override any forcing input
        }
        else {
            print_usage(argv[0]);
            return 1;
        }
    }

    if (!args->cfg_path) {
        fprintf(stderr, "RUNTIME ERROR: CFE config file required.\n");
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}

//#####################################
int check_command_line_overrides(cfe_options_struct* options, const cfe_cmdline_args_struct* args)
{
    // Apply command line verbosity override
    if (args->command_line_verbosity >= 0) {
        options->verbosity = args->command_line_verbosity;
    }
    
    // Handle forcing override
    if (args->run_without_forcing) {
        options->input_forcing_filename[0] = '\0'; // override config value
    }
    
    return 0;
}

//##########################
int open_output_files(const cfe_options_struct* options, const cfe_cmdline_args_struct* args,
                     FILE** q_fptr, FILE** Q_fptr, FILE** fluxes_fptr, FILE** storages_fptr, 
                     FILE** volbal_fptr, FILE** theta_fptr)
{
    // Initialize all file pointers to NULL
    *q_fptr = NULL;
    *Q_fptr = NULL;
    *fluxes_fptr = NULL;
    *storages_fptr = NULL;
    *volbal_fptr = NULL;
    *theta_fptr = NULL;

    char full_path[PATH_FILENAME_STRING_LENGTH];
    
    // Open files based on config settings first
    
    // Discharge (m per timestep) output
    //-------------------------------------
    if (options->output_discharge_filename[0] != '\0') {
        build_full_output_path(options->output_path_name, 
                               options->output_discharge_filename, 
                               full_path, sizeof(full_path));    
        *q_fptr = fopen(full_path, "w");
        if (!*q_fptr) {
            fprintf(stderr, "ERROR: Cannot open discharge output file from config: ");
            perror(full_path);  // Use full_path for accurate error reporting
            return 1;
        }
        if (options->verbosity > 0) printf("Writing discharge time series to: %s\n", full_path);  // Show full path so user knows where file is actually written
    }

    // Total discharge (m3/s) output
    //-------------------------------------
    if (options->output_total_discharge_m3_per_sec_filename[0] != '\0') {
        build_full_output_path(options->output_path_name, 
                               options->output_total_discharge_m3_per_sec_filename, 
                               full_path, sizeof(full_path));    
        *Q_fptr = fopen(full_path, "w");
        if (!*Q_fptr) {
            fprintf(stderr, "ERROR: Cannot open total discharge output file from config: ");
            perror(full_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Writing total discharge time series to: %s\n", full_path);
    }
    
    // Internal fluxes output
    //-------------------------------------
    if (options->output_internal_fluxes_filename[0] != '\0') {
        build_full_output_path(options->output_path_name, 
                               options->output_internal_fluxes_filename, 
                               full_path, sizeof(full_path));    
        *fluxes_fptr = fopen(full_path, "w");
        if (!*fluxes_fptr) {
            fprintf(stderr, "ERROR: Cannot open fluxes output file from config: ");
            perror(full_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Writing flux time series to: %s\n", full_path);
    }

    // Internal storages output
    //-------------------------------------
    if (options->output_internal_storages_filename[0] != '\0') {
        build_full_output_path(options->output_path_name, 
                               options->output_internal_storages_filename, 
                               full_path, sizeof(full_path));    
        *storages_fptr = fopen(full_path, "w");
        if (!*storages_fptr) {
            fprintf(stderr, "ERROR: Cannot open storages output file from config: ");
            perror(full_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Writing storage time series to: %s\n", full_path);
    }

    // Soil moisture disc theta time series
    //-------------------------------------
    if (options->output_soil_moisture_theta_filename[0] != '\0' && options->simulate_discrete_soil_moisture) {
        build_full_output_path(options->output_path_name, 
                               options->output_soil_moisture_theta_filename, 
                               full_path, sizeof(full_path));    
        *theta_fptr = fopen(full_path, "w");
        if (!*theta_fptr) {
            fprintf(stderr, "ERROR: Cannot open theta output file from config: ");
            perror(full_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Writing soil moisture theta time series to: %s\n", full_path);
    }
    
    // Output volume balance summary
    //-------------------------------------
    if (options->output_volume_balance_filename[0] != '\0') {
        build_full_output_path(options->output_path_name, 
                               options->output_volume_balance_filename, 
                               full_path, sizeof(full_path));    
        *volbal_fptr = fopen(full_path, "w");
        if (!*volbal_fptr) {
            fprintf(stderr, "ERROR: Cannot open volume balance output file from config: ");
            perror(full_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Writing volume balance summary to: %s\n", full_path);
    }

    // Apply command line overrides (close config files and open CLO files)
    
    if (args->qout_path) {
        if (*q_fptr) fclose(*q_fptr);
        *q_fptr = fopen(args->qout_path, "w");
        if (!*q_fptr) {
            fprintf(stderr, "ERROR: Cannot open command line discharge output file: ");
            perror(args->qout_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Command line override: writing discharge to: %s\n", args->qout_path);
    }

    if (args->volbal_path) {
        if (*volbal_fptr) fclose(*volbal_fptr);
        *volbal_fptr = fopen(args->volbal_path, "w");
        if (!*volbal_fptr) {
            fprintf(stderr, "ERROR: Cannot open command line volume balance output file: ");
            perror(args->volbal_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Command line override: writing volume balance to: %s\n", args->volbal_path);
    }

    if (args->fluxes_path) {
        if (*fluxes_fptr) fclose(*fluxes_fptr);
        *fluxes_fptr = fopen(args->fluxes_path, "w");
        if (!*fluxes_fptr) {
            fprintf(stderr, "ERROR: Cannot open command line fluxes output file: ");
            perror(args->fluxes_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Command line override: writing fluxes to: %s\n", args->fluxes_path);
    }

    if (args->stores_path) {
        if (*storages_fptr) fclose(*storages_fptr);
        *storages_fptr = fopen(args->stores_path, "w");
        if (!*storages_fptr) {
            fprintf(stderr, "ERROR: Cannot open command line storages output file: ");
            perror(args->stores_path);
            return 1;
        }
        if (options->verbosity > 0) printf("Command line override: writing storages to: %s\n", args->stores_path);
    }

    if (args->thetas_path) {
        if (!options->simulate_discrete_soil_moisture) {
            fprintf(stderr, "WARNING: -t option ignored because discrete soil moisture is not enabled\n");
        } else {
            if (*theta_fptr) fclose(*theta_fptr);
            *theta_fptr = fopen(args->thetas_path, "w");
            if (!*theta_fptr) {
                fprintf(stderr, "ERROR: Cannot open command line theta output file: ");
                perror(args->thetas_path);
                return 1;
            }
            if (options->verbosity > 0) printf("Command line override: writing soil moisture thetas to: %s\n", args->thetas_path);
        }
    }

    return 0;
}

//############################
void close_output_files(FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr, FILE* storages_fptr, 
                        FILE* volbal_fptr, FILE* theta_fptr, FILE* warnings_fptr, 
                        FILE* forcing_fptr)
{
    if (q_fptr)        fclose(q_fptr);
    if (Q_fptr)        fclose(Q_fptr);
    if (fluxes_fptr)   fclose(fluxes_fptr);
    if (storages_fptr) fclose(storages_fptr);
    if (theta_fptr)    fclose(theta_fptr);   
    if (volbal_fptr)   fclose(volbal_fptr);
    if (warnings_fptr) fclose(warnings_fptr);
    if (forcing_fptr)  fclose(forcing_fptr);
}

static int match_dev_column(const char *tok, aorc_cols_t *cols, int idx)
{
    if (strcmp(tok, "time") == 0)                { cols->time_idx = idx; return 1; }
    if (strcmp(tok, "APCP_surface") == 0)        { cols->apcp_idx = idx; return 1; }
    if (strcmp(tok, "precip_rate") == 0)         { cols->precip_rate_idx = idx; cols->precip_rate_is_kg_m2_s1 = 0; return 1; }
    if (strcmp(tok, "DLWRF_surface") == 0)       { cols->dlwrf_idx = idx; return 1; }
    if (strcmp(tok, "DSWRF_surface") == 0)       { cols->dswrf_idx = idx; return 1; }
    if (strcmp(tok, "PRES_surface") == 0)        { cols->pres_idx  = idx; return 1; }
    if (strcmp(tok, "SPFH_2maboveground") == 0)  { cols->spfh_idx  = idx; return 1; }
    if (strcmp(tok, "TMP_2maboveground") == 0)   { cols->tmp_idx   = idx; return 1; }
    if (strcmp(tok, "UGRD_10maboveground") == 0) { cols->ugrd_idx  = idx; return 1; }
    if (strcmp(tok, "VGRD_10maboveground") == 0) { cols->vgrd_idx  = idx; return 1; }
    return 0;
}

static int match_ngen_column(const char *tok, aorc_cols_t *cols, int idx)
{
    if (strcmp(tok, "Time") == 0)    { cols->time_idx = idx; return 1; }
    if (strcmp(tok, "RAINRATE") == 0){ cols->precip_rate_idx = idx; cols->precip_rate_is_kg_m2_s1 = 1; return 1; }
    if (strcmp(tok, "LWDOWN") == 0)  { cols->dlwrf_idx = idx; return 1; }
    if (strcmp(tok, "SWDOWN") == 0)  { cols->dswrf_idx = idx; return 1; }
    if (strcmp(tok, "PSFC") == 0)    { cols->pres_idx  = idx; return 1; }
    if (strcmp(tok, "Q2D") == 0)     { cols->spfh_idx  = idx; return 1; }
    if (strcmp(tok, "T2D") == 0)     { cols->tmp_idx   = idx; return 1; }
    if (strcmp(tok, "U2D") == 0)     { cols->ugrd_idx  = idx; return 1; }
    if (strcmp(tok, "V2D") == 0)     { cols->vgrd_idx  = idx; return 1; }
    return 0;
}

//##########################
int open_forcing_file(const char* path, FILE** fptr, aorc_cols_t* cols)
{
    if (path == NULL || path[0] == '\0' || fptr == NULL || cols == NULL) return -1;
    *fptr = fopen(path, "r");
    if (*fptr == NULL) return -1;

    // Read header line
    char header[2048];
    if (fgets(header, sizeof(header), *fptr) == NULL) {
        fclose(*fptr);
        *fptr = NULL;
        return -1;
    }

    // Initialize all column indices to -1 (not found)
    cols->time_idx = -1;
    cols->apcp_idx = -1;
    cols->precip_rate_idx = -1;
    cols->dlwrf_idx = -1;
    cols->dswrf_idx = -1;
    cols->pres_idx = -1;
    cols->spfh_idx = -1;
    cols->tmp_idx = -1;
    cols->ugrd_idx = -1;
    cols->vgrd_idx = -1;
    cols->precip_rate_is_kg_m2_s1 = 0;
    cols->detected_format = FORCING_FMT_UNKNOWN;

    int idx = 0;
    int dev_matches = 0, ngen_matches = 0;
    for (char* tok = strtok(header, ",\r\n"); tok != NULL; tok = strtok(NULL, ",\r\n"), idx++) {
        dev_matches  += match_dev_column(tok, cols, idx);
        ngen_matches += match_ngen_column(tok, cols, idx);
    }

    if (ngen_matches > 0 && dev_matches == 0)
        cols->detected_format = FORCING_FMT_NGEN;
    else if (dev_matches > 0 && ngen_matches == 0)
        cols->detected_format = FORCING_FMT_DEV;
    else
        cols->detected_format = FORCING_FMT_UNKNOWN;

    // Check for required columns
    if (cols->time_idx < 0) {
        fprintf(stderr, "ERROR: AORC header missing 'time' or 'Time' column\n");
        fclose(*fptr);
        *fptr = NULL;
        return -1;
    }
    if (cols->apcp_idx < 0 && cols->precip_rate_idx < 0) {
        fprintf(stderr, "ERROR: AORC header missing precipitation column "
                "(need 'APCP_surface', 'precip_rate', or 'RAINRATE')\n");
        fclose(*fptr);
        *fptr = NULL;
        return -1;
    }
    
    return 0;
}

//###############################
int read_next_forcing_aorc(FILE* f, const aorc_cols_t* cols, int dt_seconds,
                           cfe_forcing_struct* forcing, char* time_str)
{
    if (f == NULL || cols == NULL || forcing == NULL) return -1;
    
    char line[2048];
    if (fgets(line, sizeof(line), f) == NULL) return 0;  // EOF
    
    // Initialize all fields to 0
    memset(forcing, 0, sizeof(*forcing));
    
    // Split the line into tokens
    int idx = 0;
    int have_apcp = 0;
    int have_rate = 0;
    
    for (char* tok = strtok(line, ",\r\n"); tok != NULL; tok = strtok(NULL, ",\r\n"), idx++) {
        if (idx == cols->time_idx && time_str != NULL) {
            // Copy the time string if requested
            strncpy(time_str, tok, 63);
            time_str[63] = '\0';
            // Clean up time string (remove quotes/whitespace)
            char* start = time_str;
            while (*start && (*start == '"' || *start == ' ')) start++;
            if (start != time_str) {
                memmove(time_str, start, strlen(start) + 1);
            }
            char* end = time_str + strlen(time_str) - 1;
            while (end > time_str && (*end == '"' || *end == ' ' || *end == '\n' || *end == '\r')) {
                *end = '\0';
                end--;
            }
        }
        else if (idx == cols->apcp_idx && *tok != '\0') {
            forcing->APCP_surface = atof(tok);
            have_apcp = 1;
        }
        else if (idx == cols->precip_rate_idx && *tok != '\0') {
            forcing->precip_rate = atof(tok);
            if (cols->precip_rate_is_kg_m2_s1) {
                forcing->precip_rate /= 1000.0;
            }
            have_rate = 1;
        }
        else if (idx == cols->dlwrf_idx && *tok != '\0') {
            forcing->DLWRF_surface = atof(tok);
        }
        else if (idx == cols->dswrf_idx && *tok != '\0') {
            forcing->DSWRF_surface = atof(tok);
        }
        else if (idx == cols->pres_idx && *tok != '\0') {
            forcing->PRES_surface = atof(tok);
        }
        else if (idx == cols->spfh_idx && *tok != '\0') {
            forcing->SPFH_2maboveground = atof(tok);
        }
        else if (idx == cols->tmp_idx && *tok != '\0') {
            forcing->TMP_2maboveground = atof(tok);
        }
        else if (idx == cols->ugrd_idx && *tok != '\0') {
            forcing->UGRD_10maboveground = atof(tok);
        }
        else if (idx == cols->vgrd_idx && *tok != '\0') {
            forcing->VGRD_10maboveground = atof(tok);
        }
    }
    
    // Calculate rainfall_depth_m (convert mm to m)
    if (have_apcp) {
        forcing->rainfall_depth_m = forcing->APCP_surface / 1000.0;
    } else if (have_rate) {
        forcing->rainfall_depth_m = forcing->precip_rate * (double)dt_seconds;
    } else {
        forcing->rainfall_depth_m = 0.0;
    }

    forcing->et_potential_m = 0.0;  // calculated later if Priestley-Taylor is enabled

    forcing->day_of_year = 1;
    if (time_str != NULL && time_str[0] != '\0') {
        aorc_forcing_time_struct parsed_time = {0};
        if (parse_time_string(time_str, &parsed_time) == 0) {
            forcing->day_of_year = calculate_day_of_year(parsed_time.year,
                                                          parsed_time.month,
                                                          parsed_time.day);
        }
    }

    return 1;
}

//#########################
int parse_time_string(const char* time_str, aorc_forcing_time_struct* forcing_time)
{
    if (!time_str || !forcing_time) return -1;

    // Accept both "2015/12/01 00:00:00" (CFE dev format) and
    // "2012-10-01 00:00:00" (NextGen/AORC-standard format).
    char normalized[TIME_STRING_LENGTH];
    strncpy(normalized, time_str, sizeof(normalized) - 1);
    normalized[sizeof(normalized) - 1] = '\0';
    for (char* p = normalized; *p; p++) {
        if (*p == '-') *p = '/';
    }

    int parsed = sscanf(normalized, "%d/%d/%d %d:%d:%d",
                       &forcing_time->year,
                       &forcing_time->month,
                       &forcing_time->day,
                       &forcing_time->hour,
                       &forcing_time->minute,
                       &forcing_time->second);
    
    if (parsed != 6) {
        printf("Warning: Could not parse time string '%s', using defaults\n", time_str);
        // Set defaults
        forcing_time->year = 2015;
        forcing_time->month = 12;
        forcing_time->day = 1;
        forcing_time->hour = 0;
        forcing_time->minute = 0;
        forcing_time->second = 0;
        return -1;
    }
    
    return 0;
}

//#############################
int calculate_time_delta(const char* time1_str, const char* time2_str)
{
    aorc_forcing_time_struct time1 = {0}, time2 = {0};
    
    if (parse_time_string(time1_str, &time1) != 0 || 
        parse_time_string(time2_str, &time2) != 0) {
        return 3600;  // Default to 1 hour
    }
    
    // aorc data are hourly.   verify.
    int delta_hours = (time2.hour - time1.hour);
    if (delta_hours <= 0) delta_hours += 24;  // Handle day rollover
    
    if(delta_hours != 1) {
        fprintf(stdout,"WARNING: AORC delta t differs from 1 hour in calculate_time_delta().\n");
    }
    return delta_hours * 3600;  // Convert to seconds
}

// This function determines the date/time of the first data, number of data lines, and delta_t of the forcings data
//############################
int analyze_forcing_file(const char* forcing_filename, aorc_forcing_time_struct* forcing_time, int verbosity)
{
    if (!forcing_filename || !forcing_time) return -1;
    
    FILE* forcing_fptr = NULL;
    aorc_cols_t cols;
    
    // Use your existing function to open and parse header
    if (open_forcing_file(forcing_filename, &forcing_fptr, &cols) != 0) {
        return -1;
    }
    
    // Create temporary forcing struct for reading
    cfe_forcing_struct temp_forcing;
    
    // Read first data line to get start time
    char first_time_str[TIME_STRING_LENGTH];
    if (read_next_forcing_aorc(forcing_fptr, &cols, 3600, &temp_forcing, first_time_str) <= 0) {
        fclose(forcing_fptr);
        return -1;
    }
    
    // Parse the start time
    if (parse_time_string(first_time_str, forcing_time) != 0) {
        printf("Warning: Using default start time\n");
    }
    
    // Read second line to determine delta_t
    char second_time_str[TIME_STRING_LENGTH];
    if (read_next_forcing_aorc(forcing_fptr, &cols, 3600, &temp_forcing, second_time_str) > 0) {
        forcing_time->delta_t_seconds = calculate_time_delta(first_time_str, second_time_str);
    } else {
        forcing_time->delta_t_seconds = 3600;  // Default to hourly
    }
    
    // Count remaining valid lines
    int valid_lines = 2;  // We already read two
    while (read_next_forcing_aorc(forcing_fptr, &cols, 3600, &temp_forcing, NULL) > 0) {
        valid_lines++;
    }
    
    forcing_time->num_valid_lines = valid_lines;
    
    // Store the time format for output use
    strncpy(forcing_time->time_format, "YYYY/MM/DD HH:MM:SS", sizeof(forcing_time->time_format) - 1);
    forcing_time->time_format[sizeof(forcing_time->time_format) - 1] = '\0';
    
    fclose(forcing_fptr);
    if(verbosity > 1)
           printf("Analyzed forcing file: start=%04d/%02d/%02d %02d:%02d:%02d, dt=%d s, lines=%d\n",
           forcing_time->year, forcing_time->month, forcing_time->day,
           forcing_time->hour, forcing_time->minute, forcing_time->second,
           forcing_time->delta_t_seconds, forcing_time->num_valid_lines);
    
    return 0;
}

// Format timestamp based on the requested format
//####################
void format_timestamp(char* timestamp_str, size_t str_size, int timestep,
                     const aorc_forcing_time_struct* forcing_time,
                     const char* time_format) {
    
    if (string_compare_ignore_case(time_format, "timestep") == 0) {
        snprintf(timestamp_str, str_size, "%d", timestep);
    }
    else if (string_compare_ignore_case(time_format, "datetime") == 0 ||
             string_compare_ignore_case(time_format, "datetime_dash") == 0) {
        // Build a struct tm from the forcing_time start (assumed UTC)
        struct tm t = {0};
        t.tm_year = forcing_time->year - 1900;  // struct tm years since 1900
        t.tm_mon  = forcing_time->month - 1;    // struct tm months 0-11
        t.tm_mday = forcing_time->day;
        t.tm_hour = forcing_time->hour;
        t.tm_min  = forcing_time->minute;
        t.tm_sec  = forcing_time->second;

        // Convert to time_t as UTC, not local time
    #if defined(_WIN32) || defined(_WIN64)
        time_t start = _mkgmtime(&t);   // Windows
    #else
        time_t start = timegm(&t);      // POSIX (Linux/macOS)
    #endif

        // Advance by timestep * delta_t
        time_t now = start + (time_t)timestep * forcing_time->delta_t_seconds;

        // Convert back to UTC calendar time
    #if defined(_WIN32) || defined(_WIN64)
        struct tm curr;
        gmtime_s(&curr, &now);
    #else
        struct tm curr;
        gmtime_r(&now, &curr);
    #endif

        // "datetime" -> YYYY/MM/DD HH:MM:SS (customary CFE format)
        // "datetime_dash" -> YYYY-MM-DD HH:MM:SS (ISO-8601-like)
        if (string_compare_ignore_case(time_format, "datetime_dash") == 0) {
            snprintf(timestamp_str, str_size, "%04d-%02d-%02d %02d:%02d:%02d",
                    curr.tm_year + 1900, curr.tm_mon + 1, curr.tm_mday,
                    curr.tm_hour, curr.tm_min, curr.tm_sec);
        } else {
            snprintf(timestamp_str, str_size, "%04d/%02d/%02d %02d:%02d:%02d",
                    curr.tm_year + 1900, curr.tm_mon + 1, curr.tm_mday,
                    curr.tm_hour, curr.tm_min, curr.tm_sec);
        }
    }
    else if (string_compare_ignore_case(time_format, "juliandate") == 0) {
        // Calculate Julian date for current timestep
        int total_seconds = timestep * forcing_time->delta_t_seconds;
        int total_hours = total_seconds / 3600;
        
        int current_hour = forcing_time->hour + total_hours;
        int current_day = forcing_time->day;
        int current_month = forcing_time->month;
        int current_year = forcing_time->year;
        
        // Handle day overflow (simplified)
        while (current_hour >= 24) {
            current_hour -= 24;
            current_day++;
        }
        
        double julian_date = calculate_julian_date(current_year, current_month, current_day,
                                                  current_hour, forcing_time->minute, 
                                                  forcing_time->second);
        snprintf(timestamp_str, str_size, "%.6f", julian_date);
    }
    else {
        // Default to timestep
        snprintf(timestamp_str, str_size, "%d", timestep);
    }
}

// Write headers for output files
//##############################
void write_output_headers(const cfe_options_struct* options,
                         const char* time_format,
                         const char* delimiter,
                         FILE* q_fptr,
                         FILE* Q_fptr,
                         FILE* fluxes_fptr,
                         FILE* storages_fptr,
                         FILE* theta_fptr) {
    
    if (q_fptr) {
        fprintf(q_fptr, "# CFE Discharge Output\n");
        fprintf(q_fptr, "#%s%sdischarge_m_per_timestep\n", time_format, delimiter);
    }

    if (Q_fptr) {
        fprintf(Q_fptr, "# CFE Total Discharge Output\n");
        fprintf(Q_fptr, "#%s%sdischarge_m3_per_s\n", time_format, delimiter);
    }

    if (fluxes_fptr) {
        fprintf(fluxes_fptr, "# CFE Internal Fluxes Output\n");
        fprintf(fluxes_fptr,
                "# %s%sinput_rainfall_m%simpervious_runoff_m%spervious_runoff_m"
                "%ssurface_runoff_generated_m%ssurface_routed_m"
                "%slateral_flow_generated_m%slateral_flow_routed_m"
                "%spotential_transpiration_m%sbare_soil_evaporation_m"
                "%sactual_transpiration_m%saet_m%spercolation_to_gw_m"
                "%sbaseflow_m%stotal_outflow_m\n",
                time_format,
                delimiter, delimiter, delimiter,
                delimiter, delimiter,
                delimiter, delimiter,
                delimiter, delimiter,
                delimiter, delimiter, delimiter,
                delimiter, delimiter);
    }
    
    if (theta_fptr) {
        fprintf(theta_fptr, "# CFE Soil Moisture Theta Output\n");
        fprintf(theta_fptr, "#%s%stheta1%stheta2%stheta3%stheta4\n", 
                time_format, delimiter, delimiter, delimiter, delimiter);
    }
    
    if (storages_fptr) {
        fprintf(storages_fptr, "# CFE Internal Storages Output\n");
        if(options->simulate_discrete_soil_moisture == FALSE) {
            fprintf(storages_fptr, "#%s%sstorage_in_soil_reservoir_m%sstorage_in_gw_reservoir_m"
                    "%ssurface_routing_storage_m%ssubsurface_lateral_storage_m\n",
                    time_format, delimiter, delimiter, delimiter, delimiter);
        } else {
            fprintf(storages_fptr, "#%s", time_format);
            for (int i = 0; i < NDISC; i++) {
                fprintf(storages_fptr, "%sstorage_in_disc%d_m", delimiter, i + 1);
            }
            fprintf(storages_fptr, "%stotal_storage_in_soil_m%sstorage_in_gw_reservoir_m"
                    "%ssurface_routing_storage_m%ssubsurface_lateral_storage_m\n",
                    delimiter, delimiter, delimiter, delimiter);
        }
    }
}

// Main output writing function
//###########################
void write_all_outputs(int timestep,
                       const cfe_options_struct* options,
                       const cfe_parameters_struct* params,
                       const cfe_state_struct* state,
                       const cfe_outputs_struct* outputs,
                       const cfe_forcing_struct* forcing,
                       const aorc_forcing_time_struct* forcing_time,
                       const cfe_volbal_struct* volbal,
                       FILE* q_fptr,
                       FILE* Q_fptr,
                       FILE* fluxes_fptr,
                       FILE* storages_fptr,
                       FILE* volbal_fptr,
                       FILE* theta_fptr) {
    
    // Get the actual delimiter string
    const char* delimiter = get_delimiter_string(options->output_file_delimiter);
    
    // Format timestamp based on config option
    char timestamp_str[TIME_STRING_LENGTH];
    format_timestamp(timestamp_str, sizeof(timestamp_str), timestep, forcing_time, options->output_time_standard_format);


    // Write discharge output
    if (q_fptr) {
        fprintf(q_fptr, "%s%s", timestamp_str, delimiter);
        fprintf(q_fptr, options->output_value_format, outputs->qout_m);
        fprintf(q_fptr, "\n");
    }

    if (Q_fptr) {
        // Convert discharge from m/timestep to m�/s
        double discharge_m3_per_sec = outputs->qout_m *  params->catchment_area_km2 * 1.0e+06/ options->time_step_seconds;
        fprintf(Q_fptr, "%s%s",timestamp_str, delimiter);
        fprintf(Q_fptr, options->output_value_format, discharge_m3_per_sec);
        fprintf(Q_fptr,"\n");
    }
    
    // Write fluxes output
    if (fluxes_fptr) {
        fprintf(fluxes_fptr, "%s%s", timestamp_str, delimiter);
        fprintf(fluxes_fptr, options->output_value_format, forcing->rainfall_depth_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->impervious_runoff_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->pervious_runoff_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->surface_runoff_generated_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->surface_routed_to_outlet_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->lateral_flow_generated_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->lateral_flow_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->potential_et_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->bare_soil_evaporation_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format,
                outputs->actual_et_m - outputs->bare_soil_evaporation_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->actual_et_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->soil_to_gw_percolation_flux_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->baseflow_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, options->output_value_format, outputs->total_outflow_m);
        fprintf(fluxes_fptr, "\n");
    }

    // Write storages output
    if (storages_fptr) {
        double surface_routing_storage_m = 0.0;
        for (int i = 0; i < params->giuh_num_ordinates; i++) {
            surface_routing_storage_m += state->giuh_queue_m[i];
        }

        double subsurface_lateral_storage_m = 0.0;
        for (int i = 0; i < params->nash_subsurface_N; i++) {
            subsurface_lateral_storage_m += state->nash_subsurface_storage_m[i];
        }

        fprintf(storages_fptr, "%s%s", timestamp_str, delimiter);
        if(options->simulate_discrete_soil_moisture == FALSE) {
            fprintf(storages_fptr, options->output_value_format, state->soil_storage_m);
        } else {
            for (int i = 0; i < NDISC; i++) {
                double storage_in_disc_m = state->soil_state_out.theta_out[i] *
                                           state->soil_geometry.dz_m[i];
                fprintf(storages_fptr, options->output_value_format, storage_in_disc_m);
                fprintf(storages_fptr, "%s", delimiter);
            }
            fprintf(storages_fptr, options->output_value_format, state->soil_state_out.total_storage_m);
        }
        fprintf(storages_fptr, "%s", delimiter);
        fprintf(storages_fptr, options->output_value_format, state->gw_storage_m);
        fprintf(storages_fptr, "%s", delimiter);
        fprintf(storages_fptr, options->output_value_format, surface_routing_storage_m);
        fprintf(storages_fptr, "%s", delimiter);
        fprintf(storages_fptr, options->output_value_format, subsurface_lateral_storage_m);
        fprintf(storages_fptr, "\n");
    }

    // Write theta output (iff discrete soil moisture is enabled and file is open)
    if (theta_fptr && options->simulate_discrete_soil_moisture) {
        fprintf(theta_fptr, "%s%s", timestamp_str, delimiter);
        fprintf(theta_fptr, options->output_value_format, state->soil_state_out.theta_out[0]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, options->output_value_format, state->soil_state_out.theta_out[1]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, options->output_value_format, state->soil_state_out.theta_out[2]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, options->output_value_format, state->soil_state_out.theta_out[3]);
        fprintf(theta_fptr, "\n");
    }
}

// Initialize the volume balance terms at t=0.
//################################
void cfe_initialize_volume_balance(const cfe_parameters_struct* params,
                                   const cfe_options_struct* options, 
                                   const cfe_state_struct* state,
                                   cfe_volbal_struct* volbal)
{
    // Add up initial subsurface Nash storage
    double volstart_subsurface = 0.0;
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        volstart_subsurface += state->nash_subsurface_storage_m[i];
    }

    double volstart_surface = 0.0;

    int fred_debug = FALSE;
    if(fred_debug) {
 
        if (options->surface_routing_scheme == SURF_ROUTE_GIUH) {
            if (options->verbosity > 0) {
                printf("DEBUG: Volume balance using GIUH - giuh_num_ordinates = %d\n", params->giuh_num_ordinates);
            }
            for (int i = 0; i < params->giuh_num_ordinates; i++) {
                if (options->verbosity > 0) {
                    printf("DEBUG: Adding GIUH queue[%d] = %.6f\n", i, state->giuh_queue_m[i]);
                }
                volstart_surface += state->giuh_queue_m[i];
            }
            if (options->verbosity > 0) {
                printf("DEBUG: Total GIUH volstart_surface = %.6f\n", volstart_surface);
            }
        }
    }
    else {
        for (int i = 0; i < params->giuh_num_ordinates; i++) {
            volstart_surface += state->giuh_queue_m[i];
        }
    }

    /** These are for Global volume balance */
    if (options->simulate_discrete_soil_moisture) {
        volbal->volstart_soil = state->soil_state_in.total_storage_m;
    } else {
        volbal->volstart_soil = state->soil_storage_m;
    }

    volbal->volstart_gw = state->gw_storage_m;
    volbal->volstart_surface = volstart_surface;      // ADD THIS
    volbal->volstart_subsurface = volstart_subsurface; // ADD THIS  
    volbal->volstart = volbal->volstart_soil + volbal->volstart_gw + volstart_surface + volstart_subsurface;
    
    /* these are for storage compartment volume balance accounting (yes, redundant, but ... */               
    volbal->vol_soil_start = state->soil_storage_m;
    volbal->vol_in_gw_start = state->gw_storage_m;
}

//######################################
void write_volume_balance_summary(FILE* output_fptr,
                                  const cfe_options_struct* options,
                                  const cfe_volbal_struct* volbal,
                                  const cfe_state_struct* final_state,
                                  const cfe_parameters_struct* params,
                                  double *soil_dz)
{
    if (output_fptr == NULL) return;
    
    // Calculate final volumes
    double vol_soil_end = final_state->soil_storage_m;
    double vol_gw_end = final_state->gw_storage_m;
    double vol_nash_subsurface_end = 0.0;
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        vol_nash_subsurface_end += final_state->nash_subsurface_storage_m[i];
    }
    
    double vol_surface_end = 0.0;
    if (options->surface_routing_scheme == SURF_ROUTE_GIUH) {
        for (int i = 0; i < params->giuh_num_ordinates; i++) {
            vol_surface_end += final_state->giuh_queue_m[i];
        }
    } else {  // using Nash cascade 
        for (int i = 0; i < MAX_NUM_SURFACE_NASH_CASCADE; i++)
            vol_surface_end += final_state->nash_surface_storage_m[i];
    }

    double total_AET_vol  = volbal->vol_et_from_rain + volbal->vol_et_from_soil;
 
    double volstart = volbal->volstart_soil + volbal->volstart_gw + volbal->volstart_surface + volbal->volstart_subsurface;
   
    // Configuration-aware headers
    const char* partition_name = (options->liquid_partitioning_scheme == PARTITION_SCHAAKE) ? "Schaake" : "Xinanjiang";
    const char* surface_name = (options->surface_routing_scheme == SURF_ROUTE_GIUH) ? "GIUH" : "Nash cascade";
        
    
    // GLOBAL VOLUME BALANCE
    double aet_total  = volbal->vol_et_from_rain
                      + volbal->vol_et_from_soil;

    double qout_total = volbal->vol_out_surface
                      + volbal->vol_out_subsurf_nash         
                      + volbal->vol_from_gw;

    double volout_total = qout_total + aet_total;

    // PRECIPITATION VOLUME BALANCE  
    double direct_residual = volbal->volin - volbal->vol_runoff - volbal->vol_infilt - volbal->vol_et_from_rain;
    
    fprintf(output_fptr, "\n****************** PRECIPITATION VOLUME BALANCE *****************\n");
    fprintf(output_fptr, " Volume input                       = %8.4lf m\n", volbal->volin);
    fprintf(output_fptr, " Surface runoff generated           = %8.4lf m\n", volbal->vol_runoff);
    fprintf(output_fptr, "   Impervious runoff                = %8.4lf m\n", volbal->vol_impervious_runoff);
    fprintf(output_fptr, "   Pervious runoff                  = %8.4lf m\n", volbal->vol_pervious_runoff);
    fprintf(output_fptr, " Added to soil moisture             = %8.4lf m\n", volbal->vol_infilt);
    fprintf(output_fptr, " Volume of ET from rain             = %8.4lf m\n", volbal->vol_et_from_rain);
    fprintf(output_fptr, " Precip residual                    = %6.4e m\n", direct_residual);
    if (fabs(direct_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: DIRECT RUNOFF PARTITIONING VOLUME BALANCE CHECK FAILED\n");
    }

    // SURFACE VOLUME BALANCE


    double surface_residual = volbal->volstart_surface + volbal->vol_runoff -
                              volbal->vol_out_surface - vol_surface_end;

    int use_sci_notation = FALSE;

    fprintf(output_fptr, "\n************ SURFACE ROUTING VOLUME BALANCE *********************\n");
    if(!use_sci_notation) {
      fprintf(output_fptr, " Initial surface routing storage    = %8.4lf m\n", volbal->volstart_surface);
      fprintf(output_fptr, " Runoff into surface routing        = %8.4lf m\n", volbal->vol_runoff);
      fprintf(output_fptr, " Outflow from surface routing       = %8.4lf m\n", volbal->vol_out_surface);
      fprintf(output_fptr, " Final surface routing storage      = %8.4lf m\n", vol_surface_end);
      fprintf(output_fptr, " Surface residual                   = %6.4e m\n", surface_residual);
    } else {
      fprintf(output_fptr, " Initial surface routing storage    = %e m\n", volbal->volstart_surface);
      fprintf(output_fptr, " Runoff into surface routing        = %e m\n", volbal->vol_runoff);
      fprintf(output_fptr, " Outflow from surface routing       = %e m\n", volbal->vol_out_surface);
      fprintf(output_fptr, " Final surface routing storage      = %e m\n", vol_surface_end);
      fprintf(output_fptr, " Surface residual                   = %6.4e m\n", surface_residual);
    }
    if (fabs(surface_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: SURFACE ROUTING VOLUME BALANCE CHECK FAILED\n");
    }

    if (options->simulate_discrete_soil_moisture) {
        fprintf(output_fptr, "\n************** DISCRETE SOIL BALANCE MODEL VOLUME BALANCE **************\n");
        vol_soil_end = 0.0;
        for(int i = 0; i < NDISC; i++) {
            vol_soil_end += final_state->soil_discrete_storage_theta[i] * soil_dz[i]; 
        }
    } else {
        fprintf(output_fptr, "\n********* SOIL WATER CONCEPTUAL RESERVOIR VOLUME BALANCE ********\n");
        vol_soil_end = final_state->soil_storage_m;
    }

    double soil_residual = volbal->vol_soil_start + volbal->vol_infilt -
                           volbal->vol_soil_to_lat_flow - volbal->vol_to_gw -
                           volbal->vol_et_from_soil - vol_soil_end;

    use_sci_notation = FALSE;
    if(use_sci_notation == TRUE) {
        fprintf(output_fptr, " Initial soil vol.                  = %e m\n", volbal->vol_soil_start);
        fprintf(output_fptr, " Infiltration into soil             = %e m\n", volbal->vol_infilt);
        fprintf(output_fptr, " From soil to lat. flow             = %e m\n", volbal->vol_soil_to_lat_flow);
        fprintf(output_fptr, " From from soil to GW               = %e m\n", volbal->vol_soil_to_gw);
        fprintf(output_fptr, " Vegetated-area AET                 = %e m\n", volbal->vol_forest_aet);
        fprintf(output_fptr, " Bare-soil evaporation              = %e m\n", volbal->vol_bare_soil_evaporation);
        fprintf(output_fptr, " Total ET from soil                 = %e m\n", volbal->vol_et_from_soil);
        fprintf(output_fptr, " Final soil vol.                    = %e m\n", vol_soil_end);
        fprintf(output_fptr, " Soil residual                      = %6.4e m\n", soil_residual);
    } else {
        fprintf(output_fptr, " Initial soil vol.                  = %8.4lf m\n", volbal->vol_soil_start);
        fprintf(output_fptr, " Infiltration into soil             = %8.4lf m\n", volbal->vol_infilt);
        fprintf(output_fptr, " From soil to lat. flow             = %8.4lf m\n", volbal->vol_soil_to_lat_flow);
        fprintf(output_fptr, " From from soil to GW               = %8.4lf m\n", volbal->vol_soil_to_gw);
        fprintf(output_fptr, " Vegetated-area AET                 = %8.4lf m\n", volbal->vol_forest_aet);
        fprintf(output_fptr, " Bare-soil evaporation              = %8.4lf m\n", volbal->vol_bare_soil_evaporation);
        fprintf(output_fptr, " Total ET from soil                 = %8.4lf m\n", volbal->vol_et_from_soil);
        fprintf(output_fptr, " Final soil vol.                    = %8.4lf m\n", vol_soil_end);
        fprintf(output_fptr, " Soil residual                      = %6.4e m\n", soil_residual);
    }
    if (fabs(soil_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: SOIL CONCEPTUAL RESERVOIR VOLUME BALANCE CHECK FAILED\n");
    }

    // NASH CASCADE VOLUME BALANCE
    double nash_residual = volbal->vol_in_subsurf_nash - volbal->vol_out_subsurf_nash - vol_nash_subsurface_end;
    
    fprintf(output_fptr, "\n****** SUBSURFACE LATERAL FLOW NASH CASCADE VOLUME BALANCE ******\n");
    fprintf(output_fptr, " Volume into subsurface Nash cascade= %8.4lf m\n", volbal->vol_in_subsurf_nash);
    fprintf(output_fptr, " Volume out subsurf. Nash cascade   = %8.4lf m\n", volbal->vol_out_subsurf_nash);
    fprintf(output_fptr, " Final volume in Nash cascade       = %8.4lf m\n", vol_nash_subsurface_end);
    fprintf(output_fptr, " Nash cascade residual              = %6.4e m\n", nash_residual);
    if (fabs(nash_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: SUBSURFACE LATERAL FLOW NASH CASCADE CONCEPTUAL RESERVOIR VOLUME BALANCE CHECK FAILED\n");
    }

    // GROUNDWATER VOLUME BALANCE
    double gw_residual = volbal->vol_in_gw_start + volbal->vol_to_gw - volbal->vol_from_gw - vol_gw_end;
    
    fprintf(output_fptr, "\n********* GROUNDWATER CONCEPTUAL RESERVOIR VOLUME BALANCE *******\n");
    fprintf(output_fptr, " Initial GW storage                 = %8.4lf m\n", volbal->vol_in_gw_start);
    fprintf(output_fptr, " Volume from soil to GW             = %8.4lf m\n", volbal->vol_to_gw);
    fprintf(output_fptr, " Volume from GW to outflow          = %8.4lf m\n", volbal->vol_from_gw);
    fprintf(output_fptr, " Final GW storage                   = %8.4lf m\n", vol_gw_end);
    fprintf(output_fptr, " GW residual                        = %6.4e m\n", gw_residual);
    if (fabs(gw_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: GROUNDWATER CONCEPTUAL RESERVOIR VOLUME BALANCE CHECK FAILED\n");
    }

    double volend          = vol_soil_end   + vol_gw_end    + vol_surface_end   + vol_nash_subsurface_end;
    double global_residual = volstart       + volbal->volin - volout_total - volend;
        
    fprintf(output_fptr, "\n********************* GLOBAL VOLUME BALANCE ********************* \n");
    fprintf(output_fptr, " Partitioning scheme: %s\n", partition_name);
    fprintf(output_fptr, " Surface routing scheme: %s\n", surface_name);
    fprintf(output_fptr, " Initial soil reservoir storage     = %8.4lf m\n", volbal->volstart_soil);
    fprintf(output_fptr, " Initial g.w. reservoir storage     = %8.4lf m\n", volbal->volstart_gw);
    fprintf(output_fptr, " Initial surface routing storage    = %8.4lf m\n", volbal->volstart_surface);
    fprintf(output_fptr, " Initial subsurface routing storage = %8.4lf m\n", volbal->volstart_subsurface);
    fprintf(output_fptr, " Initial storage total ------A----- = %8.4lf m\n", volstart);
    fprintf(output_fptr, " Total liquid water INPUT ---B----- = %8.4lf m\n", volbal->volin);
    fprintf(output_fptr, " Total PET                          = %8.4lf m\n", volbal->volin_PET);
    fprintf(output_fptr, " Total AET                          = %8.4lf m\n", total_AET_vol);
    fprintf(output_fptr, " Forest/root-zone AET               = %8.4lf m\n", volbal->vol_forest_aet);
    fprintf(output_fptr, " Bare-soil evaporation              = %8.4lf m\n", volbal->vol_bare_soil_evaporation);
    fprintf(output_fptr, " Total outflow                      = %8.4lf m\n", volbal->volout);
    fprintf(output_fptr, " Sum of outflow + AET -------C----- = %8.4lf m\n", total_AET_vol + volbal->volout);
    fprintf(output_fptr, " Final soil reservoir storage       = %8.4lf m\n", vol_soil_end);
    fprintf(output_fptr, " Final g.w. reservoir storage       = %8.4lf m\n", vol_gw_end);
    fprintf(output_fptr, " Final surface routing storage      = %8.4lf m\n", vol_surface_end);
    fprintf(output_fptr, " Final subsurface routing storage   = %8.4lf m\n", vol_nash_subsurface_end);
    fprintf(output_fptr, " Final storage total --------D----- = %8.4lf m\n", volend);
    fprintf(output_fptr, " Global residual   (A+B-C-D)        = %6.4e m\n", global_residual);

    if (volbal->volin > 0.0){
        fprintf(output_fptr, " Global percent error = %6.4e percent of inputs\n", global_residual / volbal->volin * 100.0);
    } else {
        fprintf(output_fptr, " Global pct. err: %6.4e percent of initial\n", global_residual / volbal->volstart * 100.0);
    }
    
    if (fabs(global_residual) > 1.0e-12) {
        fprintf(output_fptr, "!!! WARNING: GLOBAL VOLUME BALANCE CHECK FAILED\n");
    }

    fprintf(output_fptr, "\n*********** POTENTIAL ET/TRANSPIRATION (information only) *********\n");
    fprintf(output_fptr, " Total catchment PET (pre-partition)^          = %8.4lf m\n", volbal->volin_PET);
    fprintf(output_fptr, " Total vegetated-area potential transpiration^ = %8.4lf m\n",
            volbal->volin_PET * params->catchment_vegetated_fraction);
    fprintf(output_fptr, "\n ^ Not used in the volume balance calculation (A+B-C-D); shown for\n");
    fprintf(output_fptr, "   reference only. These may differ because the entire catchment might\n");
    fprintf(output_fptr, "   not be vegetated: potential transpiration = catchment PET x\n");
    fprintf(output_fptr, "   catchment_vegetated_fraction.\n");
}

// This function writes out a new .cfg file with updated model states
// The filename follows: prefix_YYYYMMDDHHMMSS.cf3 where the prefix is specified in the original config file
// e.g. prefix=configs/hotstart  and YYYYMMDDHHMMSS are the time stamp of the last forcings used in the
// simulation.   Other than state updating, all options, parameters, and output specs remain the same
//#######################
int write_hotstart_config(const CFE_CONFIG* cfg,
                          const cfe_state_struct* state,
                          const aorc_forcing_time_struct* forcing_time)
{
    if (cfg == NULL || state == NULL || forcing_time == NULL) {
        fprintf(stderr, "write_hotstart_config: null pointer input\n");
        return -1;
    }

    // Build simulation end timestamp YYYYMMDDhhmmss
    char sim_timestamp[32];
    snprintf(sim_timestamp, sizeof(sim_timestamp),
             "%04d%02d%02d%02d%02d%02d",
             forcing_time->year,
             forcing_time->month,
             forcing_time->day,
             forcing_time->hour,
             forcing_time->minute,
             forcing_time->second);

    // Build hotstart filename
    char filename[PATH_FILENAME_STRING_LENGTH];
    
    // conceivable that something could cause it to be larger than allocated size.  Check it.
    int ret = snprintf(filename, sizeof(filename), "%s.%s.cf3",
                       cfg->output_new_config_filename, sim_timestamp);
    if (ret >= sizeof(filename)) {
        fprintf(stderr, "ERROR: Hotstart filename exceeds buffer size\n");
        return -1;
    }

    FILE* hotstart_fptr = fopen(filename, "w");
    if (hotstart_fptr == NULL) {
        fprintf(stderr, "write_hotstart_config: cannot create hotstart file %s\n", filename);
        return -1;
    }

    // Get wall-clock time
    time_t now = time(NULL);
    struct tm lt;
#if defined(_WIN32) || defined(_WIN64)
    localtime_s(&lt, &now);
#else
    struct tm* tmp = localtime(&now);
    lt = *tmp;
#endif

    // Write a complete hotstart config file from scratch
    // This preserves the original structure but writes it fresh with updated state

    fprintf(hotstart_fptr, "# Conceptual Functional Equivalent (CFE) to the WRF-Hydro Based NOAA/NWS\n");
    fprintf(hotstart_fptr, "# National Water Model (versions 3.1 and earlier) stormflow generation function]\n"); 
    fprintf(hotstart_fptr, "# Hotstart Configuration File, CFE v%3.1f %s\n", CFE_VERSION, CFE_SUBVERSION_STRING);
    fprintf(hotstart_fptr, "# Compatible with CFE v3.\n");
    fprintf(hotstart_fptr, "# Generated on: %04d-%02d-%02d %02d:%02d:%02d\n",
                                    lt.tm_year + 1900, lt.tm_mon + 1, lt.tm_mday,
                                    lt.tm_hour, lt.tm_min, lt.tm_sec);
    fprintf(hotstart_fptr, "# This hotstart file written by CFE\n");
    fprintf(hotstart_fptr, "# Including model state updated at simulation end time: %s\n", sim_timestamp);
    fprintf(hotstart_fptr, "# Run using forcing file named: %s\n",cfg->control_input_forcing_filename);
    fprintf(hotstart_fptr, "\n");

    // write the CFE version number
    fprintf(hotstart_fptr,"cfe_config_version=%.1f[]\n",CFE_VERSION);

    // Model Controls
    fprintf(hotstart_fptr, "# Model Controls\n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "control_model_timestep_h=%.4f [h]\n", cfg->timestep_h);
    fprintf(hotstart_fptr, "control_input_forcing_filename=%s\n", cfg->control_input_forcing_filename); // could be BMI or CLI
    fprintf(hotstart_fptr, "control_total_num_simulation_timesteps=%d[]\n", cfg->total_timesteps);
    fprintf(hotstart_fptr, "control_verbosity=%d[]\n", cfg->verbosity);
    fprintf(hotstart_fptr, "control_ET_simulate_Priestley_Taylor=%.3f[]\n", cfg->et_alpha_pt);
    fprintf(hotstart_fptr, "control_soil_simulate_discrete_soil_moisture_true_false=%s\n",
            cfg->control_soil_simulate_discrete_soil_moisture_true_false ? "TRUE" : "FALSE");
    fprintf(hotstart_fptr, "control_soil_simulate_soil_evaporation=%s\n",
            cfg->control_soil_simulate_soil_evaporation ? "TRUE" : "FALSE");
    fprintf(hotstart_fptr, "control_ET_deepest_root_zone_discretization=%d\n",
            cfg->control_ET_deepest_root_zone_discretization);
    fprintf(hotstart_fptr, "control_soil_use_lookup_table_num_points=%d\n",
            cfg->control_soil_use_lookup_table_num_points);
    fprintf(hotstart_fptr, "control_soil_simulate_freeze_thaw_true_false=%s\n",
            cfg->control_soil_simulate_freeze_thaw_true_false ? "TRUE" : "FALSE");
    if (strlen(cfg->simulation_start_date) > 0) {
        fprintf(hotstart_fptr, "control_simulation_start_date=%s\n", cfg->simulation_start_date);
    }
    fprintf(hotstart_fptr, "\n");

    // Catchment Characteristics
    fprintf(hotstart_fptr, "# Catchment Characteristics\n");
    fprintf(hotstart_fptr, "#===========================\n");
    if (strlen(cfg->cat_id) > 0) {
        // Check if cat_id already has quotes, if so don't add more
        if (cfg->cat_id[0] == '"') {
            fprintf(hotstart_fptr, "catchment_id=%s\n", cfg->cat_id);
        } else {
            fprintf(hotstart_fptr, "catchment_id=\"%s\"\n", cfg->cat_id);
        }
    }
    if (cfg->cat_latitude != 0.0) {
        fprintf(hotstart_fptr, "catchment_latitude_decimal_degree=%.4f[decimaldegree]\n", cfg->cat_latitude);
    }
    if (cfg->cat_longitude != 0.0) {
        fprintf(hotstart_fptr, "catchment_longitude_decimal_degree=%.4f[decimaldegree]\n", cfg->cat_longitude);
    }
    if (cfg->cat_elev != 0.0) {
        fprintf(hotstart_fptr, "catchment_elevation=%.2f[m]\n", cfg->cat_elev);
    }
    if (cfg->cat_area_km2 != 0.0) {
        fprintf(hotstart_fptr, "catchment_area_km2=%e[km2]\n", cfg->cat_area_km2);
    }
    fprintf(hotstart_fptr, "catchment_impervious_fraction_0-1=%.6f[]\n", cfg->cat_impervious_fraction);
    fprintf(hotstart_fptr, "catchment_forested_fraction_0-1=%.6f[]\n", cfg->catchment_vegetated_fraction);
    fprintf(hotstart_fptr, "\n");

    // Soil Parameters
    fprintf(hotstart_fptr, "# Soil Parameters\n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "soil_depth_m=%.4f[m]\n", cfg->soil_depth_m);
    fprintf(hotstart_fptr, "soil_Clapp_Hornberger_exponent_b=%.4f[]\n", cfg->soil_Clapp_Hornberger_exponent_b);
    fprintf(hotstart_fptr, "soil_sat_hydraulic_conductivity_cm_per_h=%.4f[cm h-1]\n", 
            cfg->soil_sat_hydraulic_conductivity_cm_per_h);
    fprintf(hotstart_fptr, "soil_sat_capillary_head_cm=%.4f[cm]\n", cfg->soil_sat_capillary_head_cm);
    fprintf(hotstart_fptr, "soil_effective_porosity=%.4f[V V-1]\n", cfg->soil_effective_porosity);
    fprintf(hotstart_fptr, "soil_field_capacity_Pcap_over_Patm_0_1=%.4f[P P-1]\n", 
            cfg->soil_field_capacity_Pcap_over_Patm_0_1);
    fprintf(hotstart_fptr, "soil_reservoir_rate_const_to_subsurface_lateral_flow=%e[h-1]\n", 
            cfg->soil_reservoir_rate_const_to_subsurface_lateral_flow);
    fprintf(hotstart_fptr, "soil_to_gw_percolation_rate_limiter_0_to_1=%e[]\n",
            cfg->soil_to_gw_percolation_rate_limiter_0_to_1);
    fprintf(hotstart_fptr, "soil_ice_content_impervious_threshold=%.4f[]\n",
            cfg->soil_ice_content_impervious_threshold);
    fprintf(hotstart_fptr, "\n");

    // Initial soil moisture storage - UPDATE WITH CURRENT STATE
    fprintf(hotstart_fptr, "#----- Initial soil moisture storage (UPDATED FROM MODEL STATE)\n");
    if (cfg->control_soil_simulate_discrete_soil_moisture_true_false) {
        fprintf(hotstart_fptr, "state_soil_reservoir_init_discrete_storage_theta=");
        for (int i = 0; i < NDISC; i++) {
            if (i > 0) fprintf(hotstart_fptr, ", ");
            fprintf(hotstart_fptr, "%.17e", state->soil_discrete_storage_theta[i]);
        }
        fprintf(hotstart_fptr, " []\n");
    } else {
        fprintf(hotstart_fptr, "state_soil_reservoir_init_storage_m=%.17e[m]\n", state->soil_storage_m);
    }
    fprintf(hotstart_fptr, "\n");

    // Groundwater Parameters
    fprintf(hotstart_fptr, "# Groundwater Parameters \n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "gw_reservoir_max_storage_m=%.5f[m]\n", cfg->gw_reservoir_max_storage_m);
    fprintf(hotstart_fptr, "gw_discharge_coeff_m_per_timestep=%e[m h-1]\n", cfg->gw_discharge_coeff_m_per_timestep);
    fprintf(hotstart_fptr, "gw_discharge_exponent=%.4f[]\n", cfg->gw_discharge_exponent);
    fprintf(hotstart_fptr, "#----- Initial groundwater storage (UPDATED FROM MODEL STATE)\n");  
    fprintf(hotstart_fptr, "state_gw_reservoir_init_storage_m=%.17e[m]\n", state->gw_storage_m);
    fprintf(hotstart_fptr, "\n");

    // Partitioning scheme
    fprintf(hotstart_fptr, "# Partitioning of liquid water at land surface\n");
    fprintf(hotstart_fptr, "# must pick one: (SCHAAKE or XINANJIANG)\n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "#\n");
    
    if (string_compare_ignore_case(cfg->partitioning_scheme_name, "schaake") == 0) {
        fprintf(hotstart_fptr, "#----- SCHAAKE  (requires no other keywords or params)\n");
        fprintf(hotstart_fptr, "partitioning_scheme_name=SCHAAKE\n");
        fprintf(hotstart_fptr, "#\n");
    } else {
        // Xinanjiang
        fprintf(hotstart_fptr, "#----- XINANJIANG (requires following 3 keywords/params\n");
        fprintf(hotstart_fptr, "partitioning_scheme_name=XINANJIANG\n");
        fprintf(hotstart_fptr, "partitioning_Xinanjiang_tension_water_inflection_point=%.6f[V V-1]\n",
                cfg->soil_Xinanjiang_tension_water_inflection_point);
        fprintf(hotstart_fptr, "partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent=%.4f[]\n",
                cfg->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent);
        fprintf(hotstart_fptr, "partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent=%.4f[]\n",
                cfg->soil_Xinanjiang_free_water_soil_moist_distrib_exponent);
    }
    fprintf(hotstart_fptr, "\n");

    // Surface Routing Parameters — GIUH only (UPDATED FROM MODEL STATE)
    fprintf(hotstart_fptr, "# Surface Routing Parameters (GIUH)\n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "surface_routing_num_giuh_ordinates=%d\n", cfg->surface_routing_num_giuh_ordinates);
    fprintf(hotstart_fptr, "surface_routing_giuh_ordinates=");
    for (int i = 0; i < cfg->surface_routing_num_giuh_ordinates; i++) {
        if (i > 0) fprintf(hotstart_fptr, ", ");
        fprintf(hotstart_fptr, "%.6f", cfg->surface_routing_giuh_ordinates[i]);
    }
    fprintf(hotstart_fptr, " [m m-1]\n");
    fprintf(hotstart_fptr, "#----- Initial GIUH convolution queue storage (UPDATED FROM MODEL STATE)\n");
    fprintf(hotstart_fptr, "state_surface_routing_init_giuh_convolution_queue_m=");
    for (int i = 0; i < cfg->surface_routing_num_giuh_ordinates; i++) {
        if (i > 0) fprintf(hotstart_fptr, ", ");
        fprintf(hotstart_fptr, "%.17e", state->giuh_queue_m[i]);
    }
    fprintf(hotstart_fptr, "[m]\n");
    fprintf(hotstart_fptr, "\n");

    // Subsurface Lateral Flow Parameters (UPDATED FROM MODEL STATE)
    fprintf(hotstart_fptr, "# Subsurface Lateral Flow Parameters\n");
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "subsurface_routing_nash_reservoir_time_constant_k=%.6f[h-1]\n", cfg->subsurface_routing_nash_K);
    fprintf(hotstart_fptr, "#----- Initial subsurface Nash cascade storage (UPDATED FROM MODEL STATE)\n");
    fprintf(hotstart_fptr, "state_subsurface_routing_init_nash_cascade_storage_m=");
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        if (i > 0) fprintf(hotstart_fptr, ",");
        fprintf(hotstart_fptr, "%.17e", state->nash_subsurface_storage_m[i]);
    }
    fprintf(hotstart_fptr, "[m]\n");
    fprintf(hotstart_fptr, "\n");

    // PET temperature state (UPDATED FROM MODEL STATE)
    fprintf(hotstart_fptr, "#===========================\n");
    fprintf(hotstart_fptr, "#----- PET temperature state (UPDATED FROM MODEL STATE)\n");
    fprintf(hotstart_fptr, "state_skin_temperature_k=%.17e[K]\n",
            state->pet_temperature_state.skin_temperature_k);
    fprintf(hotstart_fptr, "state_upper_soil_temperature_k=%.17e[K]\n",
            state->pet_temperature_state.upper_soil_temperature_k);
    fprintf(hotstart_fptr, "state_estimated_annual_air_temperature_k=%.17e[K]\n",
            state->pet_temperature_state.estimated_annual_air_temperature_k);
    fprintf(hotstart_fptr, "state_air_temperature_time_integral_k_s=%.17e[K s]\n",
            state->pet_temperature_state.air_temperature_time_integral_k_s);
    fprintf(hotstart_fptr, "state_accumulated_time_s=%.17e[s]\n",
            state->pet_temperature_state.accumulated_time_s);
    fprintf(hotstart_fptr, "state_pet_initialized=%d[]\n",
            state->pet_temperature_state.initialized);
    fprintf(hotstart_fptr, "\n");

    // Output Configuration (copy from original)
    fprintf(hotstart_fptr, "# Disk Output Configuration Options\n");
    fprintf(hotstart_fptr, "#===========================\n");
    if (strlen(cfg->output_time_standard_format) > 0) {
        fprintf(hotstart_fptr, "output_time_standard_format=%s\n", cfg->output_time_standard_format);
    }
    if (strlen(cfg->output_file_delimiter) > 0) {
        fprintf(hotstart_fptr, "output_file_delimiter=%s\n", cfg->output_file_delimiter);
    }
    if (strlen(cfg->output_value_format) > 0) {
        fprintf(hotstart_fptr, "output_value_format=\"%s\"\n", cfg->output_value_format);  // \" is literal quote
    }    
    fprintf(hotstart_fptr, "\n");

    // Optional output files
    fprintf(hotstart_fptr, "# Optional output files to write\n");
    fprintf(hotstart_fptr, "#===========================\n");
    if (strlen(cfg->output_path_name) > 0) {
        fprintf(hotstart_fptr, "output_path_name=\"%s\"\n", cfg->output_path_name);
    }
    if (strlen(cfg->output_internal_fluxes_filename) > 0) {
        fprintf(hotstart_fptr, "output_internal_fluxes_m_per_timestep_filename=\"%s\"\n", cfg->output_internal_fluxes_filename);
    }
    if (strlen(cfg->output_internal_storages_filename) > 0) {
        fprintf(hotstart_fptr, "output_internal_storages_m_per_timestep_filename=\"%s\"\n", cfg->output_internal_storages_filename);
    }
    if (strlen(cfg->output_volume_balance_filename) > 0) {
        fprintf(hotstart_fptr, "output_volume_balance_filename=\"%s\"\n", cfg->output_volume_balance_filename);
    }
    if (strlen(cfg->output_soil_moisture_theta_filename) > 0) {
        fprintf(hotstart_fptr, "output_soil_moisture_theta_filename=\"%s\"\n", cfg->output_soil_moisture_theta_filename);
    }
    if (strlen(cfg->output_discharge_filename) > 0) {
        fprintf(hotstart_fptr, "output_discharge_m_per_timestep_filename=\"%s\"\n", cfg->output_discharge_filename);
    }
    if (strlen(cfg->output_total_discharge_m3_per_sec_filename) > 0) {
        fprintf(hotstart_fptr, "output_total_discharge_m3_per_sec_filename=\"%s\"\n", cfg->output_total_discharge_m3_per_sec_filename);
    }
    fprintf(hotstart_fptr, "\n");

    // Hotstart generation
    fprintf(hotstart_fptr, "# HOTSTART\n");
    fprintf(hotstart_fptr, "#===========================\n");
    if (strlen(cfg->output_new_config_filename) > 0) {
        fprintf(hotstart_fptr, "output_new_config_filename_prefix=\"%s\"\n", cfg->output_new_config_filename);
    }

    fclose(hotstart_fptr);

    printf("Hotstart config written to %s\n", filename);
    return 0;
}

// Time utility functions
//########################
int calculate_julian_day(int year, int month, int day) {
    // Standard Julian day calculation algorithm
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    return jdn;
}

double calculate_julian_date(int year, int month, int day, int hour, int minute, int second) {
    int jdn = calculate_julian_day(year, month, day);
    
    // Add fractional day for time
    double fraction = (hour + minute/60.0 + second/3600.0) / 24.0;
    
    return (double)jdn + fraction;
}

#ifdef __cplusplus
}
#endif
