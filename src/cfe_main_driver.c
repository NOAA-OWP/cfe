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

//mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
// This file is the standalone main() CFE model driver.  It does not use BMI- it uses straight
// function calls.  It is actually really pretty well written and thought out.   It uses many 
// other functions contained in cfe_driver_utils.c and cfe_helpers.c.  It also calls the parser
// utilities in parser_helpers.c   It parses config files, analyzes forcing files, runs CFE using
// a stateless BMI-like function call, writes necessary outputs, and calculates volume balance
// summary info at the end of the simulation.  FLO 9/2025
//mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "cfe_config.h"
#include "cfe_helpers.h"
#include "cfe_context.h"               // context create/update/destroy + setters/getters
#include "cfe_types.h"                 // All data types used in cfe model, outputs, etc.
#include "cfe.h"                       // includes volbal structure definition
#include "cfe_driver_utils.h"
#include "cfe_pet_priestley_taylor.h"
#include "cfe_soil_skin_temperature.h"
#include "calculate_bare_soil_evap.h"

#define TIME_STRING_LENGTH 64

//######
int main(int argc, char* argv[])
{
    // Parse command line arguments
    cfe_cmdline_args_struct cmdline_args;
    if (parse_command_line(argc, argv, &cmdline_args) != 0) {
        return 1;
    }

    // Initialize data structures
    CFE_CONFIG            config  = {0};  //<--- this contains everything read in from a config file
    cfe_options_struct    options = {0};
    cfe_parameters_struct params  = {0};
    cfe_state_struct      state   = {0};
    cfe_forcing_struct    forcing = {0};
    cfe_outputs_struct    outputs = {0};
    cfe_volbal_struct     volbal  = {0};


    set_state_defaults(&state);

    // 1. Read and parse configuration file
    double cfg_version = read_cfe_config_version(cmdline_args.cfg_path);
    if (fabs(cfg_version) < 1.0e-04) {
        fprintf(stderr,
            "WARNING: No cfe_config_version found in: %s\n"
            "         Add cfe_config_version=3.0 to the config file to suppress this warning.\n"
            "         If this is a legacy v2 config, use cfe_migrate_config to convert it.\n",
            cmdline_args.cfg_path);
        cfg_version = 3.0;
    }

    int status_flag = parse_config_driver(cmdline_args.cfg_path, cfg_version, &config, &params, &options);
    if (status_flag != 0) {
        fprintf(stderr, "RUNTIME ERROR: Failed to parse cfe config file.\n");
        return 1;
    }

    // 2. Apply command line overrides to options
    //    Arguments coming in through the command line override arguments set through a config file
    if (check_command_line_overrides(&options, &cmdline_args) != 0) {
        return 1;
    }
    
    if (options.verbosity > 0) fprintf(stderr, "Parsed CFE version %.1f config file\n", cfg_version);
    
    // 3. Initialize CFE state
    if (cfe_initialize(&params, &options, &state) != 0) {
        fprintf(stderr, "ERROR: Failed to initialize CFE state\n");
        return 1;
    }

    if (config.state_pet_initialized) {
        state.pet_temperature_state.skin_temperature_k =
            config.state_skin_temperature_k;
        state.pet_temperature_state.upper_soil_temperature_k =
            config.state_upper_soil_temperature_k;
        state.pet_temperature_state.estimated_annual_air_temperature_k =
            config.state_estimated_annual_air_temperature_k;
        state.pet_temperature_state.air_temperature_time_integral_k_s =
            config.state_air_temperature_time_integral_k_s;
        state.pet_temperature_state.accumulated_time_s =
            config.state_accumulated_time_s;
        state.pet_temperature_state.initialized =
            config.state_pet_initialized;
    }

    // 4. Handle forcing file setup
    char forcing_file_path[PATH_FILENAME_STRING_LENGTH] = "";
    int forcing_from_cmdline = 0;
    int have_forcing = 0;
    FILE* forcing_fptr = NULL; 
    aorc_cols_t cols; 
    aorc_forcing_time_struct forcing_time = (aorc_forcing_time_struct){0};

    // 5. Determine forcing file path
    if (cmdline_args.forcing_path != NULL && strlen(cmdline_args.forcing_path) > 0) {
        strncpy(forcing_file_path, cmdline_args.forcing_path, sizeof(forcing_file_path) - 1);
        forcing_file_path[sizeof(forcing_file_path) - 1] = '\0';
        forcing_from_cmdline = 1;
    }
    else if (!cmdline_args.run_without_forcing && strlen(options.input_forcing_filename) > 0) {
        strncpy(forcing_file_path, options.input_forcing_filename, sizeof(forcing_file_path) - 1);
        forcing_file_path[sizeof(forcing_file_path) - 1] = '\0';
    }
    
    // 6. Analyze forcing file iff provided (sometimes not provide ffor testing purposes (e.g. --dryrun flag or BMI)
    have_forcing = (!cmdline_args.run_without_forcing) && (forcing_file_path[0] != '\0');
    if (have_forcing) {
        int flag = analyze_forcing_file(forcing_file_path, &forcing_time, options.verbosity);
        if (flag != 0) {
            if (forcing_from_cmdline) {
                fprintf(stderr, "ERROR: Could not analyze command-line forcing file: %s\n", forcing_file_path);
            } else {
                fprintf(stderr, "ERROR: Could not analyze forcing file from config: %s\n", forcing_file_path);
            }
            fprintf(stderr, "analyze_forcing_file() returned %d\n", flag);
            return 1;
        }
        
        if (open_forcing_file(forcing_file_path, &forcing_fptr, &cols) != 0) {
            fprintf(stderr, "ERROR: cannot open forcing file: %s\n", forcing_file_path);
            return 1; 
        } 
    }
    
    // 7.  Print entire input deck iff high verbosity
    if (options.verbosity > 1) { 
        print_cfe_input_debug(&options, &params, &state, cmdline_args.cfg_path, 
                             cmdline_args.forcing_path, cmdline_args.qout_path, cmdline_args.volbal_path); 
    }

    // 8. Open output files (config + command line overrides)
    FILE* q_fptr = NULL;       // discharge per unit area per unit timestep (m per timestep)
    FILE* Q_fptr = NULL;       // total discharge (m3/s) - instantaneous
    FILE* fluxes_fptr = NULL;
    FILE* storages_fptr = NULL;
    FILE* volbal_fptr = NULL;
    FILE* theta_fptr = NULL;
    FILE* warnings_fptr = NULL;

    if (open_output_files(&options, &cmdline_args, &q_fptr, &Q_fptr, &fluxes_fptr, 
                         &storages_fptr, &volbal_fptr, &theta_fptr) != 0) {
        // on error, close all files
        goto AA; // Gotta have at lease one ffor my Fortran homeys
    }

    // 9. Determine time step.  Since code cannot interpolate/extrapolate forcings in time,
    //    The model time step and forcing time step must be the same.
    if (have_forcing && forcing_time.delta_t_seconds > 0.0) {
        if (options.time_step_seconds <= 0.0) {
            // Config time step not set - use forcing file value
            options.time_step_seconds = forcing_time.delta_t_seconds;
            if (options.verbosity > 0) {
                fprintf(stderr, "Using delta t from forcings file: %d s\n", forcing_time.delta_t_seconds);
            }
        } else {
            // Config time step is set - enforce consistency
            if (options.time_step_seconds != forcing_time.delta_t_seconds) {
                fprintf(stderr, "ERROR: Time step mismatch between config (%d s) and forcing file (%d s)\n",
                        options.time_step_seconds, forcing_time.delta_t_seconds);
                fprintf(stderr, "CFE cannot interpolate forcings data. Forcing data and model time steps must be the same.\n");
                goto AA;  
            }
        }
    }
    if (options.time_step_seconds <= 0.0) {
        fprintf(stderr, "ERROR: time_step_seconds <= 0; cannot proceed\n");
        goto AA;
    }

    // 10. Write output headers
    int yes_write_output = FALSE;
    if (q_fptr || fluxes_fptr || storages_fptr || volbal_fptr || theta_fptr) {
        yes_write_output = TRUE;
        const char* delimiter = get_delimiter_string(options.output_file_delimiter);
        write_output_headers(&options, options.output_time_standard_format, delimiter,
                            q_fptr, Q_fptr, fluxes_fptr, storages_fptr, theta_fptr);
    }

    // 11. Initialize volume balance
    memset(&volbal, 0, sizeof(volbal));
    cfe_initialize_volume_balance(&params, &options, &state, &volbal);

    if (options.verbosity > 0) {
        printf("DEBUG: After CFE_initialize_volume_balance stored volbal.volstart_surface = %.6f\n", volbal.volstart_surface);
    }
    // 12. Main simulation loop
    int nsteps;
    if(options.num_timesteps == 0) {
        if (cmdline_args.run_without_forcing) {
            // Dryrun mode - use default timesteps
            nsteps = 120;
            if (options.verbosity > 0) {
                fprintf(stderr, "Dryrun mode: using default %d timesteps\n", nsteps);
            }
        } else if (have_forcing) {
            // Use forcing file length
            nsteps = forcing_time.num_valid_lines;
            if (options.verbosity > 0) {
                fprintf(stderr, "Using forcing file length: %d timesteps\n", nsteps);
            }
        } else {
            fprintf(stderr, "ERROR: num_timesteps=0 requires either a forcing file or --dryrun flag\n");
            goto AA;
        }
    } else {
        nsteps = options.num_timesteps; 
    }
    char forcing_time_str[TIME_STRING_LENGTH];
   
    for (int t = 0; t < nsteps; t++) {   //<------------------ TIME LOOP -----------------------<<

        // 12a. Read forcing data
        if (have_forcing) {
            read_next_forcing_aorc(forcing_fptr, &cols, options.time_step_seconds, &forcing, forcing_time_str);
        } else {
            memset(&forcing, 0, sizeof(forcing));
        }

        // 12b. Update skin temperature, PET, and bare-soil evaporation
        cfe_update_pet_and_bare_soil_evap(
            &forcing, &params, &state, &options);

        // 12c. Run CFE model step
        cfe_step(&params, &options, &state, &forcing, (double)options.time_step_seconds, &outputs, &volbal);

        // Print model configuration (first timestep only)
        if (t == 0 && options.verbosity > 0) {
            printf("\n======================== CFE MODEL v3.0 ======================\n");
            printf("Conceptual Functional Equivalent to the NWS National Water Model\n");
            printf("(versions 3.1 and earlier) stormflow generation and catchment-scale\n");
            printf("routing routines.  (C) NOAA/NWS 2025.\n");
            printf("---- Configured Options-------------------------------------\n");
            printf("First time through time loop in main() for %d timesteps.\n", nsteps);
            printf("Soil moisture module:                          %s\n", options.simulate_discrete_soil_moisture ? "Discrete soil moisture balance module" : "Conceptual linear reservoir");
            printf("Land surface liquid water partitioning scheme: %s\n", options.liquid_partitioning_scheme == PARTITION_SCHAAKE ? "Schaake" : "Xinanjiang");
            printf("Surface routing scheme:                        %s\n", options.surface_routing_scheme == SURF_ROUTE_GIUH ? "GIUH" : "Nash Cascade");
            printf("============================================================\n");
        }

        // 12d. Verbose output
        if (options.verbosity > 1) {
            printf("Step %d rain=%.2f mm et=%.2f mm Qout=%.6f m\n", t, 
                   forcing.rainfall_depth_m * 1000.0, forcing.et_potential_m * 1000.0, outputs.qout_m);
        }

        // 12e. Write outputs
        if (yes_write_output) {
            write_all_outputs(t, &options, &params, &state,  &outputs,  &forcing, &forcing_time, &volbal, 
                              q_fptr, Q_fptr,  fluxes_fptr,  storages_fptr, volbal_fptr,  theta_fptr);
        }
        
        // 12f. If end of time loop, save date/time stamp of last forcing data point
        if (t == (nsteps-1)){
            parse_time_string(forcing_time_str, &forcing_time);
        }
    }     //<--------------------------------------------------- END OF TIME LOOP --------------<<

    // 13. Write final volume balance summary
    double soil_disc_dz[4] = {0.10, 0.30, 0.60, 1.00};
    if (volbal_fptr) {
        write_volume_balance_summary(volbal_fptr, &options, &volbal, &state, &params, soil_disc_dz);
    }
    
    if (options.verbosity > 0) {
        write_volume_balance_summary(stdout, &options, &volbal, &state, &params, soil_disc_dz);
    }


    // 14. Write hotstart config file iff needed
    if (strlen(config.output_new_config_filename) > 0) {
        int result = write_hotstart_config(&config, &state, &forcing_time);
        if (result != 0) {
            fprintf(stderr, "WARNING: failed to write hotstart config file.  Reason: UNKNOWN.\n");
        }
    }

    // 15. Cleanup and exit
    
AA:
    close_output_files(q_fptr, Q_fptr, fluxes_fptr, storages_fptr, volbal_fptr, 
                      theta_fptr, warnings_fptr, forcing_fptr);
    return 0;
}


