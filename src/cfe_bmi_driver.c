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

// THIS FILE IS THE main() CFE BMI driver.
//
// FLO used claude.ai to adapt the standalone stateless CFE 
// driver program into this BMI driver program.  FLO 9/2025

// command line arguments accepted.  See cfe_driver_utils.h.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bmi_cfe.h"
#include "cfe_driver_utils.h"
#include "cfe_pet_priestley_taylor.h"

#define TIME_STRING_LENGTH 64

// Forward declarations for static functions
static void write_output_headers_simple(const char* time_format, const char* delimiter,
                                       FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr,
                                       FILE* storages_fptr, FILE* theta_fptr, 
                                       int simulate_discrete_soil_moisture);

static void write_outputs_from_bmi(int timestep, double catchment_area_km2, double dt_seconds,
                                  double qout_m, double surface_runoff_m, double lateral_flow_m,
                                  double baseflow_m, double soil_storage_m, double gw_storage_m,
                                  double soil_theta[4], int simulate_discrete_soil_moisture,
                                  FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr,
                                  FILE* storages_fptr, FILE* theta_fptr);

static void write_volume_balance_summary_bmi(FILE* output_fptr, const char* config_file,
                                            const char* forcing_file, int timesteps_run,
                                            Bmi* model, int simulate_discrete_soil_moisture,
                                            double calculated_residual);

static double total_input_m;
static double total_output_m;
static double total_storage_start_m;
static double total_storage_end_m;

//######
int main(int argc, char* argv[])
{
    // Parse command line arguments (only -c config file and -v verbosity supported)
    cfe_cmdline_args_struct cmdline_args;
    if (parse_command_line(argc, argv, &cmdline_args) != 0) {
        return 1;
    }

    // Initialize BMI model
    Bmi* model = (Bmi*)malloc(sizeof(Bmi));
    if (model == NULL) {
        fprintf(stderr, "ERROR: Failed to allocate Bmi struct\n");
        return 1;
    }
    register_bmi_cfe(model);

    // Initialize the model with config file
    if (model->initialize(model, cmdline_args.cfg_path) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: Failed to initialize BMI CFE model with config: %s\n", cmdline_args.cfg_path);
        delete_bmi_cfe(model);
        return 1;
    }

    // Apply verbosity override from command line
    if (cmdline_args.command_line_verbosity >= 0) {
        model->set_value(model, "verbosity", &cmdline_args.command_line_verbosity);
    }

    // Get model configuration information from BMI instead of parsing config file
    int verbosity = 0;
    int simulate_discrete_soil_moisture = 0;
    double catchment_area_km2 = 0.0;
    double alpha_pt = 0.0;
    
    model->get_value(model, "verbosity", &verbosity);
    model->get_value(model, "config_simulate_discrete_soil_moisture", &simulate_discrete_soil_moisture);
    model->get_value(model, "param_catchment_area_km2", &catchment_area_km2);
    
    // Try to get alpha_pt - if this fails, assume PET calculation is disabled
    // Note: You may need to add "param_alpha_pt" to the BMI interface
    int enable_pt = FALSE;
    if (model->get_value(model, "param_alpha_pt", &alpha_pt) == BMI_SUCCESS && alpha_pt > 0) {
        enable_pt = TRUE;
    }
    
    if (verbosity > 0) fprintf(stderr, "BMI CFE driver initialized successfully\n");

    // Handle forcing file setup - BMI driver still needs to analyze and read forcing files
    char forcing_file_path[PATH_FILENAME_STRING_LENGTH] = "";
    int forcing_from_cmdline = 0;
    int have_forcing = 0;
    FILE* forcing_fptr = NULL; 
    aorc_cols_t cols; 
    aorc_forcing_time_struct forcing_time = (aorc_forcing_time_struct){0};

    // Determine forcing file path - command line overrides BMI model
    if (cmdline_args.forcing_path != NULL && strlen(cmdline_args.forcing_path) > 0) {
        strncpy(forcing_file_path, cmdline_args.forcing_path, sizeof(forcing_file_path) - 1);
        forcing_file_path[sizeof(forcing_file_path) - 1] = '\0';
        forcing_from_cmdline = 1;
        
        // Override the BMI model's forcing file setting
        model->set_value(model, "forcing_file_path", forcing_file_path);
    }
    else {
        // Get forcing file from BMI model
        model->get_value(model, "forcing_file_path", forcing_file_path);
        if (!cmdline_args.run_without_forcing && strlen(forcing_file_path) == 0) {
            fprintf(stderr, "ERROR: No forcing file specified in config or command line\n");
            delete_bmi_cfe(model);
            return 1;
        }
    }
    
    // Analyze forcing file if provided (required for BMI driver to know timesteps and delta_t)
    have_forcing = (!cmdline_args.run_without_forcing) && (forcing_file_path[0] != '\0');
    if (!have_forcing) {
        fprintf(stderr, "ERROR: BMI driver requires a forcing file to be specified\n");
        delete_bmi_cfe(model);
        return 1;
    }

    int flag = analyze_forcing_file(forcing_file_path, &forcing_time, verbosity);
    if (flag != 0) {
        if (forcing_from_cmdline) {
            fprintf(stderr, "ERROR: Could not analyze command-line forcing file: %s\n", forcing_file_path);
        } else {
            fprintf(stderr, "ERROR: Could not analyze forcing file from BMI model: %s\n", forcing_file_path);
        }
        fprintf(stderr, "analyze_forcing_file() returned %d\n", flag);
        delete_bmi_cfe(model);
        return 1;
    }
    
    if (open_forcing_file(forcing_file_path, &forcing_fptr, &cols) != 0) {
        fprintf(stderr, "ERROR: cannot open forcing file: %s\n", forcing_file_path);
        delete_bmi_cfe(model);
        return 1; 
    }

    // Get model time step and verify consistency with forcing
    double model_dt_seconds;
    if (model->get_time_step(model, &model_dt_seconds) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: Failed to get model time step\n");
        close_output_files(NULL, NULL, NULL, NULL, NULL, NULL, NULL, forcing_fptr);
        delete_bmi_cfe(model);
        return 1;
    }

    // Verify forcing and model timesteps are consistent
    if (fabs(model_dt_seconds - forcing_time.delta_t_seconds) > 1.0e-3) {
        fprintf(stderr, "ERROR: Time step mismatch between model (%d s) and forcing file (%d s)\n",
                (int)model_dt_seconds, forcing_time.delta_t_seconds);
        fprintf(stderr, "CFE cannot interpolate forcings data. Forcing data and model time steps must be the same.\n");
        close_output_files(NULL, NULL, NULL, NULL, NULL, NULL, NULL, forcing_fptr);
        delete_bmi_cfe(model);
        return 1;
    }

    // Open output files based on command line overrides
    // The BMI driver uses command line args to override output files
    FILE* q_fptr = NULL;       
    FILE* Q_fptr = NULL;       
    FILE* fluxes_fptr = NULL;
    FILE* storages_fptr = NULL;
    FILE* volbal_fptr = NULL;
    FILE* theta_fptr = NULL;
    FILE* warnings_fptr = NULL;

    // Open files based on command line arguments
    if (cmdline_args.qout_path) {
        q_fptr = fopen(cmdline_args.qout_path, "w");
        if (!q_fptr) {
            fprintf(stderr, "ERROR: Cannot open discharge output file: %s (does the output directory exist?)\n", cmdline_args.qout_path);
            close_output_files(NULL, NULL, NULL, NULL, NULL, NULL, NULL, forcing_fptr);
            delete_bmi_cfe(model);
            return 1;
        }
        if (verbosity > 0) printf("Writing discharge to: %s\n", cmdline_args.qout_path);
    }

    if (cmdline_args.volbal_path) {
        volbal_fptr = fopen(cmdline_args.volbal_path, "w");
        if (!volbal_fptr) {
            fprintf(stderr, "ERROR: Cannot open volume balance output file: %s (does the output directory exist?)\n", cmdline_args.volbal_path);
            close_output_files(q_fptr, Q_fptr, fluxes_fptr, storages_fptr, NULL, theta_fptr, warnings_fptr, forcing_fptr);
            delete_bmi_cfe(model);
            return 1;
        }
        if (verbosity > 0) printf("Writing volume balance to: %s\n", cmdline_args.volbal_path);
    }

    if (cmdline_args.fluxes_path) {
        fluxes_fptr = fopen(cmdline_args.fluxes_path, "w");
        if (!fluxes_fptr) {
            fprintf(stderr, "ERROR: Cannot open fluxes output file: %s (does the output directory exist?)\n", cmdline_args.fluxes_path);
            close_output_files(q_fptr, Q_fptr, NULL, storages_fptr, volbal_fptr, theta_fptr, warnings_fptr, forcing_fptr);
            delete_bmi_cfe(model);
            return 1;
        }
        if (verbosity > 0) printf("Writing fluxes to: %s\n", cmdline_args.fluxes_path);
    }

    if (cmdline_args.stores_path) {
        storages_fptr = fopen(cmdline_args.stores_path, "w");
        if (!storages_fptr) {
            fprintf(stderr, "ERROR: Cannot open storages output file: %s (does the output directory exist?)\n", cmdline_args.stores_path);
            close_output_files(q_fptr, Q_fptr, fluxes_fptr, NULL, volbal_fptr, theta_fptr, warnings_fptr, forcing_fptr);
            delete_bmi_cfe(model);
            return 1;
        }
        if (verbosity > 0) printf("Writing storages to: %s\n", cmdline_args.stores_path);
    }

    if (cmdline_args.thetas_path) {
        if (!simulate_discrete_soil_moisture) {
            fprintf(stderr, "WARNING: -t option ignored because discrete soil moisture is not enabled\n");
        } else {
            theta_fptr = fopen(cmdline_args.thetas_path, "w");
            if (!theta_fptr) {
                fprintf(stderr, "ERROR: Cannot open theta output file: %s (does the output directory exist?)\n", cmdline_args.thetas_path);
                close_output_files(q_fptr, Q_fptr, fluxes_fptr, storages_fptr, volbal_fptr, NULL, warnings_fptr, forcing_fptr);
                delete_bmi_cfe(model);
                return 1;
            }
            if (verbosity > 0) printf("Writing soil moisture thetas to: %s\n", cmdline_args.thetas_path);
        }
    }

    // Write output headers if we have any output files
    int yes_write_output = FALSE;
    if (q_fptr || fluxes_fptr || storages_fptr || volbal_fptr || theta_fptr) {
        yes_write_output = TRUE;
        
        // Use simple defaults for BMI driver output format
        const char* delimiter = " ";  // space delimiter
        const char* time_format = "timestep";
        
        write_output_headers_simple(time_format, delimiter, q_fptr, Q_fptr, fluxes_fptr, 
                                   storages_fptr, theta_fptr, simulate_discrete_soil_moisture);
    }

    // Determine number of timesteps to run
    int nsteps = forcing_time.num_valid_lines;
    if (verbosity > 0) {
        fprintf(stderr, "Using forcing file length: %d timesteps\n", nsteps);
    }

    char forcing_time_str[TIME_STRING_LENGTH];
    cfe_forcing_struct forcing = {0};

    // Print model configuration (before time loop)
    if (verbosity > 0) {
        printf("\n=================== CFE BMI MODEL v3.0 ====================\n");
        printf("Conceptual Functional Equivalent to the NWS National Water Model\n");
        printf("(versions 3.1 and earlier) stormflow generation and catchment-scale\n");
        printf("routing routines.  (C) NOAA/NWS 2025.\n");
        printf("---- BMI Driver Configuration ---------------------------\n");
        printf("Running %d timesteps via BMI interface.\n", nsteps);
        printf("Soil moisture module:                          %s\n", simulate_discrete_soil_moisture ? "Discrete soil moisture balance module" : "Conceptual linear reservoir");
        printf("PET calculation:                               %s\n", enable_pt ? "Priestley-Taylor" : "Provided in forcing data or disabled");
        printf("=========================================================\n");
    }
   
    double timestep_storage_start_m = 0.0;
    double timestep_storage_end_m   = 0.0; 
    // Main simulation loop - BMI style
    for (int t = 0; t < nsteps; t++) {   //<------------------- START OF TIME LOOP ---------------<<

        // Read forcing data from AORC file
        if (read_next_forcing_aorc(forcing_fptr, &cols, (int)model_dt_seconds, &forcing, forcing_time_str) <= 0) {
            fprintf(stderr, "ERROR: Failed to read forcing data at timestep %d\n", t);
            break;
        }

        // BMI inputs are rates (m/s); convert depth values from forcing reader
        // P-T PET (if enabled) is calculated inside cfe_context_update()
        double rain_rate = forcing.rainfall_depth_m / model_dt_seconds;
        double pet_rate  = forcing.et_potential_m   / model_dt_seconds;

        if (model->set_value(model, "rainfall_depth_m", &rain_rate) != BMI_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set rainfall_depth_m at timestep %d\n", t);
            break;
        }

        if (model->set_value(model, "et_potential_m", &pet_rate) != BMI_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set et_potential_m at timestep %d\n", t);
            break;
        }

        model->set_value(model, "day_of_year", &forcing.day_of_year);

        // Pass AORC atmospheric forcing for Priestley-Taylor PET and skin temperature
        model->set_value(model, "TMP_2maboveground",  &forcing.TMP_2maboveground);
        model->set_value(model, "DSWRF_surface",       &forcing.DSWRF_surface);
        model->set_value(model, "DLWRF_surface",       &forcing.DLWRF_surface);
        model->set_value(model, "PRES_surface",        &forcing.PRES_surface);
        model->set_value(model, "SPFH_2maboveground",  &forcing.SPFH_2maboveground);
        model->set_value(model, "UGRD_10maboveground",  &forcing.UGRD_10maboveground);
        model->set_value(model, "VGRD_10maboveground",  &forcing.VGRD_10maboveground);

        // get from BMI or trap the initial storage iff known
        if(t==0) {
            model->get_value(model, "timestep_storage_start_m", &timestep_storage_start_m);
        } else {
            timestep_storage_start_m = timestep_storage_end_m; // Use end from previous step
        }
        
        //############################### BMI UPDATE ###########################
        // Update model via BMI
        if (model->update(model) != BMI_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to update model at timestep %d\n", t);
            break;
        }

        // Get output values via BMI
        double qout_m                   = 0.0;
        double surface_runoff_m         = 0.0;
        double lateral_flow_m           = 0.0;
        double baseflow_m               = 0.0;
        double actual_et_m              = 0.0;
        double timestep_input_m         = 0.0;
        double timestep_output_m        = 0.0;

        model->get_value(model, "discharge_m", &qout_m);
        model->get_value(model, "surface_runoff_m", &surface_runoff_m);
        model->get_value(model, "lateral_flow_m", &lateral_flow_m);
        model->get_value(model, "baseflow_m", &baseflow_m);
        model->get_value(model, "actual_et_m", &actual_et_m);
        model->get_value(model, "timestep_storage_start_m", &timestep_storage_start_m);
        model->get_value(model, "timestep_input_m", &timestep_input_m);
        model->get_value(model, "timestep_output_m", &timestep_output_m);
        model->get_value(model, "timestep_storage_end_m", &timestep_storage_end_m);



        // Accumulate ffor final volume balance
        total_input_m += timestep_input_m;
        total_output_m += timestep_output_m;
        
        // Note: storage start/end are absolute values, not accumulated
        if (t == 0) total_storage_start_m = timestep_storage_start_m;  // at beginning of time step, only save 1st one
        total_storage_end_m = timestep_storage_end_m;                  // at end of time step

        double time_step_residual = timestep_storage_start_m + timestep_input_m - timestep_output_m - timestep_storage_end_m;
        // printf("Time step: %5d  residual = %15f\n", t, time_step_residual); // DEBUG
        
        // Get state information for output writing
        double soil_storage_m = 0.0;
        double gw_storage_m = 0.0;
        double soil_theta[4] = {0};
        
        if (simulate_discrete_soil_moisture) {
            model->get_value(model, "state_soil_moisture_theta", soil_theta);
        } else {
            model->get_value(model, "state_soil_storage_m", &soil_storage_m);
        }
        model->get_value(model, "state_gw_storage_m", &gw_storage_m);

        // Verbose output
        if (verbosity > 1) {
            printf("Step %d rain=%.2f mm et=%.2f mm Qout=%.6f m AET=%.6f m\n", t, 
                   forcing.rainfall_depth_m * 1000.0, forcing.et_potential_m * 1000.0, 
                   qout_m, actual_et_m);
        }

        // Write outputs using BMI-retrieved data
        if (yes_write_output) {
            write_outputs_from_bmi(t, catchment_area_km2, model_dt_seconds, qout_m, surface_runoff_m,
                                  lateral_flow_m, baseflow_m, soil_storage_m, gw_storage_m, soil_theta,
                                  simulate_discrete_soil_moisture, q_fptr, Q_fptr, fluxes_fptr,
                                  storages_fptr, theta_fptr);
        }

        // Save date/time stamp of last forcing data point ffor potential hotstart
        if (t == (nsteps-1)){
            parse_time_string(forcing_time_str, &forcing_time);
        }
    }  // <------------------ END OF TIME LOOP ----------------<<

    // Debug the volume balance components
    double vol_in, vol_out, vol_start, vol_end;
    model->get_value(model, "vol_in", &vol_in);     
    model->get_value(model, "vol_out", &vol_out);     // <- note these are not persistent from one BMI call to the next...
    model->get_value(model, "vol_start", &vol_start);
    model->get_value(model, "vol_end", &vol_end);
    
   // NOTE: uncomment these two lines to desmonstrate that the volbal struct is not persistent from one BMI call to the next.
   //  printf("DEBUG2: vol_start=%.12e vol_in=%.12e vol_out=%.12e vol_end=%.12e\n", 
   //        vol_start, vol_in, vol_out, vol_end);

double cumulative_residual = total_storage_start_m + total_input_m - total_output_m - total_storage_end_m;
double relative_residual_pct = 0.0;

    if      (total_input_m > 0.0)         relative_residual_pct = 100.0 * cumulative_residual/total_input_m;
    else if (total_storage_start_m > 0.0) relative_residual_pct = 100.0 * cumulative_residual/total_storage_start_m;
    else                                  relative_residual_pct = 0.0; // unquantifiable
    
    printf("Global volume balance based on driver program summation of per bmi-call quantities:\n");
    printf("total_storage_start_m:         %.15e\n", total_storage_start_m);
    printf("total_input_m:                 %.15e\n", total_input_m);
    printf("total_output_m:                %.15e\n", total_output_m);
    printf("total_storage_end_m:           %.15e\n", total_storage_end_m);
    printf("cumulative_residual:           %.15e\n", cumulative_residual);
    
    if      (total_input_m > 0.0)          printf("residual (%% of input):         %.2f%%\n", relative_residual_pct);  
    else if (total_storage_start_m > 0.0)  printf("residual (%% of init. storage): %.2f%%\n", relative_residual_pct);
    else                                   printf("residual percent incalculable: %.2f\n", relative_residual_pct);                          
           
    // Write final volume balance summary if requested
    if (volbal_fptr) {
        write_volume_balance_summary_bmi(volbal_fptr, cmdline_args.cfg_path, forcing_file_path, 
                                        nsteps, model, simulate_discrete_soil_moisture, cumulative_residual);
    }
    
    if (verbosity > 0) {
        double vol_balance_residual = 0.0;
        model->get_value(model, "vol_balance_residual_m", &vol_balance_residual);
        
        printf("\n=================== BMI SIMULATION COMPLETE ===================\n");
        printf("Completed %d timesteps via BMI interface\n", nsteps);
        printf("Volume balance residual: %.6e m\n", vol_balance_residual);
        printf("==============================================================\n");
    }

    // BMI driver doesn't support hotstart file generation
    // (Could be added by getting final states via BMI and creating config)

    // Cleanup and exit
    close_output_files(q_fptr, Q_fptr, fluxes_fptr, storages_fptr, volbal_fptr, 
                      theta_fptr, warnings_fptr, forcing_fptr);
    
    // Finalize BMI model
    model->finalize(model);
    delete_bmi_cfe(model);
    
    return 0;
}

// Simple output header writing function for BMI driver
static void write_output_headers_simple(const char* time_format, const char* delimiter,
                                       FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr,
                                       FILE* storages_fptr, FILE* theta_fptr, 
                                       int simulate_discrete_soil_moisture) {
    
    if (q_fptr) {
        fprintf(q_fptr, "# CFE BMI Driver Discharge Output\n");
        fprintf(q_fptr, "#%s%sdischarge_m_per_timestep\n", time_format, delimiter);
    }

    if (Q_fptr) {
        fprintf(Q_fptr, "# CFE BMI Driver Total Discharge Output\n");
        fprintf(Q_fptr, "#%s%sdischarge_m3_per_s\n", time_format, delimiter);
    }

    if (fluxes_fptr) {
        fprintf(fluxes_fptr, "# CFE BMI Driver Internal Fluxes Output\n");
        fprintf(fluxes_fptr, "#%s%ssurface_runoff_m%slateral_flow_m%sbaseflow_m\n",
                time_format, delimiter, delimiter, delimiter);
    }
    
    if (theta_fptr) {
        fprintf(theta_fptr, "# CFE BMI Driver Soil Moisture Theta Output\n");
        fprintf(theta_fptr, "#%s%stheta1%stheta2%stheta3%stheta4\n", 
                time_format, delimiter, delimiter, delimiter, delimiter);
    }
    
    if (storages_fptr) {
        if (!simulate_discrete_soil_moisture) {
            fprintf(storages_fptr, "# CFE BMI Driver Internal Storages Output\n");
            fprintf(storages_fptr, "#%s%ssoil_storage_m%sgw_storage_m\n",
                    time_format, delimiter, delimiter);
        } else {
            fprintf(storages_fptr, "# CFE BMI Driver Internal Storages Output\n");
            fprintf(storages_fptr, "#%s%stheta1%stheta2%stheta3%stheta4%sgw_storage_m\n", 
                    time_format, delimiter, delimiter, delimiter, delimiter, delimiter);
        }       
    }
}

// Write outputs using BMI-retrieved data
static void write_outputs_from_bmi(int timestep, double catchment_area_km2, double dt_seconds,
                                  double qout_m, double surface_runoff_m, double lateral_flow_m,
                                  double baseflow_m, double soil_storage_m, double gw_storage_m,
                                  double soil_theta[4], int simulate_discrete_soil_moisture,
                                  FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr,
                                  FILE* storages_fptr, FILE* theta_fptr) {
    
    const char* delimiter = " ";
    const char* format = "%.8e";
    
    // Write discharge output
    if (q_fptr) {
        fprintf(q_fptr, "%d%s", timestep, delimiter);
        fprintf(q_fptr, format, qout_m);
        fprintf(q_fptr, "\n");
    }

    if (Q_fptr) {
        // Convert discharge from m/timestep to m�/s
        double discharge_m3_per_sec = qout_m * catchment_area_km2 * 1.0e+06 / dt_seconds;
        fprintf(Q_fptr, "%d%s", timestep, delimiter);
        fprintf(Q_fptr, format, discharge_m3_per_sec);
        fprintf(Q_fptr, "\n");
    }
    
    // Write fluxes output
    if (fluxes_fptr) {
        fprintf(fluxes_fptr, "%d%s", timestep, delimiter);
        fprintf(fluxes_fptr, format, surface_runoff_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, format, lateral_flow_m);
        fprintf(fluxes_fptr, "%s", delimiter);
        fprintf(fluxes_fptr, format, baseflow_m);
        fprintf(fluxes_fptr, "\n");
    }

    // Write storages output
    if (storages_fptr) {
        if (!simulate_discrete_soil_moisture) {
            fprintf(storages_fptr, "%d%s", timestep, delimiter);
            fprintf(storages_fptr, format, soil_storage_m);
            fprintf(storages_fptr, "%s", delimiter);
            fprintf(storages_fptr, format, gw_storage_m);
            fprintf(storages_fptr, "\n");
        } else {
            fprintf(storages_fptr, "%d%s", timestep, delimiter);
            fprintf(storages_fptr, format, soil_theta[0]);
            fprintf(storages_fptr, "%s", delimiter);
            fprintf(storages_fptr, format, soil_theta[1]);
            fprintf(storages_fptr, "%s", delimiter);
            fprintf(storages_fptr, format, soil_theta[2]);
            fprintf(storages_fptr, "%s", delimiter);
            fprintf(storages_fptr, format, soil_theta[3]);
            fprintf(storages_fptr, "%s", delimiter);
            fprintf(storages_fptr, format, gw_storage_m);
            fprintf(storages_fptr, "\n");
        }
    }

    // Write theta output (if discrete soil moisture is enabled and file is open)
    if (theta_fptr && simulate_discrete_soil_moisture) {
        fprintf(theta_fptr, "%d%s", timestep, delimiter);
        fprintf(theta_fptr, format, soil_theta[0]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, format, soil_theta[1]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, format, soil_theta[2]);
        fprintf(theta_fptr, "%s", delimiter);
        fprintf(theta_fptr, format, soil_theta[3]);
        fprintf(theta_fptr, "\n");
    }
}

// Write volume balance summary using BMI-retrieved data
static void write_volume_balance_summary_bmi(FILE* output_fptr, const char* config_file,
                                            const char* forcing_file, int timesteps_run,
                                            Bmi* model, int simulate_discrete_soil_moisture,
                                            double calc_vol_balance_residual) {
    if (output_fptr == NULL || model == NULL) return;
    
    // Get volume balance residual from BMI
    double vol_balance_residual = 0.0;
    model->get_value(model, "vol_balance_residual_m", &vol_balance_residual);
    
    // Get model parameters for context
    double catchment_area_km2 = 0.0;
    double soil_depth_m = 0.0;
    double soil_porosity = 0.0;
    
    model->get_value(model, "param_catchment_area_km2", &catchment_area_km2);
    model->get_value(model, "param_soil_depth_m", &soil_depth_m);
    model->get_value(model, "param_soil_porosity", &soil_porosity);
    
    // Get current time step information
    double time_step_seconds = 0.0;
    model->get_time_step(model, &time_step_seconds);
    
    // Write volume balance summary
    fprintf(output_fptr, "# CFE BMI Driver Volume Balance Summary\n");
    fprintf(output_fptr, "# Generated by: CFE BMI Driver v3.0\n");
    fprintf(output_fptr, "# ==========================================\n");
    fprintf(output_fptr, "# Configuration Information:\n");
    fprintf(output_fptr, "# Config file:              %s\n", config_file);
    fprintf(output_fptr, "# Forcing file:             %s\n", forcing_file);
    fprintf(output_fptr, "# Timesteps completed:      %d\n", timesteps_run);
    fprintf(output_fptr, "# Time step size:           %.0f seconds\n", time_step_seconds);
    fprintf(output_fptr, "# Soil moisture model:      %s\n", 
            simulate_discrete_soil_moisture ? "Discrete soil moisture balance" : "Conceptual linear reservoir");
    fprintf(output_fptr, "# ==========================================\n");
    fprintf(output_fptr, "# Model Parameters:\n");
    fprintf(output_fptr, "# Catchment area:           %.6f km2\n", catchment_area_km2);
    fprintf(output_fptr, "# Soil depth:               %.3f m\n", soil_depth_m);
    fprintf(output_fptr, "# Soil porosity:            %.3f [-]\n", soil_porosity);
    fprintf(output_fptr, "# ==========================================\n");
    fprintf(output_fptr, "# Volume Balance Results:\n");
    fprintf(output_fptr, "# Global volume balance residual (bmi): %.12e m\n", vol_balance_residual);
    fprintf(output_fptr, "# Calculated vol. balance residual:     %.12e m\n", calc_vol_balance_residual);
    
    // Interpret the residual
    if (fabs(vol_balance_residual) < 1.0e-12) {
        fprintf(output_fptr, "# Volume balance check: PASSED (residual < 1e-12)\n");
    } else if (fabs(vol_balance_residual) < 1.0e-10) {
        fprintf(output_fptr, "# Volume balance check: ACCEPTABLE (residual < 1e-10)\n");
    } else {
        fprintf(output_fptr, "# Volume balance check: WARNING (residual >= 1e-10)\n");
    }
    
    fprintf(output_fptr, "# =================================================================\n");
    fprintf(output_fptr, "# Note: For some reason complete volume balance details unavailable\n");
    fprintf(output_fptr, "#       through BMI get_value() calls for individual volume balance\n");
    fprintf(output_fptr, "#       components.  It seems volbal struct is not persistent.\n");
    fprintf(output_fptr, "#       Calculated vol. balance residual based on sum of individual\n");
    fprintf(output_fptr, "#       get_value() calls after each bmi call to cfe_step() (update).\n");
    fprintf(output_fptr, "# ==================================================================\n");
}
