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

#ifndef CFE_DRIVER_UTILS_H
#define CFE_DRIVER_UTILS_H

#include <stdio.h>
#include "cfe_config.h"
#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Command line argument parsing
//------------------------------
void print_usage(const char* prog);

int parse_command_line(int argc, char* argv[], cfe_cmdline_args_struct* args);

int check_command_line_overrides(cfe_options_struct* options, const cfe_cmdline_args_struct* args);

// File I/O operations
//--------------------
int open_output_files(const cfe_options_struct* options, const cfe_cmdline_args_struct* args,
                     FILE** q_fptr, FILE** Q_fptr, FILE** fluxes_fptr, FILE** storages_fptr, 
                     FILE** volbal_fptr, FILE** theta_fptr);
                     
void close_output_files(FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr, FILE* storages_fptr, 
                       FILE* volbal_fptr, FILE* theta_fptr, FILE* warnings_fptr, 
                       FILE* forcing_fptr);

// Forcing file operations
//------------------------
int open_forcing_file(const char* path, FILE** fptr, aorc_cols_t* cols);

int read_next_forcing_aorc(FILE* f, const aorc_cols_t* cols, int dt_seconds,
                          cfe_forcing_struct* forcing, char* time_str);
                          
int parse_time_string(const char* time_str, aorc_forcing_time_struct* forcing_time);

int calculate_time_delta(const char* time1_str, const char* time2_str);

int analyze_forcing_file(const char* forcing_filename, aorc_forcing_time_struct* forcing_time, int verbosity);

// Output formatting and writing
//------------------------------
void format_timestamp(char* timestamp_str, size_t str_size, int timestep,
                     const aorc_forcing_time_struct* forcing_time,
                     const char* time_format);
                     
void write_output_headers(const cfe_options_struct* options,
                         const char* time_format,
                         const char* delimiter,
                         FILE* q_fptr,
                         FILE* Q_fptr,
                         FILE* fluxes_fptr,
                         FILE* storages_fptr,
                         FILE* theta_fptr);
                         
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
                      FILE* theta_fptr);

// Volume balance operations
//--------------------------
void cfe_initialize_volume_balance(const cfe_parameters_struct* params,
                                  const cfe_options_struct* options, 
                                  const cfe_state_struct* state,
                                  cfe_volbal_struct* volbal);
                                  
void write_volume_balance_summary(FILE* output_fptr,
                                 const cfe_options_struct* options,
                                 const cfe_volbal_struct* volbal,
                                 const cfe_state_struct* final_state,
                                 const cfe_parameters_struct* params,
                                 double *soil_dz);

// Configuration file writing
//---------------------------
int write_hotstart_config(const CFE_CONFIG* cfg,
                         const cfe_state_struct* state,
                         const aorc_forcing_time_struct* forcing_time);

// Time utilities
//---------------
int calculate_julian_day(int year, int month, int day);

double calculate_julian_date(int year, int month, int day, int hour, int minute, int second);

#ifdef __cplusplus
}
#endif

#endif /* CFE_DRIVER_UTILS_H */
