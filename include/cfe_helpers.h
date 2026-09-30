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
 

#ifndef CFE_HELPERS_H
#define CFE_HELPERS_H

#include "cfe_types.h"
#include "cfe_config.h"
#include "cfe.h"
#include "cfe_soil_discrete.h"

/* Unit conversions between config/BMI units and internal SI */
static inline double cm_to_m(double x)              { return x / 100.0; }
static inline double m_to_cm(double x)              { return x * 100.0; }
static inline double cm_per_h_to_m_per_s(double x)  { return x / 360000.0; }
static inline double m_per_s_to_cm_per_h(double x)  { return x * 360000.0; }

/* Defaults */

int is_leap_year(int year);
int calculate_day_of_year(int year, int month, int day);
int day_of_year_from_epoch(double epoch_seconds);
double parse_date_to_epoch(const char *date_str);

int is_fabs_less_than_epsilon(double a,double epsilon);

void set_parameters_defaults(cfe_parameters_struct* params);

void set_options_defaults(cfe_options_struct* opts);

void set_state_defaults(cfe_state_struct* state);        // sets all state variable defaults

int normalize_config_units(CFE_CONFIG* config);

int map_config_to_cfe_structs(CFE_CONFIG* config, cfe_parameters_struct* params, cfe_options_struct* options);

/* Validate cfe required parameters (errors -> -1, warnings printed to stderr) */
int validate_required_parameters(const CFE_CONFIG* cfg, const int verbosity);   // all versions

/* Map config struct elements to parameters, options, and states (intelligently copy array elements -hopefully) */
int map_config_to_parameters_and_options(const CFE_CONFIG* cfg,              // all versions
                                             cfe_parameters_struct* params,
                                             cfe_options_struct* opts);

/* Parse, validate, and map a v3 config file into CFE model structs */
int parse_config_driver(const char* config_file, double config_file_version,
                             CFE_CONFIG* config,
                             cfe_parameters_struct* params,
                             cfe_options_struct* opts);

/* Initialize state from params/options, allocating queues as needed */
int cfe_initialize(const cfe_parameters_struct* params,
                   const cfe_options_struct* opts,
                   cfe_state_struct* state);

/* Run one step. TODO: wire to cfe() kernel inside this function. */
int cfe_step(const cfe_parameters_struct* params,
             const cfe_options_struct* opts,
             cfe_state_struct* state,
             const cfe_forcing_struct* forcing,
             double dt_seconds,
             cfe_outputs_struct* outputs,
             cfe_volbal_struct* volbal);

                                   
double cfe_get_last_qout_m(const cfe_outputs_struct* outputs);

// ADDED FOR DEBUGGING PURPOSES:
void print_cfe_input_debug(const cfe_options_struct*    o,
                           const cfe_parameters_struct* p,
                           const cfe_state_struct*      s,
                           const char*                  cfg_path,
                           const char*                  forcing_path,
                           const char*                  qout_path,
                           const char*                  volbal_path);

// ADDED FOR DEBUGGING PURPOSES:
void print_exchange_values(int timestep,
                          const cfe_parameters_struct* p,
                          const cfe_options_struct* o,
                          const cfe_state_struct* s,
                          const cfe_forcing_struct* forcing,
                          double dt_seconds,
                          const cfe_outputs_struct* outputs);

// ADDED FOR DEBUGGING PURPOSES:
void check_dsbm_local_volume_balance(
    TimestepSoilVolbal *soil_volbal,
    SoilFluxes *soil_fluxes,
    SoilStateIn *soil_state_in,
    SoilStateOut *soil_state_out,
    SoilGeometry *soil_geometry,
    double infiltration_depth_m,
    double actual_et_from_soil_m,
    double balance_tolerance
);

/* Output writing helper functions */
const char* get_delimiter_string(const char* delimiter_name);

/* Re-derive catchment_bare_soil_fraction from vegetated + impervious fractions.
 * Call after any BMI set_value on catchment_vegetated_fraction or impervious. */
int update_catchment_land_cover_fractions(cfe_parameters_struct* parameters);

/* Resync derived quantities and cached copies after calibration parameter changes.
 * Call once before the first update after any set_value on a calibration parameter. */
int cfe_resync_derived_params(cfe_parameters_struct* params,
                              const cfe_options_struct* opts,
                              cfe_state_struct* state);

/* Cleanup any allocations made in state */
int cfe_finalize(cfe_state_struct* s);


#endif
