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

// THIS source file contains a bunch of helper functions that support the CFE model.  Some very important 
// things happen in the code in this file including:

//  1. defaultt parameter value assignment
//  2. parameter validation
//  3. parameter mapping from the config structt that is read from the config files into the CFE model
//     options, parameters, state, and output structs.
//  4. the cfe_step() function is called by the update() bmi wrapper function.
//  FLO 9/2025

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "cfe_helpers.h"
#include "parser_helpers.h"
#include "cfe.h"  // brings in NWM_SOIL_PARAMETERS_STRUCTURE, CONCEPTUAL_RESERVOIR_STRUCTURE, NASH_CASCADE_PARAMETERS_STRUCTURE, rainfall_partitioning_parameters_structure, evapotranspiration_structure, volbal_struct
#include "nash_cascade.h"  // MAX_NUM_SURFACE_NASH_CASCADE, MAX_NUM_SUBSURFACE_NASH_CASCADE
#include "giuh.h"          // MAX_NUM_GIUH_ORDINATES
#include "soil_helpers.h"  // Brings in codes needed for discrete soil moisture simulation
#include "soil_config.h"   // THETA_MIN

int is_leap_year(int year)
{
    return ((year % 4 == 0 && year % 100 != 0) ||
            (year % 400 == 0));
}

int calculate_day_of_year(int year, int month, int day)
{
    static const int days_before_month[12] = {
        0, 31, 59, 90, 120, 151,
        181, 212, 243, 273, 304, 334
    };

    if (month < 1 || month > 12 || day < 1) {
        return 1;
    }

    int day_of_year = days_before_month[month - 1] + day;

    if (month > 2 && is_leap_year(year)) {
        day_of_year++;
    }

    return day_of_year;
}

int day_of_year_from_epoch(double epoch_seconds)
{
    time_t t = (time_t)epoch_seconds;
    struct tm utc;
    if (gmtime_r(&t, &utc) == NULL)
        return 1;
    return utc.tm_yday + 1;
}

double parse_date_to_epoch(const char *date_str)
{
    if (date_str == NULL || date_str[0] == '\0')
        return 0.0;
    int y, m, d;
    if (sscanf(date_str, "%d-%d-%d", &y, &m, &d) != 3)
        return 0.0;
    struct tm tm_val = {0};
    tm_val.tm_year = y - 1900;
    tm_val.tm_mon  = m - 1;
    tm_val.tm_mday = d;
    time_t t = timegm(&tm_val);
    if (t == (time_t)-1)
        return 0.0;
    return (double)t;
}

// Helper functions

static void trim_inplace(char* s)
{
    if (s == NULL) return;
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) {
        s[--n] = '\0';
    }
    size_t i = 0;
    while (s[i] && isspace((unsigned char)s[i])) i++;
    if (i > 0) memmove(s, s + i, n - i + 1);
}

static int equals_ic(const char* a, const char* b)
{
    return string_compare_ignore_case(a, b) == 0;
}

// is a parameter p out of bounds? 
int is_oob(double p, double lower_lim, double upper_lim) {
    return (p < lower_lim || p > upper_lim);
}

int is_fabs_less_than_epsilon(double a,double epsilon)  // returns true if fabs(a)<epsilon
{
  if(fabs(a)<epsilon) return(TRUE);
  else                return(FALSE);
}

// ----- scheme name helper functions
//##########################
static const char* surf_scheme_name(cfe_surface_route_t scheme)
{
    switch (scheme) {
        case SURF_ROUTE_GIUH:         return "GIUH";
        case SURF_ROUTE_NASH_CASCADE: return "Nash Cascade";
        default:                      return "Unknown";
    }
}
static const char* partition_scheme_name(cfe_partition_scheme_t scheme)
{
    switch (scheme) {
        case PARTITION_SCHAAKE:     return "Schaake";
        case PARTITION_XINANJIANG:  return "Xinanjiang";
        default:                    return "Unknown";
    }
}

// ---------------- Defaults ----------------
//##############################
/* NOTE: cfe_parameters_struct uses fixed-size arrays only (no pointers).
 * This enables safe memset/calloc initialization. If pointers are ever
 * added, this function and cfe_finalize must be updated accordingly. */
void set_parameters_defaults(cfe_parameters_struct* p)
{
    if (p == NULL) return;
    memset(p, 0, sizeof(*p));
    p->soil_depth_m = 2.0;
    p->nash_subsurface_N = 2;

    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++)
        p->nash_subsurface_init_storage_m[i] = 0.0;

    for (int i = 0; i < MAX_NUM_GIUH_ORDINATES; i++) {
        p->giuh_ordinates[i]   = 0.0;
        p->giuh_init_queue_m[i] = 0.0;
    }

}

//###########################
void set_options_defaults(cfe_options_struct* o)
{
    if (o == NULL) return;
    memset(o, 0, sizeof(*o));
    o->liquid_partitioning_scheme = PARTITION_SCHAAKE;
    o->surface_routing_scheme = SURF_ROUTE_GIUH;
    o->time_step_seconds = 3600.0;
}
//#####################
void set_state_defaults(cfe_state_struct* s)
{
    if (s == NULL) return;
    memset(s, 0, sizeof(*s));
}

// Sanity checks on parsed config values
//########################
int normalize_config_units(CFE_CONFIG* config) {

    // Basic sanity checks on critical values
    if (config->soil_depth_m <= 0.0) {
        fprintf(stderr, "ERROR: soil_depth_m must be > 0 after parsing\n");
        return -1;
    }

    if (config->soil_effective_porosity <= 0.0 || config->soil_effective_porosity > 1.0) {
        fprintf(stderr, "ERROR: soil_effective_porosity must be in (0,1] after parsing\n");
        return -1;
    }

    return 0;
}

// ---------------- Validation (strict + warnings) ----------------
// Unified parameter validation ffor all CFE versions
// Uses human-interpretable units (cm/h, cm, etc.) as stored in the temporary config structt
// Converts from human-interpretable units m, m/s before stuffing into CFE model structs
//###############################
int validate_required_parameters(const CFE_CONFIG* cfg, const int verbosity) {
    
    if(verbosity > 0) fprintf(stderr, "Validating CFE version %.1f config parameters.\n", cfg->version);
    
    //============================
    // Parameter bounds (human-readable units)   FLO bounds with input from Xia Feng
    //============================
    const double Ksat_LOW_cm_per_h      = 0.01;   // 0.07 used a lower limit in NWM calibration
    const double Ksat_HIGH_cm_per_h     = 510.0;  // ridiculously high, but used in NWM calibration     
    const double satpsi_LOW_cm          = 0.5;    // 3.6 used in NWM calibration  
    const double satpsi_HIGH_cm         = 200.0;  // 95.5 used in NWM calibration
    const double claphorn_b_LOW         = 2.0;    // limit used in NWM calibrations
    const double claphorn_b_HIGH        = 15.0;   // limit used in NWM calibrations       
    const double alpha_fc_LOW           = 0.14;   // 0.15 typ. ffor sands
    const double alpha_fc_HIGH          = 0.35;   // 0.36 typ. ffor other soil textures
    const double porosity_LOW           = 0.05;   // 0.16 used in NWM calibration
    const double porosity_HIGH          = 0.7;    // 0.58 used in NWM calibration 0.7 is pretty darn high (volcanic clays)
    const double gw_max_storage_LOW_m   = 0.01;   // 0.01 used in NWM calibration
    const double gw_max_storage_HIGH_m  = 3.0;    // 0.25 used in NWM calibration        
    const double gw_init_LOW_m          = 0.0;
    const double cgw_LOW_per_ts         = 1.8e-06; // used in NWM calibration
    const double cgw_HIGH_per_ts        = 1.8e-03; // used in NWM calibration
    const double gw_exp_LOW             = 1.0;    // used in NWM calibration
    const double gw_exp_HIGH            = 8.0;    // used in NWM calibration
    const double k_lf_LOW               = 0.0;    // should never be negative
    const double k_lf_HIGH              = 1.0;    // should never be > 1
    const double nashK_subsurf_LOW_per_ts       = 0.0;
    const double nashK_subsurf_HIGH_per_ts      = 1.0;
    const double timestep_LOW_h         = 0.25;     
    const double timestep_HIGH_h        = 24.0;       
    
    const int has_discrete_moisture = cfg->control_soil_simulate_discrete_soil_moisture_true_false;
    
    // Derived constraints
    const double max_soil_storage_m = cfg->soil_effective_porosity * cfg->soil_depth_m;
    
    //============================
    // Core parameter validation
    //============================
    
    // Control parameters
    // NOTE: When running under BMI (e.g. ngen), the forcing filename is typically
    // set to "BMI" by the config file — forcings arrive via set_value(), not from
    // this file path. The standalone driver uses this path to open a CSV forcing file.
    if (strlen(cfg->control_input_forcing_filename) == 0) {
        fprintf(stderr, "ERROR: Missing forcing filename\n");
        return -1;
    }
    
    if (is_oob(cfg->timestep_h, timestep_LOW_h, timestep_HIGH_h)) {
        fprintf(stderr, "ERROR: timestep_h %.3e out of bounds [%.3e, %.3e]\n",
                cfg->timestep_h, timestep_LOW_h, timestep_HIGH_h);
        return -1;
    }
    
    if (cfg->total_timesteps < 0) {
        fprintf(stderr, "ERROR: total_timesteps %d must be >= 0\n", cfg->total_timesteps);
        return -1;
    }
    
    // Soil parameters
    if (cfg->soil_depth_m <= 0.0) {
        fprintf(stderr, "ERROR: soil_depth_m must be > 0\n");
        return -1;
    }
    
    if (is_oob(cfg->soil_effective_porosity, porosity_LOW, porosity_HIGH)) {
        fprintf(stderr, "ERROR: soil_effective_porosity %.3e out of bounds [%.3e, %.3e]\n",
                cfg->soil_effective_porosity, porosity_LOW, porosity_HIGH);
        return -1;
    } 
    
    if (is_oob(cfg->soil_Clapp_Hornberger_exponent_b, claphorn_b_LOW, claphorn_b_HIGH)) {
        fprintf(stderr, "ERROR: soil_Clapp_Hornberger_exponent_b %.3e out of bounds [%.3e, %.3e]\n",
                cfg->soil_Clapp_Hornberger_exponent_b, claphorn_b_LOW, claphorn_b_HIGH);
        return -1;
    }
    
    // Hydraulic conductivity (both parsers normalize to cm/h)
    if (is_oob(cfg->soil_sat_hydraulic_conductivity_cm_per_h, Ksat_LOW_cm_per_h, Ksat_HIGH_cm_per_h)) {
        fprintf(stderr, "ERROR: soil_sat_hydraulic_conductivity_cm_per_h %.3e out of bounds [%.3e, %.3e]\n",
                cfg->soil_sat_hydraulic_conductivity_cm_per_h, Ksat_LOW_cm_per_h, Ksat_HIGH_cm_per_h);
        return -1;
    }
    
    // Capillary head (expected in cm in config)
    if (cfg->soil_sat_capillary_head_cm > 0) {
        if (is_oob(cfg->soil_sat_capillary_head_cm, satpsi_LOW_cm, satpsi_HIGH_cm)) {
            fprintf(stderr, "ERROR: soil_sat_capillary_head_cm %.3e out of bounds [%.3e, %.3e]\n",
                    cfg->soil_sat_capillary_head_cm, satpsi_LOW_cm, satpsi_HIGH_cm);
            return -1;
        }
    }
    
    if (is_oob(cfg->soil_field_capacity_Pcap_over_Patm_0_1, alpha_fc_LOW, alpha_fc_HIGH)) {
        fprintf(stderr, "ERROR: soil_field_capacity_Pcap_over_Patm_0_1 %.3e out of bounds [%.3e, %.3e]\n",
                cfg->soil_field_capacity_Pcap_over_Patm_0_1, alpha_fc_LOW, alpha_fc_HIGH);
        return -1;
    }
    
    if (is_oob(cfg->soil_to_gw_percolation_rate_limiter_0_to_1, 0.0, 1.0)) {
        fprintf(stderr, "ERROR: soil_to_gw_percolation_rate_limiter_0_to_1 %.3e out of bounds [0.0, 1.0]\n",
                cfg->soil_to_gw_percolation_rate_limiter_0_to_1);
        return -1;
    }
    
    //============================
    // Soil storage validation (version-specific)
    //============================
    if (has_discrete_moisture) {
        // v3.0+ discrete moisture validation
        int has_discrete_theta = TRUE;
        for (int i = 0; i < NDISCS; i++) {
            if (cfg->soil_reservoir_init_discrete_storage_theta[i] < THETA_MIN) {
                has_discrete_theta = FALSE;
                break;
            }
        }

        if (cfg->control_soil_use_lookup_table_num_points < 0) {
            fprintf(stderr, "ERROR: control_soil_use_lookup_table_num_points must be >= 0\n");
            return -1;
        }
        if (cfg->control_soil_use_lookup_table_num_points > MAX_LOOKUP_TABLE_POINTS) {
            fprintf(stderr, "ERROR: control_soil_use_lookup_table_num_points entered as: %d but cannot exceed %d\n", 
                             cfg->control_soil_use_lookup_table_num_points, MAX_LOOKUP_TABLE_POINTS);
            return -1;
        }

        
        int has_total_storage = (cfg->soil_reservoir_init_storage_m > 0.0);
        
        if (!has_discrete_theta && !has_total_storage) {
            fprintf(stderr, "ERROR: When simulating discretized dynamic soil moisture, specify either:\n");
            fprintf(stderr, "  - soil_reservoir_init_storage_m, OR\n");
            fprintf(stderr, "  - soil_reservoir_init_discrete_storage_theta (all %d values)\n", NDISCS);
            return -1;
        }
        
        if (has_discrete_theta) {
            for (int i = 0; i < NDISCS; i++) {
                double theta = cfg->soil_reservoir_init_discrete_storage_theta[i];
                if (is_oob(theta, THETA_MIN, cfg->soil_effective_porosity)) {
                    fprintf(stderr, "ERROR: soil_reservoir_init_discrete_storage_theta[%d] = %.4f\n", i, theta);
                    fprintf(stderr, "       Must be between %.4f and soil_effective_porosity (%.4f)\n", 
                            THETA_MIN, cfg->soil_effective_porosity);
                    return -1;
                }
            }
        }
        
        if (has_total_storage) {
            double min_storage = (double)THETA_MIN * cfg->soil_depth_m;
            if (is_oob(cfg->soil_reservoir_init_storage_m, min_storage, max_soil_storage_m)) {
                fprintf(stderr, "ERROR: soil_reservoir_init_storage_m %.6f m out of bounds [%.6f, %.6f]\n", 
                        cfg->soil_reservoir_init_storage_m, min_storage, max_soil_storage_m);
                return -1;
            }
        }
    } else {
        // Standard soil storage (absolute storage in meters)
        if (is_oob(cfg->soil_reservoir_init_storage_m, 0.0, max_soil_storage_m)) {
            fprintf(stderr, "ERROR: soil_reservoir_init_storage_m %.3e m out of bounds [0, %.3e]\n",
                    cfg->soil_reservoir_init_storage_m, max_soil_storage_m);
            return -1;
        }
    }
    
    //============================
    // Groundwater validation
    //============================
    if (is_oob(cfg->gw_reservoir_max_storage_m, gw_max_storage_LOW_m, gw_max_storage_HIGH_m)) {
        fprintf(stderr, "ERROR: gw_reservoir_max_storage_m %.3e out of bounds [%.3e, %.3e]\n",
                cfg->gw_reservoir_max_storage_m, gw_max_storage_LOW_m, gw_max_storage_HIGH_m);
        return -1;
    }
    
    if (is_oob(cfg->gw_reservoir_init_storage_m, gw_init_LOW_m, cfg->gw_reservoir_max_storage_m)) {
        fprintf(stderr, "ERROR: gw_reservoir_init_storage_m %.3e must be in [%.3e, %.3e]\n",
                cfg->gw_reservoir_init_storage_m, gw_init_LOW_m, cfg->gw_reservoir_max_storage_m);
        return -1;
    }
    
    if (is_oob(cfg->gw_discharge_coeff_m_per_timestep, cgw_LOW_per_ts, cgw_HIGH_per_ts)) {
        fprintf(stderr, "ERROR: gw_discharge_coeff_m_per_timestep %.3e out of bounds [%.3e, %.3e]\n",
                cfg->gw_discharge_coeff_m_per_timestep, cgw_LOW_per_ts, cgw_HIGH_per_ts);
        return -1;
    }
    
    if (is_oob(cfg->gw_discharge_exponent, gw_exp_LOW, gw_exp_HIGH)) {
        fprintf(stderr, "ERROR: gw_discharge_exponent %.3e out of bounds [%.3e, %.3e]\n",
                cfg->gw_discharge_exponent, gw_exp_LOW, gw_exp_HIGH);
        return -1;
    }
    
    //============================
    // Additional parameters (if present)
    //============================
    if (cfg->soil_reservoir_rate_const_to_subsurface_lateral_flow >= 0) {
        if (is_oob(cfg->soil_reservoir_rate_const_to_subsurface_lateral_flow, k_lf_LOW, k_lf_HIGH)) {
            fprintf(stderr, "ERROR: soil_reservoir_rate_const_to_subsurface_lateral_flow %.3e out of bounds [%.3e, %.3e]\n",
                    cfg->soil_reservoir_rate_const_to_subsurface_lateral_flow, k_lf_LOW, k_lf_HIGH);
            return -1;
        }
    }
    
    if (cfg->subsurface_routing_nash_K > 0) {
        if (is_oob(cfg->subsurface_routing_nash_K, nashK_subsurf_LOW_per_ts, nashK_subsurf_HIGH_per_ts)) {
            fprintf(stderr, "ERROR: subsurface_routing_nash_K %.3e out of bounds [%.3e, %.3e]\n",
                    cfg->subsurface_routing_nash_K, nashK_subsurf_LOW_per_ts, nashK_subsurf_HIGH_per_ts);
            return -1;
        }
    }
    
    //============================
    // Scheme-specific validation
    //============================

    // Partitioning scheme (exactly one)
    int is_schaake = (string_compare_ignore_case(cfg->partitioning_scheme_name, "schaake") == 0);
    int is_xinan   = (string_compare_ignore_case(cfg->partitioning_scheme_name, "xinanjiang") == 0);
    if (!(is_schaake ^ is_xinan)) {
        fprintf(stderr, "ERROR: Must specify exactly one partitioning scheme: 'schaake' or 'xinanjiang'\n");
        return -1;
    }

    // Xinanjiang parameters (if selected)
    if (is_xinan) {
        if (is_oob(cfg->soil_Xinanjiang_tension_water_inflection_point, -0.49, 0.49)) {
            fprintf(stderr, "ERROR: soil_Xinanjiang_tension_water_inflection_point %.3e out of bounds [-0.49, 0.49]\n",
                    cfg->soil_Xinanjiang_tension_water_inflection_point);
            return -1;
        }
        if (is_oob(cfg->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent, 0.0, 1.0)) {
            fprintf(stderr, "ERROR: soil_Xinanjiang_tension_water_soil_moist_distrib_exponent %.3e out of bounds [0.0, 1.0]\n",
                    cfg->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent);
            return -1;
        }
        if (is_oob(cfg->soil_Xinanjiang_free_water_soil_moist_distrib_exponent, 0.0, 1.0)) {
            fprintf(stderr, "ERROR: soil_Xinanjiang_free_water_soil_moist_distrib_exponent %.3e out of bounds [0.0, 1.0]\n",
                    cfg->soil_Xinanjiang_free_water_soil_moist_distrib_exponent);
            return -1;
        }
    }

    // GIUH validation - NEW: Added bounds checking ffor ordinates count
    
    if (cfg->surface_routing_num_giuh_ordinates < 1) {
        fprintf(stderr, "ERROR: surface_routing_num_giuh_ordinates must be >= 1\n");
        return -1;
    }
    if (cfg->surface_routing_num_giuh_ordinates > MAX_NUM_GIUH_ORDINATES) {
        fprintf(stderr, "ERROR: surface_routing_num_giuh_ordinates cannot exceed %d\n", MAX_NUM_GIUH_ORDINATES);
        return -1;
    }
    
    // Deepest root zone discretization bounds
    if (has_discrete_moisture) {
        if (cfg->control_ET_deepest_root_zone_discretization < 1 ||
            cfg->control_ET_deepest_root_zone_discretization > NDISC) {
            fprintf(stderr, "ERROR: control_ET_deepest_root_zone_discretization must be 1..%d (got %d)\n",
                    NDISC, cfg->control_ET_deepest_root_zone_discretization);
            return -1;
        }
    }

    //============================
    // Warnings ffor non-critical issues
    //============================
    if (cfg->soil_depth_m != 2.0) {
        fprintf(stderr, "WARN: soil_depth_m is %.6f (expected 2.0 m)\n", cfg->soil_depth_m);
    }
    if (cfg->soil_effective_porosity > 0.6) {
        fprintf(stderr, "WARN: soil_effective_porosity = %.6f (> 0.6 is unusually high)\n",
                cfg->soil_effective_porosity);
    }
    if (cfg->soil_Clapp_Hornberger_exponent_b < 4.0 || cfg->soil_Clapp_Hornberger_exponent_b > 12.0) {
        fprintf(stderr, "WARN: soil_Clapp_Hornberger_exponent_b = %.6f (typical range ~4..12)\n",
                cfg->soil_Clapp_Hornberger_exponent_b);
    }
    
    return 0;  // Success
}  // <---------------------------- End of parameter validation code




int update_catchment_land_cover_fractions(cfe_parameters_struct* parameters)
{
    double land_cover_fraction_sum;

    if (parameters == NULL) {
        fprintf(stderr,
                "ERROR: NULL parameters passed to "
                "update_catchment_land_cover_fractions().\n");
        return -1;
    }

    if (parameters->catchment_impervious_fraction < 0.0 ||
        parameters->catchment_impervious_fraction > 1.0) {
        fprintf(stderr,
                "ERROR: catchment_impervious_fraction must be between "
                "0 and 1; got %.17g\n",
                parameters->catchment_impervious_fraction);
        return -1;
    }

    if (parameters->catchment_vegetated_fraction < 0.0 ||
        parameters->catchment_vegetated_fraction > 1.0) {
        fprintf(stderr,
                "ERROR: catchment_vegetated_fraction must be between "
                "0 and 1; got %.17g\n",
                parameters->catchment_vegetated_fraction);
        return -1;
    }

    parameters->catchment_bare_soil_fraction = 0.0;

    land_cover_fraction_sum =
        parameters->catchment_impervious_fraction +
        parameters->catchment_vegetated_fraction;

    if (land_cover_fraction_sum <= 1.0) {
        parameters->catchment_bare_soil_fraction =
            1.0 - land_cover_fraction_sum;
    }
    else {
        fprintf(stderr,
                "WARNING: catchment_vegetated_fraction (%.6f) + "
                "catchment_impervious_fraction (%.6f) exceeds 1.0. "
                "Reducing vegetated fraction to %.6f and setting "
                "bare_soil_fraction to 0.0.\n",
                parameters->catchment_vegetated_fraction,
                parameters->catchment_impervious_fraction,
                1.0 -
                    parameters->catchment_impervious_fraction);

        parameters->catchment_vegetated_fraction =
            1.0 - parameters->catchment_impervious_fraction;

        parameters->catchment_bare_soil_fraction = 0.0;
    }

    return 0;
}

// ---------------- Mapper ----------------
// This function maps values from the cfe config struct (cfg) into stateless cfe arrays
// Map config structt to CFE model structs with final unit conversions
// Convert config entries from human-interpretable config units (cm/h, cm) to CFE model units (m, s)
//######################################
int map_config_to_parameters_and_options(const CFE_CONFIG* cfg,
                                             cfe_parameters_struct* p,
                                             cfe_options_struct* o)
{
    if (cfg == NULL || p == NULL || o == NULL) return -1;
    
    set_parameters_defaults(p);
    set_options_defaults(o);

    o->cfe_version = cfg->version;
    // options
    //-----------------------------------------------------
    // Partitioning
    char part[64];
    strncpy(part, cfg->partitioning_scheme_name, sizeof(part) - 1);
    part[sizeof(part) - 1] = '\0';
    trim_inplace(part);

    if (equals_ic(part, "Xinanjiang")) {
        o->liquid_partitioning_scheme = PARTITION_XINANJIANG;
    } else if (equals_ic(part, "Schaake")) {
        o->liquid_partitioning_scheme = PARTITION_SCHAAKE;
    } else {
        fprintf(stderr, "ERROR: Unrecognized partitioning scheme '%s' (expected 'Schaake' or 'Xinanjiang')\n", part);
        return -1;
    }

    // Surface routing
    char surf[64];
    strncpy(surf, cfg->surface_routing_scheme_name, sizeof(surf) - 1);
    surf[sizeof(surf) - 1] = '\0';
    trim_inplace(surf);

//
//
//    if (equals_ic(surf, "NASH_CASCADE") || equals_ic(surf, "NASHCASCADE") || equals_ic(surf, "NASH")) {
//        o->surface_routing_scheme = SURF_ROUTE_NASH_CASCADE;
//    } else if (equals_ic(surf, "GIUH")) {
//        o->surface_routing_scheme = SURF_ROUTE_GIUH;
//    } else {
//        o->surface_routing_scheme = SURF_ROUTE_GIUH;  // default
//    }

    o->surface_routing_scheme = SURF_ROUTE_GIUH;  // default

    //  DEBUG:
    if(o->verbosity > 0) fprintf(stderr,"In map_config_to_parameters_and_options(): o->surface_routing_scheme=%d\n", 
                                 o->surface_routing_scheme);

    o->time_step_seconds               = (int)(cfg->timestep_h * 3600.0 + 0.5);
    o->num_timesteps                   = cfg->total_timesteps;
    o->verbosity                       = cfg->verbosity;
    o->epoch_start_seconds             = parse_date_to_epoch(cfg->simulation_start_date);
    o->enable_ET_Priestley_Taylor      = (cfg->et_alpha_pt > 1.0e-03) ? TRUE : FALSE;   // iff alpha_pt not zero or tiny.
    o->enable_freeze_thaw              = cfg->control_soil_simulate_freeze_thaw_true_false;
    o->simulate_discrete_soil_moisture = cfg->control_soil_simulate_discrete_soil_moisture_true_false;
    o->simulate_soil_evaporation       = cfg->control_soil_simulate_soil_evaporation;
    o->deepest_root_zone_disc          = cfg->control_ET_deepest_root_zone_discretization;
    o->use_soil_lookup_table           = (cfg->control_soil_use_lookup_table_num_points > 0) ? TRUE : FALSE;
    if(o->verbosity > 1) printf("DEBUG: Mapped use_soil_lookup_table = %d\n", o->use_soil_lookup_table);    
    
    strncpy(o->input_forcing_filename, cfg->control_input_forcing_filename, sizeof(o->input_forcing_filename) - 1);
    o->input_forcing_filename[sizeof(o->input_forcing_filename) - 1] = '\0';
    
    // ET using Priestley-Taylor method
    if(o->enable_ET_Priestley_Taylor == TRUE) 
        p->alpha_pt = cfg->et_alpha_pt;
    else
        p->alpha_pt = 0.0; 

    // Catchment
    p->catchment_area_km2 = cfg->cat_area_km2;

    // Soil
    //-----------------------------------------------------
    p->soil_depth_m                            = cfg->soil_depth_m;
    p->soil_b                                  = cfg->soil_Clapp_Hornberger_exponent_b;
    p->ksat_m_per_s                            = cm_per_h_to_m_per_s(cfg->soil_sat_hydraulic_conductivity_cm_per_h);
    p->sat_capillary_head_m                    = cm_to_m(cfg->soil_sat_capillary_head_cm);
    p->effective_porosity                      = cfg->soil_effective_porosity;
    p->field_capacity_Pcap_over_Patm           = cfg->soil_field_capacity_Pcap_over_Patm_0_1;
    p->soil_init_storage_m                     = cfg->soil_reservoir_init_storage_m;
    p->soil_k_lateral_per_h                    = cfg->soil_reservoir_rate_const_to_subsurface_lateral_flow;
    p->soil_ice_imperv_threshold               = cfg->soil_ice_content_impervious_threshold;
    p->refkdt                                  = 3.0;  // it is a CONSTANT in NWM
    p->schaake_magic_constant                  = p->refkdt * p->ksat_m_per_s / 2.0e-06;
    p->soil_to_gw_percolation_rate_limiter_0_1 = cfg->soil_to_gw_percolation_rate_limiter_0_to_1;
    p->soil_init_storage_m                     = cfg->soil_reservoir_init_storage_m;
   
    if (o->simulate_discrete_soil_moisture) {
    p->lut_n_points = cfg->control_soil_use_lookup_table_num_points;
    p->lut_theta_min = LOOKUP_TABLE_THETA_MIN;
        // Copy the initial theta values
        for (int i = 0; i < NDISC; i++) {
            p->soil_discrete_init_theta[i] = cfg->soil_reservoir_init_discrete_storage_theta[i];
        }
        
    }
    // Calculate field capacity moisture content and storage
    double psi_atm_m = STANDARD_ATM_PRESS_Pa /
                       (GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2 * WATER_LIQUID_DENSITY_kg_per_m3);
    double arg = (p->field_capacity_Pcap_over_Patm * psi_atm_m/p->sat_capillary_head_m);
    p->field_capacity_moisture_content = p->effective_porosity * pow(arg, (-1.0/p->soil_b));
    p->field_capacity_storage_m = p->field_capacity_moisture_content * p->soil_depth_m;

    // Wilting point at 15 atm suction via Clapp-Hornberger
    double psi_15atm_m = 15.0 * psi_atm_m;
    p->wilting_point = p->effective_porosity * pow(psi_15atm_m / p->sat_capillary_head_m, -1.0 / p->soil_b);

    // Xinanjiang
    p->xj_tension_inflection_point               = cfg->soil_Xinanjiang_tension_water_inflection_point;
    p->xj_tension_b                            = cfg->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent;
    p->xj_free_b                               = cfg->soil_Xinanjiang_free_water_soil_moist_distrib_exponent;

    // Land-cover fractions
    p->catchment_impervious_fraction = cfg->cat_impervious_fraction;

    if (o->simulate_soil_evaporation == TRUE) {
        p->catchment_vegetated_fraction = cfg->catchment_vegetated_fraction;
    }
    else {
        p->catchment_vegetated_fraction = 1.0;
    }

    if (update_catchment_land_cover_fractions(p) != 0) {
        return -1;
    }

    p->bare_soil_rsurf_exp = CFE_BARE_SOIL_RSURF_EXP;

    // Groundwater
    p->gw_max_storage_m                        = cfg->gw_reservoir_max_storage_m;
    p->gw_init_storage_m                       = cfg->gw_reservoir_init_storage_m;
    p->gw_discharge_coeff_m_per_timestep              = cfg->gw_discharge_coeff_m_per_timestep;
    p->gw_discharge_exponent                   = cfg->gw_discharge_exponent;

    // Surface routing
    //-----------------------------------------------------
    double giuh_sum = 0.0;
    if (o->surface_routing_scheme == SURF_ROUTE_GIUH) {
    
        p->giuh_num_ordinates = cfg->surface_routing_num_giuh_ordinates;
        
        if(p->giuh_num_ordinates > MAX_NUM_GIUH_ORDINATES) {
            fprintf(stderr, "ERROR: giuh_num_ordinates (%d) cannot exceed %d\n", 
            p->giuh_num_ordinates, MAX_NUM_GIUH_ORDINATES);
            return -1;
        }
        if (p->giuh_num_ordinates > 0) {
            for (int i = 0; i < p->giuh_num_ordinates; i++) {
                p->giuh_ordinates[i] = cfg->surface_routing_giuh_ordinates[i];  // copy ordinates
                giuh_sum += p->giuh_ordinates[i];
            }
            if(fabs(giuh_sum - 1.0) > 1.0e-6) {
                fprintf(stderr, "WARNING: Sum of %d GIUH ordinates = %.6f, sum should equal 1.0 - normalizing\n",
                         p->giuh_num_ordinates, giuh_sum);
                if (giuh_sum > 0.0) {
                    double lambda = 1.0 / giuh_sum;
                    for (int i = 0; i < p->giuh_num_ordinates; i++) {
                        p->giuh_ordinates[i] *= lambda;
                    }
                    giuh_sum = 1.0;
                }
            }
            for (int i = 0; i < p->giuh_num_ordinates; i++) {  // copy the input convolution queue too
                p->giuh_init_queue_m[i] = cfg->surface_routing_init_giuh_convolution_queue_m[i];
            }
        }
    }

    // subsurface nash
    p->nash_subsurface_N = 2;
    p->nash_subsurface_K_per_h = cfg->subsurface_routing_nash_K; // per hour in kernel; convert in adapter if needed
    if (p->nash_subsurface_N > 0) {
        for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
            p->nash_subsurface_init_storage_m[i] = cfg->subsurface_routing_nash_cascade_init_storage_m[i];
        }
    }
    
    // Output options

    if (cfg->output_discharge_filename[0] != '\0') {
        strncpy(o->output_discharge_filename, cfg->output_discharge_filename, sizeof(o->output_discharge_filename) - 1);
        o->output_discharge_filename[sizeof(o->output_discharge_filename) - 1] = '\0';
    } else {
        o->output_discharge_filename[0] = '\0';
    }
    if (cfg->output_total_discharge_m3_per_sec_filename[0] != '\0') {
        strncpy(o->output_total_discharge_m3_per_sec_filename, 
                cfg->output_total_discharge_m3_per_sec_filename, 
                sizeof(o->output_total_discharge_m3_per_sec_filename) - 1);
        o->output_total_discharge_m3_per_sec_filename[sizeof(o->output_total_discharge_m3_per_sec_filename) - 1] = '\0';
    } else {
        o->output_total_discharge_m3_per_sec_filename[0] = '\0';
    }
    if (cfg->output_status_warnings_filename[0] != '\0') {
        strncpy(o->output_status_warnings_filename, cfg->output_status_warnings_filename, sizeof(o->output_status_warnings_filename) - 1);
        o->output_status_warnings_filename[sizeof(o->output_status_warnings_filename) - 1] = '\0';
    } else {
        o->output_status_warnings_filename[0] = '\0';
    }

    if (cfg->output_internal_fluxes_filename[0] != '\0') {
        strncpy(o->output_internal_fluxes_filename, cfg->output_internal_fluxes_filename, sizeof(o->output_internal_fluxes_filename) - 1);
        o->output_internal_fluxes_filename[sizeof(o->output_internal_fluxes_filename) - 1] = '\0';
    } else {
        o->output_internal_fluxes_filename[0] = '\0';
    }

    if (cfg->output_internal_storages_filename[0] != '\0') {
        strncpy(o->output_internal_storages_filename, cfg->output_internal_storages_filename, sizeof(o->output_internal_storages_filename) - 1);
        o->output_internal_storages_filename[sizeof(o->output_internal_storages_filename) - 1] = '\0';
    } else {
        o->output_internal_storages_filename[0] = '\0';
    }

    if (cfg->output_volume_balance_filename[0] != '\0') {
        strncpy(o->output_volume_balance_filename, cfg->output_volume_balance_filename, sizeof(o->output_volume_balance_filename) - 1);
        o->output_volume_balance_filename[sizeof(o->output_volume_balance_filename) - 1] = '\0';
    } else {
        o->output_volume_balance_filename[0] = '\0';
    }

    if (cfg->output_soil_moisture_theta_filename[0] != '\0') {
        strncpy(o->output_soil_moisture_theta_filename, cfg->output_soil_moisture_theta_filename, sizeof(o->output_soil_moisture_theta_filename) - 1);
        o->output_soil_moisture_theta_filename[sizeof(o->output_soil_moisture_theta_filename) - 1] = '\0';
    } else {
        o->output_soil_moisture_theta_filename[0] = '\0';
    }

    if (cfg->output_new_config_filename[0] != '\0') {
        strncpy(o->output_new_config_filename, cfg->output_new_config_filename, sizeof(o->output_new_config_filename) - 1);
        o->output_new_config_filename[sizeof(o->output_new_config_filename) - 1] = '\0';
    } else {
        o->output_new_config_filename[0] = '\0';
    }

    // Validate and set time format
    if (validate_time_format(cfg->output_time_standard_format)) {
        strncpy(o->output_time_standard_format, cfg->output_time_standard_format, sizeof(o->output_time_standard_format) - 1);
        o->output_time_standard_format[sizeof(o->output_time_standard_format) - 1] = '\0';
    } else {
        if(o->verbosity > 0) printf("Warning: Invalid time format '%s', using 'datetime'\n", cfg->output_time_standard_format);
        strncpy(o->output_time_standard_format, "datetime", sizeof(o->output_time_standard_format) - 1);
        o->output_time_standard_format[sizeof(o->output_time_standard_format) - 1] = '\0';
    }

    // Set delimiter (validation happens when we use get_delimiter_string() later)
    strncpy(o->output_file_delimiter, cfg->output_file_delimiter, sizeof(o->output_file_delimiter) - 1);
    o->output_file_delimiter[sizeof(o->output_file_delimiter) - 1] = '\0';

    // Set output format option (one of either "%.Ne" or "%.Nf", where N is the number of digits of precision following the decimal place
    // Validation occurs when we attempt to use it.  If invalid, we don't use it and defaultt to something reasonable
    strncpy(o->output_value_format, cfg->output_value_format, sizeof(o->output_value_format));
    
    // set output path name  
    if(strlen(cfg->output_path_name) == 0) {
        strncpy(o->output_path_name, "./", sizeof(o->output_path_name) - 1);
    } else {
        strncpy(o->output_path_name, cfg->output_path_name, sizeof(o->output_path_name) - 1);
    }
    o->output_path_name[sizeof(o->output_path_name) - 1] = '\0';  // Ensure null termination

    return 0;

}

// ---------------- Parse + validate + map parameters ----------------


// Unified CFE config file parser driver to work with all CFE config file versions
//#####################
int parse_config_driver(const char* cfg_path,
                        double cfg_version,
                        CFE_CONFIG* config,
                        cfe_parameters_struct* params,
                        cfe_options_struct* options)
{
    // Parse config (v3 keyword format)
    int status_flag = parse_cfe_config(cfg_path, config);
    if (status_flag != 0) {
        fprintf(stderr, "ERROR: Failed to parse config file: %s\n", cfg_path);
        return status_flag;
    }

    // Step 2: Normalize and validate
    status_flag = normalize_config_units(config);
    if (status_flag != 0) {
        fprintf(stderr, "ERROR: Failed to normalize config units\n");
        return status_flag;
    }

    // Step 3: Unified parameter validation 
    status_flag = validate_required_parameters(config, options->verbosity);
    if (status_flag != 0) {
        fprintf(stderr, "ERROR: Config validation failed\n");
        return status_flag;
    }

    // Step 4: Map all parameters from config struct to CFE model structs 
    status_flag = map_config_to_parameters_and_options(config, params, options);
    if (status_flag != 0) {
        fprintf(stderr, "ERROR: Failed to map config to CFE structs\n");
        return status_flag;
    }

    return 0;  // Success
}



// ---------------- Initialize / Step / Finalize ----------------

//#################
int cfe_initialize(const cfe_parameters_struct* p,
                   const cfe_options_struct*    o,
                   cfe_state_struct*            s)
{
    
    if (p == NULL || o == NULL || s == NULL) {
        fprintf(stderr, "ERROR in cfe_initialize(): null pointer\n");
        return -1;
    }

    // write outputs to ./ by defaultt
    
    set_state_defaults(s);

    // Basic stores
    s->soil_storage_m = p->soil_init_storage_m;
    s->gw_storage_m   = p->gw_init_storage_m;

    // Clamp counts to fixed-size arrays
    int Nsub = p->nash_subsurface_N;
    if (Nsub < 0) Nsub = 0;
    if (Nsub > MAX_NUM_SUBSURFACE_NASH_CASCADE) Nsub = MAX_NUM_SUBSURFACE_NASH_CASCADE;

    int giuh_num_ords = p->giuh_num_ordinates;
    if (giuh_num_ords < 0) giuh_num_ords = 0;
    if (giuh_num_ords > MAX_NUM_GIUH_ORDINATES) giuh_num_ords = MAX_NUM_GIUH_ORDINATES;


    // ---- Subsurface Nash storage (fixed arrays; no malloc) ----
    // Nsub == MAX_NUM_SUBSURFACE_NASH_CASCADE (2) in practice, so the
    // ternary always takes the init branch; the zero branch is defensive.
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        s->nash_subsurface_storage_m[i] = (i < Nsub) ? p->nash_subsurface_init_storage_m[i] : 0.0;
    }

    // ---- GIUH queue (fixed array; no malloc) ----
    // Only meaningful iff using GIUH; otherwise zero it.
    if (o->surface_routing_scheme == SURF_ROUTE_GIUH && giuh_num_ords > 0) {
        for (int i = 0; i < giuh_num_ords; i++) {
            s->giuh_queue_m[i] = (i < giuh_num_ords) ? p->giuh_init_queue_m[i] : 0.0;
        }
    } else {
        for (int i = 0; i < MAX_NUM_GIUH_ORDINATES; i++) {
            s->giuh_queue_m[i] = 0.0;
        }
    }

    
    if (o->simulate_discrete_soil_moisture) {
    
        // copy parameters from CFE NWM soil parameters with needed unit conversions
        s->soil_parameters.theta_r = 0.0;  // Clapp-Hornberger assumes residual = 0
        s->soil_parameters.theta_sat = p->effective_porosity;
        s->soil_parameters.theta_fc = p->field_capacity_moisture_content;
        s->soil_parameters.theta_wp = p->wilting_point;
        s->soil_parameters.theta_aet_eq_pet = p->field_capacity_moisture_content; // reasonable default
        s->soil_parameters.K_sat_cm_per_h = p->ksat_m_per_s * 360000.0; // m/s to cm/h  
        s->soil_parameters.phi_sat_cm = p->sat_capillary_head_m * 100.0; // m to cm
        s->soil_parameters.b_exp = p->soil_b;
        s->soil_parameters.perc_limiter_0_to_1 = p->soil_to_gw_percolation_rate_limiter_0_1;
        s->soil_parameters.klf_per_h = p->soil_k_lateral_per_h; // linear reservoir rate constant

        // SoilControl - simulation parameters
        s->soil_control.ndisc = NDISC;
        s->soil_control.deepest_root_disc = o->deepest_root_zone_disc;
        if (s->soil_control.deepest_root_disc < 1) s->soil_control.deepest_root_disc = 1;
        if (s->soil_control.deepest_root_disc > NDISC) s->soil_control.deepest_root_disc = NDISC;
        s->soil_control.use_ch_lookup_table = o->use_soil_lookup_table ? 1 : 0;
        s->soil_control.is_sft_coupled = o->enable_freeze_thaw ? 1 : 0;
        s->soil_control.dt_hours = 1.0;  // CFE uses hourly timesteps

        if(o->verbosity > 1) {
           if(s->soil_control.use_ch_lookup_table) printf("DEBUG: Building lookup table...\n");
           else  printf("DEBUG: NOT building lookup table...\n");
        }
        
        // SoilGeometry - discretization thicknesses (from Noah-MP):
        s->soil_geometry.dz_m[0] = 0.1;
        s->soil_geometry.dz_m[1] = 0.3;
        s->soil_geometry.dz_m[2] = 0.6;
        s->soil_geometry.dz_m[3] = 1.0;
        
        s->soil_geometry.depth_m = p->soil_depth_m;
        
        // make sure that the soil_depth_m is equal to the sum of the dz_m
        double depth_check_m = 0.0;
        for(int i = 0; i < NDISC; i++) depth_check_m += s->soil_geometry.dz_m[i];
        if(!is_fabs_less_than_epsilon(depth_check_m - s->soil_geometry.depth_m, 1.0e-3) && o->verbosity > 0) {
             printf("Initialize WARNING: while setting up discrete soil moisture balance module\n");
             printf("input soil depth %.4f m differs from sum of disc thicknesses: %.4f m\n",
                     s->soil_geometry.depth_m, depth_check_m);
        }
           
        // depth from land surfact to the center of each disc
        s->soil_geometry.zc_m[0] = 0.05;
        s->soil_geometry.zc_m[1] = 0.25;
        s->soil_geometry.zc_m[2] = 0.70;
        s->soil_geometry.zc_m[3] = 1.50;
  
        double zwt_out = -999.9;
        // Determine which soil moisture initialization is provided (total storage, or theta in four discs)

        // If discrete values provided, check that ALL theta values are provided and > THETA_MIN 
        int has_discrete_theta = TRUE;
        for (int i = 0; i < NDISC; i++) {
            if (p->soil_discrete_init_theta[i] < THETA_MIN) {  // these should be zero iff not provided in config file
                has_discrete_theta = FALSE;
                break;
            }
        }

        if (has_discrete_theta) {
            s->soil_state_in.total_storage_m = 0.0;
            // Use provided discrete theta values
            for (int i = 0; i < NDISC; i++) {
                s->soil_discrete_storage_theta[i] = p->soil_discrete_init_theta[i];
                s->soil_state_in.total_storage_m += s->soil_discrete_storage_theta[i] * s->soil_geometry.dz_m[i];
            }
            if (o->verbosity > 0) {
                printf("Using provided discrete theta values: ");
                for (int i = 0; i < NDISC; i++) {
                    printf("%.3f ", s->soil_discrete_storage_theta[i]);
                }
                printf("\n");
            }
        } else {
            // Use hydrostatic initialization from total storage
            double soil_depth = p->soil_depth_m;
            double dz_m[NDISC], zc_m[NDISC];

            for (int i = 0; i < NDISC; i++) {
                dz_m[i] =  s->soil_geometry.dz_m[i];
                zc_m[i] =  s->soil_geometry.zc_m[i];
            }


            initialize_hydrostatic_from_storage(
                soil_depth,
                p->effective_porosity,          // theta_sat
                p->sat_capillary_head_m * 100.0, // phi_sat_cm (convert m to cm)
                p->soil_b,                      // b_exp
                zc_m,                            // disc centers
                p->soil_init_storage_m,         // target storage
                &zwt_out,                       // computed water table depth
                s->soil_discrete_storage_theta  // output theta array
            );

            if (o->verbosity > 0) {
                printf("Maximum soil moisture storage:                                %.6f m\n", p->soil_depth_m * p->effective_porosity);
                printf("Input initial total soil moisture storage:                    %.6f m\n", p->soil_init_storage_m);
                printf("Computed hydrostatic profile from total soil moisture storage %.6f m\n", p->soil_init_storage_m);
                printf("Computed fictitious depth to water table:                     %.6f m\n", zwt_out);
            }
        }
    
        // Debug output: always show final theta values and how they were determined
        if(o->verbosity > 1) {
            if (has_discrete_theta) {
                printf("DEBUG: Discrete soil moisture initialized from PROVIDED theta values:\n");
            } else {
                printf("DEBUG: Discrete soil moisture initialized from CALCULATED hydrostatic profile:\n");
                printf("       (from total storage %.6f m, computed fictitious water table depth %.6f m)\n", 
                       p->soil_init_storage_m, zwt_out);  // Note: zwt_out only available in else block
            }
            printf("       Final theta values: ");
            for (int i = 0; i < NDISC; i++) {
                printf("disc_%d=%.4f ", i+1, s->soil_discrete_storage_theta[i]);
            }
            printf("\n");
        }
            
        // Build lookup table if needed
        if (o->use_soil_lookup_table) {
            int rc = soil_build_ch_lut(
                p->lut_n_points,                  
                p->lut_theta_min,                 
                0.0,                              // theta_r (residual water content, 0.0 for Clapp-Hornberger)
                p->effective_porosity,            // theta_sat
                p->soil_b,                        // b_exp
                p->sat_capillary_head_m * 100.0,  // phi_sat_cm (convert m to cm)
                p->ksat_m_per_s * 360000.0,       // K_sat_cm_per_h (convert m/s to cm/h)
                &s->ch_lookup_tables
            );
            if (rc != 0) {
                fprintf(stderr, "ERROR: Failed to build Clapp-Hornberger lookup table in cfe_initialize\n");
                return -1;
            }
            if (o->verbosity > 0) {
                printf("Built Clapp-Hornberger lookup table with %d points\n", p->lut_n_points);
            }
        } else {
            // Initialize lookup table struct to zeros if not using it
            memset(&s->ch_lookup_tables, 0, sizeof(SoilLookupTables));
        }
        

        // Initialize SoilStateIn/Out structures
        if(!has_discrete_theta) s->soil_state_in.total_storage_m = 0.0;
        for (int i = 0; i < NDISC; i++) {
            s->soil_state_in.theta_in[i] = s->soil_discrete_storage_theta[i];  // copy initial theta
            s->soil_state_in.psi_in[i] = 0.0;           // will be computed each timestep
            s->soil_state_in.K_in[i] = 0.0;             // will be computed each timestep  
            s->soil_state_in.ch_lut_hint_in[i] = -1;    // initial hint
            s->soil_state_out.ch_lut_hint_out[i] = -1;
            if(!has_discrete_theta)
                s->soil_state_in.total_storage_m += s->soil_state_in.theta_in[i] * s->soil_geometry.dz_m[i];
        }

        // Initialize flux structure to zeros
        memset(&s->soil_fluxes, 0, sizeof(SoilFluxes));

        if (o->verbosity > 0) {
            printf("Built discrete soil simulation structures:\n");
            printf("  Discretizations: ");
            for(int i = 0; i < NDISC; i++) printf("%d:%.2f m ",i+1, s->soil_geometry.dz_m[i]);
            printf("\n");
            printf("  Root zone extends to and includes disc: %d\n", s->soil_control.deepest_root_disc);
            printf("  Lookup table used?: %s\n", s->soil_control.use_ch_lookup_table ? "Yes" : "No");
            if(s->soil_control.use_ch_lookup_table) printf("   consisting of %d points\n",p->lut_n_points);

            printf("Initialized discrete soil moisture with %d discretizations\n", NDISC);
            printf("Initial theta values: ");
            for (int i = 0; i < NDISC; i++) {
                printf("disc: %d %.3f ",i+1, s->soil_discrete_storage_theta[i]);
            }
            printf("\n");
        }
    }
 
    s->current_time_step = 0;
    return 0;
}

// Adapter: builds kernel structs from params/state/forcing and calls cfe()
//#######################################
int cfe_step(const cfe_parameters_struct* p,
             const cfe_options_struct* o,
             cfe_state_struct* s,
             const cfe_forcing_struct* forcing,
             double dt_seconds,
             cfe_outputs_struct* out,
             cfe_volbal_struct* volbal)
{

    if (p == NULL || o == NULL || s == NULL || forcing == NULL || out == NULL) {
        fprintf(stderr, "ERROR: in cfe_step(): null pointer\n");
        return -1;
    }

    
    // Keep deficits updated for use by kernel  
    double soil_max_storage_m = p->effective_porosity * p->soil_depth_m;
    s->soil_storage_deficit_m = soil_max_storage_m - s->soil_storage_m;
    if (s->soil_storage_deficit_m < 0.0) s->soil_storage_deficit_m = 0.0;

    s->gw_storage_deficit_m = p->gw_max_storage_m - s->gw_storage_m;
    if (s->gw_storage_deficit_m < 0.0) s->gw_storage_deficit_m = 0.0;

    /* TODO: These structs are rebuilt every step from constant parameters.
     * Consider persisting them in the state struct to avoid per-step copies.
     * Also: cfe() takes ~30 parameters — encapsulating into fewer structs
     * would improve the interface. Both are future refactoring targets. */

    // 1) Exchange soil parameters (NWM_SOIL_PARAMETERS_STRUCTURE)
    struct NWM_SOIL_PARAMETERS_STRUCTURE nwm = {0};
    nwm.smcmax   = p->effective_porosity;
    nwm.wltsmc   = p->wilting_point;
    nwm.satdk    = p->ksat_m_per_s;
    nwm.satpsi   = p->sat_capillary_head_m;
    nwm.bb       = p->soil_b;
    nwm.slop     = p->soil_to_gw_percolation_rate_limiter_0_1;  // bottom boundary factor 0..1
    nwm.D        = p->soil_depth_m;
    nwm.alpha_fc = p->field_capacity_Pcap_over_Patm;
    nwm.refkdt   = p->refkdt;
    nwm.soil_storage = s->soil_storage_m;
    nwm.wilting_point_m = nwm.wltsmc * nwm.D;

    // 2) Soil reservoir state/params
    struct CONCEPTUAL_RESERVOIR_STRUCTURE soil_res = {0};   // this is the soil linear reservoir (not routing)
    soil_res.is_exponential                = FALSE;  
    soil_res.storage_m                     = s->soil_storage_m;
    soil_res.storage_max_m                 = p->effective_porosity * p->soil_depth_m;
    soil_res.coeff_primary                 = nwm.satdk * nwm.slop * o->time_step_seconds;
    soil_res.exponent_primary              = 1.0;
    soil_res.storage_threshold_primary_m   = p->field_capacity_storage_m;  // calculated in parameter mapping func.

    soil_res.storage_threshold_secondary_m = p->field_capacity_storage_m;
    soil_res.coeff_secondary               = p->soil_k_lateral_per_h;     // this one is really important
    soil_res.exponent_secondary            = 1.0;  // linear
    soil_res.is_sft_coupled                = o->enable_freeze_thaw ? 1 : 0;
    soil_res.ice_fraction_schaake          = forcing->ice_fraction;
    soil_res.ice_fraction_xinanjiang       = forcing->ice_fraction;
    // NOTE: IN THIS ABOVE STRUCTURE IF I'M SIMULATING DISCRETE SOIL MOISTURE, THE ONLY THING NEEDED FROM IT
    //       IS THE soil_res.coeff_primary  THE REST IS NOT USED.

    
    // 3) GW reservoir state/params
    struct CONCEPTUAL_RESERVOIR_STRUCTURE gw_res = {0};
    gw_res.storage_max_m                   = p->gw_max_storage_m;
    gw_res.storage_m                       = s->gw_storage_m;
    gw_res.is_exponential                  = TRUE;  // CFE uses exponential groundwater discharge
    gw_res.coeff_primary                   = p->gw_discharge_coeff_m_per_timestep;
    gw_res.exponent_primary                = p->gw_discharge_exponent;
    gw_res.storage_threshold_primary_m     = 0.0;  // No threshold ffor primary outlet

    // Secondary outlet not used ffor groundwater
    gw_res.coeff_secondary                 = 0.0;
    gw_res.exponent_secondary              = 0.0;
    gw_res.storage_threshold_secondary_m   = 0.0;

    // 4) Partitioning parameters
    rainfall_partitioning_parameters_structure rp = {0};
    rp.surface_water_partitioning_scheme =
        (o->liquid_partitioning_scheme == PARTITION_XINANJIANG) ? PARTITION_XINANJIANG : PARTITION_SCHAAKE;
    // Schaake magic constant is normally derived; pick a conservative defaulty unless you already have one
    rp.Schaake_adjusted_magic_constant_by_soil_type = p->schaake_magic_constant;
    rp.a_Xinanjiang_inflection_point_parameter      = p->xj_tension_inflection_point;
    rp.b_Xinanjiang_shape_parameter                 = p->xj_tension_b;
    rp.x_Xinanjiang_shape_parameter                 = p->xj_free_b;
    rp.urban_decimal_fraction                       = p->catchment_impervious_fraction;
    rp.ice_content_threshold                        = p->soil_ice_imperv_threshold;

    // 5) Subsurface Nash parameters (use per hour in kernel)
   struct NASH_CASCADE_PARAMETERS_STRUCTURE nash_sub     = {0};
    

    nash_sub.N_nash                    = (p->nash_subsurface_N > 0) ? p->nash_subsurface_N : 0;
    nash_sub.K_nash                    = p->nash_subsurface_K_per_h;      // if this is per hour already, leave as-is
    nash_sub.nsubsteps                 = 1;
    nash_sub.nash_storage              = s->nash_subsurface_storage_m;

    nash_sub.is_riparian_gw = 0;

    // 6) ET structure from forcing
    evapotranspiration_structure et = {0};
    double forest_pet_m = forcing->forest_pet_m;

    if (forest_pet_m == 0.0 && forcing->bare_soil_aet_m == 0.0 &&
        forcing->et_potential_m > 0.0) {
        forest_pet_m = forcing->et_potential_m;
    }

    et.potential_et_m_per_timestep = forest_pet_m;
    et.potential_et_m_per_s = forest_pet_m / ((double)o->time_step_seconds);
    et.actual_bare_soil_evaporation_m_per_timestep = forcing->bare_soil_aet_m;


    // 7) Outputs and flux pointers expected by kernel
    double infiltration_excess_m = 0.0;
    double infiltration_depth_m  = 0.0;
    double flux_perc_m           = 0.0;
    double flux_lat_m            = 0.0;
    double flux_from_deep_gw_to_chan_m = 0.0;
    double flux_direct_runoff_to_channel_m = 0.0;
    double flux_nash_subsurface_lateral_runoff_m = 0.0;
    double Qout_m = 0.0;

    // GIUH arrays
    int num_giuh_ordinates = p->giuh_num_ordinates;
    double* giuh_ords = (double*)p->giuh_ordinates;
    double* giuh_queue = s->giuh_queue_m; 

    // Propagate external ice fraction to DSBM soil state
    s->soil_state_in.ice_fraction = forcing->ice_fraction;

    // add input rainfall and PET to the volbal structt this time step.
    volbal->volin += forcing->rainfall_depth_m;
    volbal->volin_PET += forcing->et_potential_m;

    // Snapshot volbal accumulators before kernel call so we can extract
    // per-timestep impervious/pervious runoff from the delta.
    const double vol_impervious_runoff_start = volbal->vol_impervious_runoff;
    const double vol_pervious_runoff_start = volbal->vol_pervious_runoff;

    // 8) Call cfe model
    cfe(
        &s->soil_storage_deficit_m,
        nwm,
        &soil_res,
        ((double)o->time_step_seconds) / 3600.0,  // timestep_h
        rp,
        forcing->rainfall_depth_m,
        &infiltration_excess_m,
        &infiltration_depth_m,
        &flux_perc_m,
        &flux_lat_m,
        &s->gw_storage_deficit_m,
        &gw_res,
        &flux_from_deep_gw_to_chan_m,
        &flux_direct_runoff_to_channel_m,
        num_giuh_ordinates,
        giuh_ords,
        giuh_queue,
        &flux_nash_subsurface_lateral_runoff_m,
        &nash_sub,
        &et,
        &Qout_m,
        volbal,
        (double)o->time_step_seconds,
        o->surface_routing_scheme,
        o->simulate_discrete_soil_moisture,  // flag
        &s->soil_control,                    // discrete soil control
        &s->soil_geometry,                   // discrete soil geometry 
        &s->soil_parameters,                 // Similar to NWM_Soil_Params struct, but different in key ways. 
        &s->soil_state_in,                   // discrete soil input state
        &s->soil_state_out,                  // discrete soil output state
        &s->soil_fluxes,                     // discrete soil fluxes
        &s->ch_lookup_tables,                // lookup tables
        s->soil_discrete_storage_theta,      // current theta values
        o->verbosity                         // 0 = quiet, >0 increasing degrees of stdout status reporting
    );

    // 9) Update state storages from the reservoir structs if kernel updated them
    s->soil_storage_m = soil_res.storage_m;
    s->gw_storage_m   = gw_res.storage_m;

    // 10) Write outputs

    out->impervious_runoff_m =
        volbal->vol_impervious_runoff - vol_impervious_runoff_start;
    out->pervious_runoff_m =
        volbal->vol_pervious_runoff - vol_pervious_runoff_start;
    out->surface_runoff_generated_m = infiltration_excess_m;
    out->surface_routed_to_outlet_m = flux_direct_runoff_to_channel_m;
    out->lateral_flow_generated_m = flux_lat_m;
    out->lateral_flow_m   = flux_nash_subsurface_lateral_runoff_m;
    out->baseflow_m       = flux_from_deep_gw_to_chan_m;
    out->total_outflow_m  = flux_direct_runoff_to_channel_m +
                            flux_nash_subsurface_lateral_runoff_m +
                            flux_from_deep_gw_to_chan_m;
    out->qout_m           = Qout_m;
    out->actual_et_m      = et.actual_et_m_per_timestep;
    out->bare_soil_evaporation_m = et.actual_bare_soil_evaporation_m_per_timestep;
    out->potential_et_m   = et.potential_et_m_per_timestep;
    out->giuh_outflow_m   = flux_direct_runoff_to_channel_m;
    out->soil_to_gw_percolation_flux_m = flux_perc_m;
    
    // Advance step count
    s->current_time_step += 1;
    (void)dt_seconds;
    return 0;
}

//########################
double cfe_get_last_qout_m(const cfe_outputs_struct* outputs) {
    return outputs ? outputs->qout_m : 0.0;
}


//##############
int cfe_resync_derived_params(cfe_parameters_struct* p,
                              const cfe_options_struct* o,
                              cfe_state_struct* s)
{
    if (p == NULL || o == NULL || s == NULL) return -1;

    // Recompute Schaake magic constant (depends on ksat)
    p->schaake_magic_constant = p->refkdt * p->ksat_m_per_s / 2.0e-06;

    // Recompute field capacity moisture content (depends on porosity, b, capillary head, fc fraction)
    double psi_atm_m = STANDARD_ATM_PRESS_Pa /
                       (GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2 * WATER_LIQUID_DENSITY_kg_per_m3);
    double arg = (p->field_capacity_Pcap_over_Patm * psi_atm_m / p->sat_capillary_head_m);
    p->field_capacity_moisture_content = p->effective_porosity * pow(arg, (-1.0 / p->soil_b));
    p->field_capacity_storage_m = p->field_capacity_moisture_content * p->soil_depth_m;

    // Resync DSBM soil_parameters copy
    if (o->simulate_discrete_soil_moisture) {
        s->soil_parameters.theta_sat           = p->effective_porosity;
        s->soil_parameters.theta_fc            = p->field_capacity_moisture_content;
        s->soil_parameters.theta_wp            = p->wilting_point;
        s->soil_parameters.theta_aet_eq_pet    = p->field_capacity_moisture_content;
        s->soil_parameters.K_sat_cm_per_h      = p->ksat_m_per_s * 360000.0;
        s->soil_parameters.phi_sat_cm          = p->sat_capillary_head_m * 100.0;
        s->soil_parameters.b_exp               = p->soil_b;
        s->soil_parameters.perc_limiter_0_to_1 = p->soil_to_gw_percolation_rate_limiter_0_1;
        s->soil_parameters.klf_per_h           = p->soil_k_lateral_per_h;

        // Rebuild Clapp-Hornberger lookup table if active
        if (o->use_soil_lookup_table) {
            soil_free_ch_lut(&s->ch_lookup_tables);
            soil_build_ch_lut(
                p->lut_n_points,
                p->lut_theta_min,
                0.0,
                p->effective_porosity,
                p->soil_b,
                p->sat_capillary_head_m * 100.0,
                p->ksat_m_per_s * 360000.0,
                &s->ch_lookup_tables
            );
        }
    }

    return 0;
}

//##############
int cfe_finalize(cfe_state_struct* s)
{
    if (s && s->ch_lookup_tables.lnpsi != NULL) {
        soil_free_ch_lut(&s->ch_lookup_tables);
    }
    (void)s; // nothing to free; arrays are static
    return 0;
}


//########################
void print_cfe_input_debug(const cfe_options_struct*    o,
                           const cfe_parameters_struct* p,
                           const cfe_state_struct*      s,
                           const char*                  cfg_path,
                           const char*                  forcing_path,
                           const char*                  qout_path,
                           const char*                  volbal_path)
{
    printf("===== CFE INPUT SUMMARY =====\n");

    // Run control & outputs
    printf("Config file:                    %s\n", cfg_path     ? cfg_path     : "(none)");
    printf("Forcing file:                   %s\n", forcing_path ? forcing_path : "(none)");
    printf("Qout file:                      %s\n", qout_path    ? qout_path    : "(none)");
    printf("VolBal file:                    %s\n", volbal_path  ? volbal_path  : "(none)");
    if (o) {
        printf("Time step:                      %d [s] (%.2f [h])\n", o->time_step_seconds, o->time_step_seconds/3600.0);
        printf("Num timesteps:                  %d\n", o->num_timesteps);
        printf("Partitioning scheme:            %s\n", partition_scheme_name(o->liquid_partitioning_scheme));
        printf("Surface routing scheme:         %s\n", surf_scheme_name(o->surface_routing_scheme));
    };
    printf("\n");

    // Construct NWM structure from parameters for display
    struct NWM_SOIL_PARAMETERS_STRUCTURE nwm = {0};
    nwm.smcmax   = p->effective_porosity;
    nwm.wltsmc   = p->wilting_point;
    nwm.satdk    = p->ksat_m_per_s;
    nwm.satpsi   = p->sat_capillary_head_m;
    nwm.bb       = p->soil_b;
    nwm.slop     = p->soil_to_gw_percolation_rate_limiter_0_1;
    nwm.D        = p->soil_depth_m;
    nwm.alpha_fc = p->field_capacity_Pcap_over_Patm;
    nwm.refkdt   = p->refkdt;
    nwm.soil_storage = p->soil_init_storage_m;
    nwm.wilting_point_m = nwm.wltsmc * nwm.D;

    printf("\n-- Noah-MP Soil Parameters --\n");
    printf("Depth D:                        %.3f [m]\n",   nwm.D);
    printf("Clapp-Hornberger b:             %.3f [-]\n",   nwm.bb);
    double Ksat_cm_per_h = nwm.satdk * 360000.0;
    printf("Ksat:                           %.2f [cm h-1], %.6e [m s^-1]\n", Ksat_cm_per_h, nwm.satdk);
    printf("Saturation capillary head:      %.3f [m]\n",   nwm.satpsi);
    printf("Porosity (smcmax):              %.3f [V/V]\n", nwm.smcmax);
    printf("Wilting point (wltsmc):         %.3f [V/V]\n", nwm.wltsmc);
    printf("Field capacity alpha_fc:        %.3f [-]\n",   nwm.alpha_fc);
    printf("Calculated field capacity:      %.3f [V/V] (%.3f [m])\n", 
             p->field_capacity_moisture_content, p->field_capacity_storage_m);
    printf("refkdt:                         %.3f [-]\n",   nwm.refkdt);
    printf("percolation limiter (0-1):      %.3f [-]\n",   nwm.slop);
    printf("soil_storage (init est):        %.6f [m]\n",   nwm.soil_storage);
    // Initial storages (from state)
    if (s) {
        printf("\n-- Initial Storages --\n");
        printf("Soil storage:                   %.6f [m]\n", s->soil_storage_m);
        printf("Groundwater storage:            %.6f [m]\n", s->gw_storage_m);
        printf("\n");
    }

    if (!p) { printf("====================================\n\n"); return; }

    // Rainfall partitioning
    printf("-- Surface Partitioning --\n");
    // If you added the limiter parameter (soil->GW)
    #ifdef HAVE_TO_GW_PERC_LIMITER
    printf("Soil->GW perc limiter (0..1):   %.3f [-]\n", p->to_gw_perc_limiter_0_1);
    #endif
    printf("\n");

    // GIUH
    if (o && o->surface_routing_scheme == SURF_ROUTE_GIUH) {
        int G = p->giuh_num_ordinates;
        if (G < 0) G = 0;
        if (G > MAX_NUM_GIUH_ORDINATES) G = MAX_NUM_GIUH_ORDINATES;

        printf("-- GIUH --\n");
        printf("Num GIUH ordinates:             %d\n", G);
        if (G > 0) {
            printf("GIUH ordinates:                ");
            for (int i = 0; i < G; i++) printf(" %.6g", p->giuh_ordinates[i]);
            printf("\nInitial GIUH queue (m):        ");
            for (int i = 0; i < G; i++) printf(" %.6g", p->giuh_init_queue_m[i]);
            printf("\n");
        }
        printf("\n");
    }


    // Subsurface Nash (always show, since you hard-coded N=2)
    {
        int Nsub = p->nash_subsurface_N;
        if (Nsub < 0) Nsub = 0;
        if (Nsub > MAX_NUM_SUBSURFACE_NASH_CASCADE) Nsub = MAX_NUM_SUBSURFACE_NASH_CASCADE;

        printf("-- Subsurface Nash --\n");
        printf("N (reservoirs):                 %d\n", Nsub);
        printf("K (per hour):                   %.6f [h^-1]\n", p->nash_subsurface_K_per_h);

        if (Nsub > 0) {
            printf("Init storages (m):             ");
            for (int i = 0; i < Nsub; i++) printf(" %.6f", p->nash_subsurface_init_storage_m[i]);
            printf("\n");
        }
        if (s && Nsub > 0) {
            printf("Current storages (m):          ");
            for (int i = 0; i < Nsub; i++) printf(" %.6f", s->nash_subsurface_storage_m[i]);
            printf("\n");
        }
        printf("\n");
    }

    printf("\n-- Groundwater Parameters --\n");
    printf("GW max storage:                 %.6f [m]\n", p->gw_max_storage_m);
    printf("GW discharge coeff:             %.6e [m h^-1]\n", p->gw_discharge_coeff_m_per_timestep);
    printf("GW discharge exponent:          %.3f [-]\n", p->gw_discharge_exponent);
    printf("\n");

    // Soil & GW initial depths (from params, since they come from config)
    printf("-- Soil/GW Initial Conditions --\n");
    printf("Soil init storage:              %.6f [m]\n", p->soil_init_storage_m);
    printf("GW init storage:                %.6f [m]\n",   p->gw_init_storage_m);

    printf("-- Initial Storages --\n");
    printf("Soil (config / state):          %.6f / %.6f [m]\n",
           p ? p->soil_init_storage_m : -1.0,
           s ? s->soil_storage_m      : -1.0);
    printf("GW   (config / state):          %.6f / %.6f [m]\n",
           p ? p->gw_init_storage_m : -1.0,
           s ? s->gw_storage_m      : -1.0);
    printf("\n");

    // If your params include NWM soil properties, add them here using your actual field names.
    // Example (uncomment and correct names if present):
    //printf("Depth D:                        %.3f [m]\n",   p->D);
    //printf("Clapp-Hornberger b:             %.3f [-]\n",   p->bb);
    //printf("Ksat:                           %.6e [m s^-1]\n", p->satdk);
    //printf("Sat capillary head:             %.3f [m]\n",   p->satpsi);
    //printf("Porosity (smcmax):              %.3f [V/V]\n", p->smcmax);
    //printf("Wilting point (wltsmc):         %.3f [V/V]\n", p->wltsmc);
    //printf("Field capacity alpha_fc:        %.3f [-]\n",   p->alpha_fc);

    printf("====================================\n\n");
}

void print_exchange_values(int timestep,
                          const cfe_parameters_struct* p,
                          const cfe_options_struct* o,
                          const cfe_state_struct* s,
                          const cfe_forcing_struct* forcing,
                          double dt_seconds,
                          const cfe_outputs_struct* outputs)
{
    printf("\n========== EXCHANGE VALUES DEBUG - TIMESTEP %d ==========\n", timestep);
    
    // FORCING INPUTS
    printf("-- Forcing Inputs --\n");
    printf("Rainfall depth:                 %.6f [m]\n", forcing ? forcing->rainfall_depth_m : -999.0);
    printf("Potential ET:                   %.6f [m]\n", forcing ? forcing->et_potential_m : -999.0);
    printf("Time step:                      %.1f [s]\n", dt_seconds);
    
    // CURRENT STATE 
    printf("\n-- Current State --\n");
    if (s) {
        printf("Soil storage:                   %.6f [m]\n", s->soil_storage_m);
        printf("Soil deficit:                   %.6f [m]\n", s->soil_storage_deficit_m);
        printf("GW storage:                     %.6f [m]\n", s->gw_storage_m);
        printf("GW deficit:                     %.6f [m]\n", s->gw_storage_deficit_m);
        printf("Current timestep:               %d\n", s->current_time_step);
        
        // Nash surface storage
        if (o && o->surface_routing_scheme == SURF_ROUTE_NASH_CASCADE) {
            printf("Nash surface storage:          ");
            for (int i = 0; i < MAX_NUM_SURFACE_NASH_CASCADE; i++) {
                printf(" %.6f", s->nash_surface_storage_m[i]);
            }
            printf(" [m]\n");
        }
        
        // Nash subsurface storage  
        printf("Nash subsurface storage:       ");
        for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
            printf(" %.6f", s->nash_subsurface_storage_m[i]);
        }
        printf(" [m]\n");
        
        // GIUH queue
        if (o && o->surface_routing_scheme == SURF_ROUTE_GIUH) {
            printf("GIUH queue:                    ");
            for (int i = 0; i < 5; i++) {  // Show first 5 values
                printf(" %.6f", s->giuh_queue_m[i]);
            }
            printf(" [m]\n");
        }
    }
    
    // KEY PARAMETERS (subset)
    printf("\n-- Key Parameters --\n");
    if (p) {
        printf("Soil depth:                     %.3f [m]\n", p->soil_depth_m);
        printf("Effective porosity:             %.3f [V/V]\n", p->effective_porosity);
        printf("Field capacity:                 %.3f [V/V] (%.6f [m])\n", 
               p->field_capacity_moisture_content, p->field_capacity_storage_m);
        printf("Wilting point:                  %.3f [V/V]\n", p->wilting_point);
        printf("Ksat:                           %.6e [m/s]\n", p->ksat_m_per_s);
        printf("GW max storage:                 %.6f [m]\n", p->gw_max_storage_m);
        printf("GW discharge coeff:             %.6e [m/s]\n", p->gw_discharge_coeff_m_per_timestep);
        printf("GW discharge exponent:          %.3f [-]\n", p->gw_discharge_exponent);
    }
    
    // OPTIONS
    printf("\n-- Options --\n");
    if (o) {
        printf("Partitioning scheme:            %s\n", partition_scheme_name(o->liquid_partitioning_scheme));
        printf("Surface routing:                %s\n", surf_scheme_name(o->surface_routing_scheme));
        printf("Time step seconds:              %d [s]\n", o->time_step_seconds);
    }
    
    // OUTPUTS (from previous step)
    printf("\n-- Previous Outputs --\n");
    if (outputs) {
        printf("Surface runoff:                 %.6f [m]\n", outputs->surface_runoff_generated_m);
        printf("Lateral subsurface flow:        %.6f [m]\n", outputs->lateral_flow_m);
        printf("Baseflow:                       %.6f [m]\n", outputs->baseflow_m);
        printf("Total Qout:                     %.6f [m]\n", outputs->qout_m);
    }
    
    printf("============================================================\n\n");
}

void check_dsbm_local_volume_balance(
    TimestepSoilVolbal *soil_volbal,
    SoilFluxes *soil_fluxes,
    SoilStateIn *soil_state_in,
    SoilStateOut *soil_state_out,
    SoilGeometry *soil_geometry,
    double infiltration_depth_m,
    double actual_et_from_soil_m,
    double balance_tolerance
) {
    // DSBM input
    double dsbm_input_total = infiltration_depth_m;
    
    // DSBM outputs (all in meters per timestep)
    double dsbm_perc_out = soil_volbal->perc_m;
    double dsbm_lateral_out = soil_volbal->lateral_m;
    double dsbm_excess_out = soil_fluxes->rain_excess_m;
    double dsbm_et_out = actual_et_from_soil_m;
    
    // Calculate storage change
    double dsbm_storage_initial = 0.0;
    double dsbm_storage_final = 0.0;
    
    for(int i = 0; i < NDISC; i++) {
        dsbm_storage_initial += soil_state_in->theta_in[i] * soil_geometry->dz_m[i];
        dsbm_storage_final += soil_state_out->theta_out[i] * soil_geometry->dz_m[i];
    }
    double dsbm_storage_change = dsbm_storage_final - dsbm_storage_initial;
    
    // Volume balance: INPUT = OUTPUT + STORAGE_CHANGE
    // NOTE: ET is already accounted ffor in storage_change, so don't double-count 
    
    double dsbm_output_total = dsbm_perc_out + dsbm_lateral_out + dsbm_excess_out + dsbm_et_out;
 //   printf("DEBUG: perc=%e, lat=%e, excess=%e, sum=%e, et=%e\n", 
 //      dsbm_perc_out, dsbm_lateral_out, dsbm_excess_out, dsbm_output_total, dsbm_et_out);

    double dsbm_balance_error = dsbm_input_total - dsbm_output_total - dsbm_storage_change;
    
    // Check for volume balance closure
    if (fabs(dsbm_balance_error) > balance_tolerance) {
        printf("*** DSBM LOCAL VOLUME BALANCE ERROR ***\n");
        printf("  Balance error: %e m (tolerance: %e)\n", dsbm_balance_error, balance_tolerance);
        printf("  Input total:   %e m\n", dsbm_input_total);
        printf("  Output total:  %e m (includes ET)\n", dsbm_output_total);
        printf("    - Percolation:     %e m\n", dsbm_perc_out);
        printf("    - Lateral flow:    %e m\n", dsbm_lateral_out);
        printf("    - Excess runoff:   %e m\n", dsbm_excess_out);
        printf("    - ET from soil:    %e m\n", dsbm_et_out);
        printf("  Storage change:%e m\n", dsbm_storage_change);
        printf("    - Initial:         %e m\n", dsbm_storage_initial);
        printf("    - Final:           %e m\n", dsbm_storage_final);
        printf("  Volume balance: %.6e = %.6e + %.6e\n", 
               dsbm_input_total, dsbm_output_total, dsbm_storage_change);
        printf("******************************************\n");
    }
    
    // Check for negative fluxes
    if (dsbm_perc_out < 0.0 || dsbm_lateral_out < 0.0 || dsbm_excess_out < 0.0 || dsbm_et_out < 0.0) {
        printf("*** WARNING: DSBM produced negative flux(es) ***\n");
        printf("  Percolation: %e, Lateral: %e, Excess: %e, ET: %e\n", 
               dsbm_perc_out, dsbm_lateral_out, dsbm_excess_out, dsbm_et_out);
    }
    
    // Check physical reasonableness
    if (dsbm_output_total > dsbm_input_total && dsbm_storage_change >= 0.0) {
        printf("*** WARNING: DSBM output exceeds input with no storage decrease ***\n");
        printf("  Output/Input ratio: %.6f\n", 
               dsbm_input_total > 0.0 ? dsbm_output_total / dsbm_input_total : 0.0);
    }
}

