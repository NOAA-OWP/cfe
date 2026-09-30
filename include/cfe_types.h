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
 

#ifndef CFE_TYPES_H
#define CFE_TYPES_H

#include <stddef.h>

#include "nash_cascade.h"
#include "giuh.h"
#include "cfe_soil_discrete.h"
#include "cfe_config.h"


#include <math.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define CLAMP(val, lo, hi) fmin(fmax((val), (lo)), (hi))

// used in parsing command line args
typedef struct {
    char* cfg_path;
    char* forcing_path;
    char* qout_path;
    char* volbal_path;
    char* fluxes_path;
    char* stores_path;
    char* thetas_path;
    int command_line_verbosity;
    int run_without_forcing;
} cfe_cmdline_args_struct;

/* struct to store the temporal properties of the forcing file */
typedef struct {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    int delta_t_seconds;
    int num_valid_lines;
    char time_format[64];  // Store the original time format from forcing file
} aorc_forcing_time_struct;

//###########
typedef struct {
    int    ndisc;                 // must equal NDISC
    int    deepest_root_disc;     // 1..ndisc
    int    use_ch_lookup_table;   // 1 => use LUT; 0 => analytic CH
    int    is_sft_coupled;        // TRUE iff coupled to soil freeze-thaw 
    double dt_hours;              // usually 1.0
} SoilControl;

//###########
typedef struct {
    double dz_m[NDISC];
    double zc_m[NDISC];
    double depth_m;
} SoilGeometry;

//###########
typedef struct {
    double theta_r;               // residual saturation (m3/m3)
    double theta_sat;             // saturation (m3/m3)
    double theta_fc;              // field capacity (m3/m3)
    double theta_wp;              // wilting point (m3/m3)
    double theta_aet_eq_pet;      // AET = PET at or above this (m3/m3)

    double K_sat_cm_per_h;        // cm/h
    double phi_sat_cm;            // cm
    double b_exp;                 // Clapp-Hornberger exponent

    double perc_limiter_0_to_1;   // 0..1 bottom drainage limiter
    double klf_per_h;           // lateral removal rate constant (m/h)
} SoilParameters;

// CH lookup tables over Theta=(theta-theta_r)/(theta_sat-theta_r)
//###########
typedef struct {
    int    n;
    double lnTheta_min;
    double dlnTheta;
    double inv_dlnTheta;
    double *lnpsi;                // ln(psi[m]) length n
    double *lnK;                  // ln(K[m/h]) length n
} SoilLookupTables;

//###########
typedef struct {
     double theta_in[NDISC];
     double psi_in[NDISC];   // m
     double K_in[NDISC];     // m/h
     int    ch_lut_hint_in[NDISC];
     double ice_fraction;
     double total_storage_m;
     double storage_deficit_m;
} SoilStateIn;

//###########
typedef struct {
    double theta_out[NDISC];
    int    ch_lut_hint_out[NDISC];
    double total_storage_m;
    double storage_deficit_m;
} SoilStateOut;

//###########
typedef struct {
    double rain_mm_per_h;         // mm/h
    double pet_mm_per_h;          // mm/h
} SoilForcing;

//###########
typedef struct {
    // Step-integrated exchanges (m)
    double AET_by_disc_m[NDISC];
    double lateral_by_disc_m[NDISC];
    double percolation_to_gw_m;
    double rain_into_soil_m;
    double rain_excess_m;

    // Internal vertical exchanges: [0..NDISC-2] interfaces i to i+1; [NDISC-1] bottom perc
    double interface_vol_m[NDISC];
    double interface_rate_m_per_h[NDISC];

    int    n_substeps_used;
} SoilFluxes;

typedef struct {
    double in_rain_m;         // infiltrated
    double excess_m;          // rejected
    double perc_m;
    double AET_m;
    double lateral_m;
    double delta_storage_m;   // Sum(theta_out-theta_in)*dz
    double residual_m;        // in - (outs) - Delta_S
} TimestepSoilVolbal;

/* Partitioning and routing options */
//###########
typedef enum {
    PARTITION_SCHAAKE = 0,
    PARTITION_XINANJIANG = 1
} cfe_partition_scheme_t;

//############
typedef enum {
    SURF_ROUTE_GIUH = 0,
    SURF_ROUTE_NASH_CASCADE = 1
} cfe_surface_route_t;

/* Options: switches and run controls (from config) */
//############
typedef struct {
    double cfe_version;           // from cfe_config_version key; defaults to 3.0
    cfe_partition_scheme_t liquid_partitioning_scheme;
    cfe_surface_route_t surface_routing_scheme;
    int time_step_seconds;
    int num_timesteps;
    int verbosity;
    int enable_ET_Priestley_Taylor;
    int enable_freeze_thaw;
    int simulate_discrete_soil_moisture;
    int simulate_soil_evaporation;
    int deepest_root_zone_disc;       // 1..NDISC; used by DSBM ET extraction
    int use_soil_lookup_table;
    char input_forcing_filename[PATH_FILENAME_STRING_LENGTH];

    // Output configuration
    char output_discharge_filename[PATH_FILENAME_STRING_LENGTH];
    char output_total_discharge_m3_per_sec_filename[PATH_FILENAME_STRING_LENGTH];
    char output_status_warnings_filename[PATH_FILENAME_STRING_LENGTH];
    char output_internal_fluxes_filename[PATH_FILENAME_STRING_LENGTH];
    char output_internal_storages_filename[PATH_FILENAME_STRING_LENGTH];
    char output_volume_balance_filename[PATH_FILENAME_STRING_LENGTH];
    char output_soil_moisture_theta_filename[PATH_FILENAME_STRING_LENGTH];
    char output_time_standard_format[OPTION_STRING_LENGTH];
    char output_file_delimiter[OPTION_STRING_LENGTH];
    char output_value_format[OPTION_STRING_LENGTH];
    char output_path_name[PATH_FILENAME_STRING_LENGTH];
    char output_new_config_filename[PATH_FILENAME_STRING_LENGTH];
    double epoch_start_seconds;   // Unix epoch of simulation start; 0 = not configured
} cfe_options_struct;

/* Parameters: physical/empirical constants and arrays (from config) and initial states */
//############
typedef struct {
    /* Catchment */
    double catchment_area_km2;
    
    /* Soil */
    double soil_depth_m;
    double soil_b;
    double ksat_m_per_s;
    double sat_capillary_head_m;
    double effective_porosity;
    double wilting_point;
    double field_capacity_Pcap_over_Patm;
    double field_capacity_moisture_content; 
    double field_capacity_storage_m;        
    double soil_init_storage_m;
    double soil_k_lateral_per_h;
    double soil_ice_imperv_threshold;
    double refkdt;  // strictly not a parameter, it is a constant = 3.0
    double schaake_magic_constant;
    double soil_to_gw_percolation_rate_limiter_0_1;
    
    /* for discrete soil simulation */
    double soil_discrete_init_theta[NDISC];
    int    lut_n_points;
    double lut_theta_min;
    
    /* for using Priestley-Taylor method to calculate PET  */
    double alpha_pt;   // orginally 1.26, but often smaller for deserts and can be larger

    /* Xinanjiang, read even if not used */
    double xj_tension_inflection_point;
    double xj_tension_b;
    double xj_free_b;

    /* Groundwater */
    double gw_max_storage_m;
    double gw_init_storage_m;
    double gw_discharge_coeff_m_per_timestep;
    double gw_discharge_exponent;

    /* Surface routing */
    int giuh_num_ordinates;
    double giuh_ordinates[MAX_NUM_GIUH_ORDINATES];
    double giuh_init_queue_m[MAX_NUM_GIUH_ORDINATES];


    /* Subsurface Nash */
    int nash_subsurface_N;
    double nash_subsurface_K_per_h;
    double nash_subsurface_init_storage_m[2];

    /* Land-cover fractions */
    double catchment_vegetated_fraction;
    double catchment_impervious_fraction;
    double catchment_bare_soil_fraction;

    /* Bare-soil evaporation */
    double bare_soil_rsurf_exp;

} cfe_parameters_struct;

/* Priestley-Taylor soil temperature state */
typedef struct {
    double skin_temperature_k;
    double upper_soil_temperature_k;
    double estimated_annual_air_temperature_k;
    double air_temperature_time_integral_k_s;
    double accumulated_time_s;
    int initialized;
} cfe_pet_temperature_state_struct;

/* State: storages and internal queues that evolve over time.
 *
 * SERIALIZATION NOTE — the fields marked [SERIALIZED] below are the
 * minimum set required to checkpoint and restore a running model.
 * They are packed into a binary buffer by cfe_serialize_create()
 * and restored by cfe_serialize_deserialize().
 *
 * If you add, remove, or resize a serialized field, follow the
 * developer checklist in cfe_serialize.h.  The byte layout diagram
 * is in cfe_serialize.c.
 *
 * Remaining fields (deficit values, lookup tables, DSBM sub-structs,
 * current_time_step) are either recomputed each timestep or set once
 * during initialization and do not change — they are NOT serialized.
 */
//############
typedef struct {
    double soil_storage_m;              /* [SERIALIZED] conceptual soil reservoir */
    double soil_storage_deficit_m;
    double gw_storage_m;                /* [SERIALIZED] groundwater reservoir */
    double gw_storage_deficit_m;

    double nash_surface_storage_m[MAX_NUM_SURFACE_NASH_CASCADE];
    double nash_subsurface_storage_m[MAX_NUM_SUBSURFACE_NASH_CASCADE]; /* [SERIALIZED] */
    double giuh_queue_m[MAX_NUM_GIUH_ORDINATES];                      /* [SERIALIZED] */
    double soil_discrete_storage_theta[NDISC];                        /* [SERIALIZED] */

    // Lookup tables for discrete soil moisture (calculated once during init)
    SoilLookupTables ch_lookup_tables;

    // DSBM sub-structs — initialized once, updated each timestep
    SoilControl      soil_control;
    SoilGeometry     soil_geometry;
    SoilParameters   soil_parameters;
    SoilStateIn      soil_state_in;     /* theta_in synced from soil_discrete_storage_theta on restore */
    SoilStateOut     soil_state_out;
    SoilFluxes       soil_fluxes;

    int current_time_step;

    // Priestley-Taylor soil temperature state
    cfe_pet_temperature_state_struct pet_temperature_state;  /* [SERIALIZED] */
} cfe_state_struct;


/*
 * DLWRF_surface [W m-2] is a real observed/modeled quantity and is
 * physically always positive.  These constants let the energy-balance
 * code distinguish "genuinely supplied and physically plausible" from
 * "never set / garbage / sensor dropout" and fall back to a synthetic
 * clear-sky estimate in the latter case.
 */
#define CFE_DLWRF_SURFACE_UNINITIALIZED_SENTINEL (-9999.0)
#define CFE_DLWRF_SURFACE_MIN_VALID_W_PER_M2      (50.0)
#define CFE_DLWRF_SURFACE_MAX_VALID_W_PER_M2      (600.0)

/*
 * Aerodynamic resistance [s/m] to atmosphere<->surface transfer.  Used
 * for sensible-heat transfer to the skin (cfe_soil_skin_temperature.c)
 * and vapor transfer to the soil surface (bare-soil evaporation).
 */
#define CFE_AERODYNAMIC_RESISTANCE_S_PER_M (100.0)

/*
 * Noah-MP RSURF_EXP: Sakaguchi-Zeng dry-layer shape exponent for
 * bare-soil surface resistance.  Common default is 5.0.
 */
#define CFE_BARE_SOIL_RSURF_EXP (5.0)

/* Forcing: inputs per step */
//############
typedef struct {
    double rainfall_depth_m;      // PUT THIS FIRST - CFE expects it here!
    double et_potential_m;        // PUT THIS SECOND - CFE expects it here!
    double APCP_surface;          // mm (converted to m for rain_m)
    double DLWRF_surface;         // W/m^2
    double DSWRF_surface;         // W/m^2
    double PRES_surface;          // Pa
    double SPFH_2maboveground;    // kg/kg
    double TMP_2maboveground;     // K
    double UGRD_10maboveground;   // m/s
    double VGRD_10maboveground;   // m/s
    double precip_rate;           // m/s
    double ice_fraction;          // 0-1, from external SFT module when coupled
    double forest_pet_m;          // forest-fraction PET for root-zone AET pathway
    double bare_soil_aet_m;       // area-weighted bare-soil evaporation [m/timestep]
    int day_of_year;              // 1 through 366
    int day_of_year_set_externally; // TRUE if set via BMI set_value; suppresses internal calc
} cfe_forcing_struct;

typedef enum {
    FORCING_FMT_UNKNOWN = 0,
    FORCING_FMT_DEV,        // Fred's standalone dev format (APCP_surface, DLWRF_surface, ...)
    FORCING_FMT_NGEN        // ngen/AORC-standard format (RAINRATE, LWDOWN, SWDOWN, ...)
} forcing_format_t;

typedef struct {
    int time_idx;
    int apcp_idx;
    int precip_rate_idx;
    int dlwrf_idx;
    int dswrf_idx;
    int pres_idx;
    int spfh_idx;
    int tmp_idx;
    int ugrd_idx;
    int vgrd_idx;
    int precip_rate_is_kg_m2_s1;
    forcing_format_t detected_format;
} aorc_cols_t;

//typedef struct { <- old pre 3.0
//    double rainfall_depth_m;
//    double et_potential_m;
//} cfe_forcing_struct;



// In cfe_types.h, add this BEFORE the CFE_Model_Context definition:

typedef struct {
    double volstart            ;
    double volstart_soil       ;
    double volstart_gw         ;
    double volstart_surface    ;
    double volstart_subsurface ;
    double vol_direct_runoff   ;  // this is water that won't fit into the soil suring a timestep because it is full
                                  // particularly in the discretized soil situation, where the upper disc fills
    double vol_impervious_runoff;
    double vol_pervious_runoff  ;
    double vol_runoff          ;
    double vol_infilt          ;
    double vol_out_surface     ;
    double vol_end_surface     ;
    double vol_to_gw           ;
    double vol_in_gw_start     ;
    double vol_in_gw_end       ;
    double vol_from_gw         ;
    double vol_in_subsurf_nash ;   
    double vol_out_subsurf_nash;   
    double vol_soil_start      ;
    double vol_to_soil         ;
    double vol_soil_to_lat_flow;
    double vol_soil_to_gw      ;
    double vol_soil_end        ;
    double vol_et_from_soil    ;
    double vol_forest_aet      ;
    double vol_bare_soil_evaporation;
    double vol_et_from_rain    ;
    double vol_et_to_atm       ;
    double volin               ;
    double volin_PET           ;
    double volout              ;
    double volend_soil         ;
    double volend_gw           ;
    double volend_surface      ;
    double volend_subsurface   ;
    double volend              ;


    /* ngen BMI mass balance protocol fields */
    double cumulative_vol;      /* volstart + volin (ngen::mass_in) */
    double volume_in_domain;    /* total storage at end of timestep (ngen::mass_stored) */
    double leakage;             /* future deep-GW / boundary losses (ngen::mass_leaked) */
} cfe_volbal_struct;

/* output configuration helper structure */
typedef struct {
    char time_format[64];      // "timestep", "datetime", "juliandate"
    char delimiter[8];         // ",", " ", "\t", "|"
    int current_timestep;
    aorc_forcing_time_struct* forcing_time;
} cfe_output_config_struct;

/* Outputs: results per step */
//############
typedef struct {
    double impervious_runoff_m;
    double pervious_runoff_m;
    double surface_runoff_generated_m;
    double surface_routed_to_outlet_m;
    double lateral_flow_generated_m;
    double lateral_flow_m;
    double baseflow_m;
    double total_outflow_m;
    double qout_m;
    double actual_et_m;
    double bare_soil_evaporation_m;
    double potential_et_m;
    double giuh_outflow_m;
    double soil_to_gw_percolation_flux_m;
} cfe_outputs_struct;

/* CFE model context: the whole model instance bundled together */
//############
typedef struct {
    CFE_CONFIG              config;        /* everything read in from a config file in one package */
    cfe_options_struct      options;       /* switches and run control */
    cfe_parameters_struct   parameters;    /* physical parameters */
    cfe_state_struct        state;         /* model storages and internal state */
    cfe_forcing_struct      forcing;       /* current step forcing */
    cfe_outputs_struct      last_outputs;  /* results from last step */
    cfe_volbal_struct       volbal;        /* the volume balance struct */
    double timestep_storage_start_m;       /* four per-timestep volume balance measures*/
    double timestep_input_m;
    double timestep_output_m;
    double timestep_storage_end_m;
    double vol_balance_residual_m;  /* cached: volstart + volin - volout - volend */
    /* Dirty flag: set by BMI set_value on calibration params, cleared after resync */
    int params_dirty;
    /* Serialization protocol buffer (ngen::serialization_*) */
    char  *serialized_state;
    size_t serialized_size;
} CFE_Model_Context;

#endif
