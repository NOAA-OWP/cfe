#ifndef CFE_BMI_CFE_H
#define CFE_BMI_CFE_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "./cfe.h"
#include "bmi.h"

//--------------------------------------------------
// Experiment to simplify BMI implementation (SDP)
//--------------------------------------------------
typedef struct Variable{
    unsigned int index;
    char name[80];
    char type[80];
    unsigned int size;
    // char units[80];
    // char role[80];
    // bool is_pointer;
} Variable;

/** Read number of lines in file and max line length, returning -1 if it does not exist or could not be read. */
int read_file_line_counts_cfe(const char* file_name, int* line_count, int* max_line_length);

/*  DEFINE A STRUCTURE TO HOLD THE ENTIRE MODEL STATE  */
//DATA STRUCTURE TO HOLD AORC FORCING DATA
struct aorc_forcing_data_cfe
{
// struct NAME                          DESCRIPTION                                            ORIGINAL AORC NAME
//____________________________________________________________________________________________________________________
    double precip_kg_per_m2;                // Surface precipitation "kg/m^2"                         | APCP_surface
    double incoming_longwave_W_per_m2 ;     // Downward Long-Wave Rad. Flux at 0m height, W/m^2       | DLWRF_surface
    double incoming_shortwave_W_per_m2;     // Downward Short-Wave Radiation Flux at 0m height, W/m^2 | DSWRF_surface
    double surface_pressure_Pa;             // Surface atmospheric pressure, Pa                       | PRES_surface
    double specific_humidity_2m_kg_per_kg;  // Specific Humidity at 2m height, kg/kg                  | SPFH_2maboveground
    double air_temperature_2m_K;            // Air temparture at 2m height, K                         | TMP_2maboveground
    double u_wind_speed_10m_m_per_s;        // U-component of Wind at 10m height, m/s                 | UGRD_10maboveground
    double v_wind_speed_10m_m_per_s;        // V-component of Wind at 10m height, m/s                 | VGRD_10maboveground
    double latitude;                        // degrees north of the equator.  Negative south          | latitude
    double longitude;                       // degrees east of prime meridian. Negative west          | longitude
    long int time; //TODO: type?           // seconds since 1970-01-01 00:00:00.0 0:00               | time
} ;
typedef struct aorc_forcing_data_cfe aorc_forcing_data_cfe;

struct vol_tracking_struct{
    double vol_sch_runoff;
    double vol_sch_infilt;
    double vol_to_soil;
    double vol_to_gw;
    double vol_soil_to_gw;
    double vol_soil_to_lat_flow;
    double volstart;
    double volout;
    double volin;
    double vol_from_gw;
    double vol_out_giuh;
    double vol_in_nash;
    double vol_out_nash;
    double vol_in_gw_start;
    double vol_soil_start;
};
typedef struct vol_tracking_struct vol_tracking_struct;

struct cfe_state_struct {

    // *************************************
    double timestep_rainfall_input_m;
    double soil_reservoir_storage_deficit_m;
    double infiltration_depth_m;
    double gw_reservoir_storage_deficit_m;
    double timestep_h;

    // ***********************************************************
    // ***************** Non-dynamic allocations *****************
    // ***********************************************************
    // These structs are defined in cfe.h.
    struct conceptual_reservoir soil_reservoir;
    struct conceptual_reservoir gw_reservoir;
    struct NWM_soil_parameters NWM_soil_params;
    struct evapotranspiration_structure et_struct;
    struct vol_tracking_struct vol_struct;

    // Epoch-based start time (BMI start time is considered 0.0)
    long epoch_start_time;
    int num_timesteps;
    int current_time_step;
    int time_step_size;
  
    // an integer to flag the forcings coming as "set_value" from BMI
    int is_forcing_from_bmi;

    char* forcing_file;

    double Schaake_adjusted_magic_constant_by_soil_type;
    int num_lateral_flow_nash_reservoirs;

    double K_lf;
    double K_nash;

    int num_giuh_ordinates;

    // ***********************************************************
    // ******************* Dynamic allocations *******************
    // ***********************************************************
    struct aorc_forcing_data_cfe aorc;
    double* forcing_data_precip_kg_per_m2;
    long* forcing_data_time;

    //result_fluxes* fluxes;

    double* giuh_ordinates;
    double* nash_storage;
    double* runoff_queue_m_per_timestep;

    // These are likely only single values, but should be allocated as pointers so the pointer can be returned
//    double* flux_overland_m;   NOT NEEDED, redundant with flux_Schaake_output_runoff_m
    double* flux_Schaake_output_runoff_m;
    double* flux_giuh_runoff_m;
    double* flux_nash_lateral_runoff_m;
    double* flux_from_deep_gw_to_chan_m;
    double* flux_perc_m;
    double* flux_lat_m;
    double* flux_Qout_m;

    int verbosity;

};
typedef struct cfe_state_struct cfe_state_struct;

extern double greg_2_jul(long year, long mon, long day, long h, long mi,
                         double se);
extern void calc_date_cfe(double jd, long *y, long *m, long *d, long *h, long *mi,
                      double *sec);

extern void itwo_alloc_cfe( int ***ptr, int x, int y);
extern void dtwo_alloc_cfe( double ***ptr, int x, int y);
extern void d_alloc_cfe(double **var,int size);
extern void i_alloc_cfe(int **var,int size);

extern void parse_aorc_line_cfe(char *theString,long *year,long *month, long *day,long *hour,
                            long *minute, double *dsec, struct aorc_forcing_data_cfe *aorc);

extern void get_word_cfe(char *theString,int *start,int *end,char *theWord,int *wordlen);

int read_init_config_cfe(const char* config_file, cfe_state_struct* model, double* alpha_fc, double* soil_storage,
                     int* is_soil_storage_ratio);

extern void init_soil_reservoir(cfe_state_struct* cfe_ptr, double alpha_fc, double max_storage, double storage,
                     int is_storage_ratios);

extern double init_reservoir_storage(int is_ratio, double amount, double max_amount);

extern void initialize_volume_trackers(cfe_state_struct* cfe_ptr);

extern void print_cfe_flux_header();
extern void print_cfe_flux_at_timestep(cfe_state_struct* cfe_ptr);
extern void mass_balance_check(cfe_state_struct* cfe_ptr);

// Bmi* register_bmi(Bmi *model);

Bmi* register_bmi_cfe(Bmi *model);

extern void run_cfe(cfe_state_struct* model);

cfe_state_struct * new_bmi_cfe(void);

#if defined(__cplusplus)
}
#endif

#endif
