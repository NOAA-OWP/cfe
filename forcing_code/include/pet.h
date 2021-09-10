#ifndef PET_H
#define PET_H

#include <stdio.h>
#include <math.h>
#include <string.h>

#define TRUE  1
#define FALSE 0

#define CP  1.006e+03  //  specific heat of air at constant pressure, J/(kg K), a physical constant.
#define KV2 0.1681     //  von Karman's constant squared, equal to 0.41 squared, unitless
#define TK  273.15     //  temperature in Kelvin at zero degree Celcius
#define SB  5.67e-08   //  stefan_boltzmann_constant in units of W/m^2/K^4

// NOTE: SET YOUR EDIT WINDOW TO 120 CHARACTER WIDTH TO READ THIS CODE IN ITS ENTIRETY.

//#####################################
// evapotranspiration (ET) module, 
// Version 1.0 by Fred L. Ogden, NOAA-NWS-OWP, May, 2020.
//
// includes five different methods to clculate ET, from Chow, Maidment & Mays Textbook, and UNFAO Penman-Monteith:
// 1. energy balance method
// 2. aerodynamic method
// 3. combination method, which combines 1 & 2.
// 4. Priestley-Taylor method, which assumes the ratio between 1 & 2, and only calculates 1.
// 5. Penman-Monteith method, which requires a value of canopy resistance term, and does not rely on 1 or 2.
// 
// This subroutine requires a considerable amount of meteorological data as input.
// a. temperature and (relative-humidity or specific humidity) and the heights at which they are measured.
// b. near surface wind speed measurement and the height at which it was measured.
// c. the ambient atmospheric temperature lapse rate
// d. the fraction of the sky covered by clouds
// e. (optional) the height above ground to the cloud base. If not provided, then assumed.
// f. the day of the year (1-366) and time of day (UTC only!)
// g. the skin temperature of the earth's surface, TODO: should come from another module to calc. soil or veg. temp.
// h. the zero-plane roughness height of the atmospheric boundary layer assuming log-law behavior (from land cover)
// i. the average root zone soil temperature, or near-surface water temperature in the case of lake evaporation.
// j. the incoming solar (shortwave) radiation.  If not provided it is computed from d,e,f, using an
//    updated method similar to the one presented in Bras, R.L. Hydrology.  Requires value of the Linke atmospheric
//    turbidity factor, which varies from 2 for clear mountain air to 5 for smoggy air.  According to Hove & Manyumbu
//    2012, who calculated values over Zimbabwe that varied from 2.14 to 3.71.  Other values exist in the literature.
//    TODO: This turbidity factor could be calculated from satellite obs. or maybe NOAA already does this?
//
// All radiation calculations needed for 1, 3, 4, and 5 require net radiation calculations at the land surface.
// the net radiation is calculated using a, c, d, e, f, g, j, plus the Linke turbidity factor, which can be estimate
// from satellite observations.
//

// NOTE THE VALUE OF evapotranspiration_params.zero_plane_displacement_height COMES FROM LAND COVER DATA.
// Taken from:    https://websites.pmc.ucsc.edu/~jnoble/wind/extrap/
//
//Roughness Roughness  Landscape Type
// Class    Length (m)	
//----------------------------------------------------------------------------------------------------------------
// 0         0.0002     Smooth water surface
// 0.2       0.0005     Inlet water
// 0.5       0.0024     Completely open terrain, smooth surface, e.g. concrete runways in airports, mowed grass, etc.
// 1         0.03       Open agricultural area without fences and hedgerows and very scattered buildings. Only softly
//                        rounded hills
// 1.5       0.055      Agricultural land with some houses and 8 metre tall sheltering hedgerows with a distance of 
//                        approximately 1250 metres
// 2         0.1        Agricultural land with some houses and 8 metre tall sheltering hedgerows with a distance of 
//                        approximately 500 metres
// 2.5       0.2        Agricultural land with many houses, shrubs and plants, or 8 metre tall sheltering hedgerows 
//                        with a distance of approximately 250 metres
// 3         0.4        Villages, small towns, agricultural land with many or tall sheltering hedgerows, forests 
//                        and very rough and uneven terrain
// 3.5       0.8        Larger cities with tall buildings
// 4         1.6        Very large cities with tall buildings and skyscrapers
// Roughness definitions according to the European Wind Atlas. 
//
// According to the UN FAO Penman-Monteith example here: 
//        http://www.fao.org/3/X0490E/x0490e06.htm#aerodynamic%20resistance%20(ra)
// the zero plane roughness length,"d" can be approximated as 2/3 of the vegetation height (H): d=2/3*H.
// the momentum roughness height "zom" can be estimated as 0.123*H.
// the heat transfer roughness height "zoh" can be approximated as 0.1 * zom.
//-----------------------------------------------------------------------------------------------------------------

//DATA STRUCTURE TO HOLD AORC FORCING DATA
struct aorc_forcing_data_pet
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
  double time; //TODO: type?           // seconds since 1970-01-01 00:00:00.0 0:00               | time
};
typedef struct aorc_forcing_data_pet aorc_forcing_data_pet;

struct pevapotranspiration_options  // these determine which method is applied to calculate ET.
{
  // element NAME                       DESCRIPTION
  //____________________________________________________________________________________________________________________
  int yes_aorc;                     // set to TRUE if forcing data come from AORC
  int shortwave_radiation_provided; // set to TRUE if yes_aorc==TRUE or if solar radiation inputs are provided

  // NOTE: these following options are exclusive.  Only one should be true for a particular catchment.
  int use_energy_balance_method;    // set to TRUE if using just the energy balance method for calculating PET
  int use_aerodynamic_method;       // set to TRUE if using just the aerodynamic method for calculating PET
  int use_combination_method;       // set to TRUE if using just the combination method for calculating PET
  int use_priestley_taylor_method;  // set to TRUE if using just the Priestley-Taylor method for calculating PET
  int use_penman_monteith_method;   // set to TRUE if using just the Penman Monteith method for calculating PET
};

struct pevapotranspiration_params
{
  // element NAME                                    DESCRIPTION                                                
  //____________________________________________________________________________________________________________________
  double wind_speed_measurement_height_m;          // set to 0.0 if unknown, will default =2.0 [m] 
  double humidity_measurement_height_m;            // set to 0.0 if unknown, will default =2.0 [m]
  double vegetation_height_m;                      // TODO this should come from land cover data and a veg. height model
  double zero_plane_displacement_height_m;         // depends on surface roughness [m],
  double momentum_transfer_roughness_length_m;     // poorly defined.  If unknown, pass down 0.0 as a default
  double heat_transfer_roughness_length_m;         // poorly defined.  If unknown, pass down 0.0 as a default
  double latitude;                                 // could be used to adjust canopy resistance for seasonality
  double longitude;                                // could be used to adjust canopy resistance for seasonality
  int    day_of_year;                              // could be used to adjust canopy resistance for seasonality
};

struct pevapotranspiration_forcing
{
  // element NAME                          DESCRIPTION                                                                  
  //___________________________________________________________________________________________________________________
  double net_radiation_W_per_sq_m;       // NOTE: ground heat flux is subtracted out in calculation subroutine
  double air_temperature_C;
  double relative_humidity_percent;      // this and specific_humidity_2m_kg_per_kg are redundant, so-  
  double specific_humidity_2m_kg_per_kg; // specify the missing one using a negative number, the other will be used
  double air_pressure_Pa;  
  double wind_speed_m_per_s;             // this is the value measured 2 m above the canopy
  double canopy_resistance_sec_per_m;    // depends on vegetation type and point in growing season 
  double water_temperature_C;            // if >100, assume 15 C.  Used to calculate latent heat of vaporization
  double ground_heat_flux_W_per_sq_m;    // from a model or assumed =0.  typically small on average
};

struct surface_radiation_params
{
  // element NAME                          DESCRIPTION                                                                  
  //____________________________________________________________________________________________________________________
  double surface_longwave_emissivity;  // dimensionless (0-1) < 1.0 for water and snow, all other surfaces = 1.0
  double surface_shortwave_albedo;     // dimensionless (0-1) from land cover.  Dynamic seasonally
};

struct surface_radiation_forcing
{
  // element NAME                          DESCRIPTION                                                                  
  //____________________________________________________________________________________________________________________
  double incoming_shortwave_radiation_W_per_sq_m;  // TODO could be calculated if unavailable, but not now.
  double incoming_longwave_radiation_W_per_sq_m;  // set to a large negative number if unknown (e.g. -1.0e-5)
  double air_temperature_C;            // usually value at 2.0 m., maybe 30 m if from WRF
  double relative_humidity_percent;    // usually value at 2.0 m., maybe 30 m if from WRF
  double surface_skin_temperature_C;   // from a model or assumed...  smartly.  could be soil/rock, veg., snow, water
  double ambient_temperature_lapse_rate_deg_C_per_km;  // This is a standard WRF output.  Typ. 6.49 K/km ICAO std. atm.
  double cloud_cover_fraction;         // dimensionless (0-1).  This should be a WRF output.
  double cloud_base_height_m;          // the height from ground to bottom of clouds in m.  From WRF output
  double atmospheric_turbidity_factor; // Linke turbidity factor needed iff pet_options.shortwave_radiation_provided=FALSE
  int    day_of_year;
  double zulu_time;                    // (0.0-23.999999) hours
};

struct solar_radiation_options
{
  // element NAME                          DESCRIPTION
  //____________________________________________________________________________________________________________________
  int cloud_base_height_known;   // set this to TRUE to use the default values from the Bras textbook.
};

struct solar_radiation_parameters
{
  // element NAME                          DESCRIPTION                                                                  
  //____________________________________________________________________________________________________________________
  double latitude_degrees;       // positive north of the equator, negative south
  double longitude_degrees;      // negative west of prime meridian, positive east
  double site_elevation_m;       // elevation of the observer, m
};

struct solar_radiation_results
{
  // element NAME                          DESCRIPTION                                                                  
  //____________________________________________________________________________________________________________________
  double solar_radiation_flux_W_per_sq_m;            // on a plane perpendicular to the earth-sun line
  double solar_radiation_horizontal_flux_W_per_sq_m; // on a horizontal plane tangent to earth
  double solar_radiation_cloudy_flux_W_per_sq_m;     // on a plane perpendicular to the earth-sun line
  double solar_radiation_horizontal_cloudy_flux_W_per_sq_m; // on a horizontal plane tangent to earth considering clouds
  double solar_elevation_angle_degrees;              // height of the sun above (+) or below (-) horizon, degrees.
  double solar_azimuth_angle_degrees;                // azimuth pointing towards the sun, degrees (0-360)
  double solar_local_hour_angle_degrees;             // local hour angle (deg.) to the sun, negative=a.m., positive=p.m.
};

struct intermediate_vars
{
  // element NAME                       DESCRIPTION
  //____________________________________________________________________________________________________________________
  double liquid_water_density_kg_per_m3;       // rho_w
  double water_latent_heat_of_vaporization_J_per_kg;    // eqn 2.7.6 Chow etal., // aka 'lambda'
  double air_saturation_vapor_pressure_Pa;
  double air_actual_vapor_pressure_Pa;
  double vapor_pressure_deficit_Pa;            // VPD
  double moist_air_gas_constant_J_per_kg_K;    // R_a
  double moist_air_density_kg_per_m3;          // rho_a
  double slope_sat_vap_press_curve_Pa_s;       // delta
  //double water_latent_heat_of_vaporization_J_per_kg;
  double psychrometric_constant_Pa_per_C;      // gamma
};
struct bmi
{
  /*    
      JMFRAME: Fred suggested changing time step size to make the units of seconds explicit...
               Now we can only run PET on a timestep of seconds.
               But what if we want to run at a timestep different than seconds? 
               We probably would never do that. The output is in units of m_per_s, so this is probably good.
  */
  int time_step_size_s;
  long int num_timesteps;
  double current_time_step;   // this is the actual time of the run.
  long int current_step;        // this is a sequential value to find the correct row from forcing file
  double current_time;        // this should be in "Seconds since 1970", so should be start time plus current time step
  int verbose;
  int run_unit_tests;
  int is_forcing_from_bmi;
};

struct pet_model{
  
  // FLAGS
  int yes_aorc; // if TRUE then using AORC forcing data- if FALSE then we must calculate incoming short/longwave rad.
  int yes_wrf;  // if TRUE then we get radiation winds etc. from WRF output.  TODO not implemented.
  int pet_method;
  double pet_m_per_s;
  char* forcing_file;
  // ***********************************************************
  // ******************* Dynamic allocations *******************
  // ***********************************************************
  //aorc_forcing_data* forcings;
  double* forcing_data_precip_kg_per_m2;
  double* forcing_data_surface_pressure_Pa;
  double* forcing_data_time;
  double* forcing_data_incoming_longwave_W_per_m2 ;     // Downward Long-Wave Rad. Flux at 0m height, W/m^2       | DLWRF_surface
  double* forcing_data_incoming_shortwave_W_per_m2;     // Downward Short-Wave Radiation Flux at 0m height, W/m^2 | DSWRF_surface
  double* forcing_data_specific_humidity_2m_kg_per_kg;  // Specific Humidity at 2m height, kg/kg                  | SPFH_2maboveground
  double* forcing_data_air_temperature_2m_K;            // Air temparture at 2m height, K                         | TMP_2maboveground
  double* forcing_data_u_wind_speed_10m_m_per_s;        // U-component of Wind at 10m height, m/s                 | UGRD_10maboveground
  double* forcing_data_v_wind_speed_10m_m_per_s;        // V-component of Wind at 10m height, m/s                 | VGRD_10maboveground

  struct aorc_forcing_data_pet aorc;

  struct pevapotranspiration_options pet_options;
  struct pevapotranspiration_params  pet_params;
  struct pevapotranspiration_forcing pet_forcing;
  struct intermediate_vars inter_vars;

  struct surface_radiation_params   surf_rad_params;
  struct surface_radiation_forcing  surf_rad_forcing;

  struct solar_radiation_options    solar_options;
  struct solar_radiation_parameters solar_params;
  struct solar_radiation_results    solar_results;

  struct bmi bmi;

};
typedef struct pet_model pet_model;

extern void alloc_pet_model(pet_model *model);

extern void free_pet_model(pet_model *model);

extern int run_pet(pet_model* model);

void pet_setup(pet_model* model);
void pet_unit_tests(pet_model* model);

/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS JUST UTILITY MEMORY AND TIME FUNCTION CODE */
/**************************************************************************/
extern void parse_aorc_line_pet(char *theString,long *year,long *month, long *day,long *hour,
                            long *minute, double *dsec, struct aorc_forcing_data_pet *aorc);

extern void get_word_pet(char *theString,int *start,int *end,char *theWord,int *wordlen);
extern void itwo_alloc_pet( int ***ptr, int x, int y);
extern void dtwo_alloc_pet( double ***ptr, int x, int y);
extern void d_alloc_pet(double **var,int size);
extern void i_alloc(int **var,int size);
extern double greg_2_jul_pet(long year, long mon, long day, long h, long mi,
                         double se);
extern void calc_date_pet(double jd, long *y, long *m, long *d, long *h, long *mi,
                      double *sec);
#endif