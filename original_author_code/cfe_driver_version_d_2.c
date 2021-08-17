#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "cfe.c"

#define TRUE 1
#define FALSE 0
#define DEBUG 1
#define MAX_NUM_GIUH_ORDINATES 10
#define MAX_NUM_NASH_CASCADE    3
#define MAX_NUM_RAIN_DATA 720

// t-shirt approximation of the hydrologic routing funtionality of the National Water Model v 1.2, 2.0, and 2.1
// This code was developed to test the hypothesis that the National Water Model runoff generation, vadose zone
// dynamics, and conceptual groundwater model can be greatly simplified by acknowledging that it is truly a 
// conceptual model.  The hypothesis is supported by a number of observations made during a 2017-2018 deep dive
// into the NWM code.  Thesed are:
//
// 1. Rainfall/throughfall/melt partitioning in the NWM is based on a simple curve-number like approach that
//    was developed by Schaake et al. (1996) and which is very similar to the Probability Distributed Moisture (PDM)
//    function by Moore, 1985.   The Schaake function is a single valued function of soil moisture deficit,
//    predicts 100% runoff when the soil is saturated, like the curve-number method, and is fundamentally simple.
// 2. Run-on infiltration is strictly not calculated.  Overland flow routing applies the Schaake function repeatedly
//    to predict this phenomenon, which violates the underlying assumption of the PDM method that only rainfall 
//    inputs affect soil moisture.
// 3. The water-content based Richards' equation, applied using a coarse-discretization, can be replaced with a simple
//    conceptual reservoir because it never allows saturation or infiltration-excess runoff unless deactivated by
//    assuming no-flow lower boundary condition.  Since this form of Richards' equation cannot simulate heterogeneous
//    soil layers, it can be replaced with a conceptual reservoir.
// 4. The lateral flow routing function in the NWM is purely conceptual.  It is activated whenever the soil water
//    content in one or more of the four Richards-equation discretizations reaches the wilting point water content.
//    This activation threshold is physically unrealistic, because in most soils lateral subsurface flow is not
//    active until pore water pressures become positive at some point in the soil profile.  Furthermore, the lateral
//    flow hydraulic conductivity is assumed to be the vertical hydraulic conductivity multiplied by a calibration
//    factor "LKSATFAC" which is allowed to vary between 10 and 10,000 during calibration, resulting in an anisotropy
//    ratio that varies over the same range, without correlation with physiographic characteristics or other support.
//
//    This code implements these assumptions using pure conceptualizations.  The formulation consists of the following:
//
//    1. Rainfall is partitioned into direct runoff and soil moisture using the Schaake function.
//    2. Rainfall that becomes direct runoff is routed to the catchment outlet using a geomorphological instantanteous
//       unit hydrograph (GIUH) approach, eliminating the 250 m NWM routing grid, and the incorrect use of the Schaake
//       function to simulate run-on infiltration.
//    3. Water partitioned by the Schaake function to be soil moisture is placed into a conceptual linear reservoir
//       that consists of two outlets that apply a minimum storage activation threshold.   This activation threshold
//       is identical for both outlets, and is based on an integral solution of the storage in the soil assuming
//       Clapp-Hornberger parameters equal to those used in the NWM to determine that storage corresponding to a
//       soil water content 0.5 m above the soil column bottom that produces a soil suction head equal to -1/3 atm,
//       which is a commonly applied assumption used to estimate the field capacity water content.
//       The first outlet calculates vertical percolation of water to deep groundwater using the saturated hydraulic
//       conductivity of the soil multiplied by the NWM "slope" parameter, which when 1.0 indicates free drainage and
//       when 0.0 indicates a no-flow lower boundary condition.   The second outlet is used to calculate the flux to
//       the soil lateral flow path, using a conceptual LKSATFAC-like calibration parameter.
//    4. The lateral flow is routed to the catchment outlet using a Nash-cascade of reservoirs to produce a mass-
//       conserving delayed response, and elminates the need for the 250 m lateral flow routing grid.
//    5. The groundwater contribution to base flow is modeled using either (a) an exponential nonlinear reservoir
//       identical to the one in the NWM formulation, or (b) a nonlinear reservoir forumulation, which can also be
//       made linear by assuming an exponent value equal to 1.0.
//
//    This code was written entirely by Fred L. Ogden, May 22-24, 2020, in the service of the NOAA-NWS Office of Water
//    Prediction, in Tuscaloosa, Alabama.


// define data structures
//--------------------------

//struct conceptual_reservoir
//{
//// this data structure describes a nonlinear reservoir having two outlets, one primary with an activation
//// threshold that may be zero, and a secondary outlet with a threshold that may be zero
//// this will also simulate a linear reservoir by setting the exponent parameter to 1.0 iff is_exponential==FALSE
//// iff is_exponential==TRUE, then it uses the exponential discharge function from the NWM V2.0 forumulation
//// as the primary discharge with a zero threshold, and does not calculate a secondary discharge.
////--------------------------------------------------------------------------------------------------
//int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
//double storage_max_m;   // maximum storage in this reservoir
//double storage_m;       // state variable.
//double coeff_primary;    // the primary outlet
//double exponent_primary;
//double storage_threshold_primary_m;
//double storage_threshold_secondary_m;
//double coeff_secondary;
//double exponent_secondary;
//};
//
//struct NWM_soil_parameters
//{
//// using same variable names as used in NWM.  <sorry>
//double smcmax;  // effective porosity [V/V]
//double wltsmc;  // wilting point soil moisture content [V/V]
//double satdk;   // saturated hydraulic conductivity [m s-1]
//double satpsi;	// saturated capillary head [m]
//double bb;      // beta exponent on Clapp-Hornberger (1978) soil water relations [-]
//double mult;    // the multiplier applied to satdk to route water rapidly downslope
//double slop;   // this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow.
//double D;       // soil depth [m]
//};


struct conceptual_reservoir soil_reservoir;
struct conceptual_reservoir gw_reservoir;
struct NWM_soil_parameters NWM_soil_params;


//DATA STRUCTURE TO HOLD AORC FORCING DATA
struct aorc_forcing_data
{
// struct NAME                          DESCRIPTION                                            ORIGINAL AORC NAME     
//____________________________________________________________________________________________________________________
float precip_kg_per_m2;                // Surface precipitation "kg/m^2"                         | APCP_surface
float incoming_longwave_W_per_m2 ;     // Downward Long-Wave Rad. Flux at 0m height, W/m^2       | DLWRF_surface
float incoming_shortwave_W_per_m2;     // Downward Short-Wave Radiation Flux at 0m height, W/m^2 | DSWRF_surface
float surface_pressure_Pa;             // Surface atmospheric pressure, Pa                       | PRES_surface
float specific_humidity_2m_kg_per_kg;  // Specific Humidity at 2m height, kg/kg                  | SPFH_2maboveground
float air_temperature_2m_K;            // Air temparture at 2m height, K                         | TMP_2maboveground
float u_wind_speed_10m_m_per_s;        // U-component of Wind at 10m height, m/s                 | UGRD_10maboveground
float v_wind_speed_10m_m_per_s;        // V-component of Wind at 10m height, m/s                 | VGRD_10maboveground
float latitude;                        // degrees north of the equator.  Negative south          | latitude
float longitude;                       // degrees east of prime meridian. Negative west          | longitude
long int time; //TODO: type?           // seconds since 1970-01-01 00:00:00.0 0:00               | time
} ;

extern double greg_2_jul(long year, long mon, long day, long h, long mi,
                         double se);
extern void calc_date(double jd, long *y, long *m, long *d, long *h, long *mi,
               double *sec);

extern void itwo_alloc( int ***ptr, int x, int y);
extern void dtwo_alloc( double ***ptr, int x, int y);
extern void d_alloc(double **var,int size);
extern void i_alloc(int **var,int size);
extern void l_alloc(long **var,int size);

extern void parse_aorc_line(char *theString,long *year,long *month, long *day,long *hour,
                long *minute, double *dsec, struct aorc_forcing_data *aorc);

extern void get_word(char *theString,int *start,int *end,char *theWord,int *wordlen);

//####################################################################
//####################################################################
//####################  MAIN PROGRAM   ###############################
//####################################################################
//####################################################################

int main()
{
FILE *in_fptr;
FILE *out_fptr;
FILE *out_debug_fptr;

//local variables
//-------------------------------------------------------------------
int i;
int tstep;
double upper_lim;
double lower_lim;
char theString[513];   // dangerously hard coded string size... TODO fixme.
long year,month,day,hour,minute;
double dsec;
double jdate_start=0.0;

double catchment_area_km2=15.617167355002097;            // in the range of our desired size
double drainage_density_km_per_km2=3.5;   // this is approx. the average blue line drainage density for CONUS

int    num_timesteps;
int    num_rain_dat;
double timestep_h;
long *sim_year=NULL;
long *sim_month=NULL;
long *sim_day=NULL;
long *sim_hour=NULL;
long *sim_minute=NULL;
long *sim_second=NULL;

// forcing
double *rain_rate=NULL;
double timestep_rainfall_input_m;
int yes_aorc=TRUE;                  // change to TRUE to read in entire AORC precip/met. forcing data set.

//Schaake partitioning function parameters

double Schaake_adjusted_magic_constant_by_soil_type;
double Schaake_output_runoff_m;
double infiltration_depth_m;

// GIUH state & parameters

double *giuh_ordinates = NULL;                  // assumed GIUH hydrograph ordinates
double *runoff_queue_m_per_timestep =NULL;
int    num_giuh_ordinates;
double giuh_runoff_m;
double soil_reservoir_storage_deficit_m;        // the available space in the soil conceptual reservoir

// groundwater storage parameters and state

// lateral flow function parameters
double lateral_flow_threshold_storage_m;
double field_capacity_storage_threshold_m;
double lateral_flow_linear_reservoir_constant;
int num_lateral_flow_nash_reservoirs;
double K_nash;  // lateral_flow_nash_cascade_reservoir_constant;
double *nash_storage = NULL; 
double assumed_near_channel_water_table_slope; // [L/L]
double nash_lateral_runoff_m;

// calculated flux variables
double flux_overland_m;                // flux of surface runoff that goes through the GIUH convolution process
double flux_perc_m=0.0;                // flux from soil to deeper groundwater reservoir
double flux_lat_m;                     // lateral flux in the subsurface to the Nash cascade
double flux_from_deep_gw_to_chan_m;    // flux from the deep reservoir into the channels
double gw_reservoir_storage_deficit_m; // the available space in the conceptual groundwater reservoir
//double primary_flux,secondary_flux;    // temporary vars.    //jmframe: move to CFE for conceptual_reservoir_flux_calc

double field_capacity_atm_press_fraction; // [-]
double soil_water_content_at_field_capacity; // [V/V]
double atm_press_Pa;
double unit_weight_water_N_per_m3; 
double H_water_table_m;  // Hwt in NWM/t-shirt parameter equiv. doc  discussed between Eqn's 4 and 5.
                         // this is the distance down to a fictitious water table from the point there the 
                         // soil water content is assumed equal to field capacity at the middle of lowest discretization
double trigger_z_m;      // this is the distance up from the bottom of the soil that triggers lateral flow when 
                         // the soil water content equals field capacity.   0.5 for center of bottom discretization
double Omega;            // The limits of integration used in Eqn. 5 of the parameter equivalence document.

double Qout_m;           // the sum of giuh, nash-cascade, and base flow outputs m per timestep

struct aorc_forcing_data aorc_data;

double refkdt=3.0;   // in Sugar Creek WRF-Hydro calibration this value was 0.15  3.0 is the default from noah-mp.

//mass balance variables.  These all store cumulative discharges per unit catchment area [m3/m2] or [m]
//-----------------------

double volstart =0.0;
double vol_sch_runoff   =0.0;
double vol_sch_infilt   =0.0;

double vol_out_giuh     =0.0;
double vol_end_giuh     =0.0;

double vol_to_gw       =0.0;
double vol_in_gw_start =0.0;
double vol_in_gw_end   =0.0;
double vol_from_gw     =0.0;

double vol_in_nash      =0.0;
double vol_in_nash_end  =0.0;  // note the nash cascade is empty at start of simulation.
double vol_out_nash     =0.0;

double vol_soil_start       =0.0;
double vol_to_soil          =0.0;
double vol_soil_to_lat_flow =0.0;
double vol_soil_to_gw       =0.0;  // this should equal vol_to_gw
double vol_soil_end         =0.0;
// note, vol_from_soil_to_gw is same as vol_to_gw.

double volin    =0.0;
double volout   =0.0;
double volend   =0.0;


if((out_fptr=fopen("test.out","w"))==NULL)
  {printf("Can't open output file\n");exit(0);}

#ifdef DEBUG
        if((out_debug_fptr=fopen("debug.out","w"))==NULL)
                {printf("Can't open output file\n");exit(0);}
#endif

d_alloc(&giuh_ordinates, MAX_NUM_GIUH_ORDINATES);  // allocate memory to store the GIUH ordinates 
d_alloc(&runoff_queue_m_per_timestep, MAX_NUM_GIUH_ORDINATES+1); // allocate memory to store convolution queue

// catchment properties
//------------------------
catchment_area_km2=15.617167355002097;

//initialize simulation constants
//---------------------------
num_timesteps=10000;
timestep_h=1.0;
atm_press_Pa=101325.0;
unit_weight_water_N_per_m3=9810.0;

//initialize NWM soil parameters, using LOAM soils from SOILPARM.TBL.
//--------------------------------------------------------------------
NWM_soil_params.smcmax=0.439;   // [V/V]
NWM_soil_params.wltsmc=0.066;   // [V/V]
NWM_soil_params.satdk=3.38e-06; // [m per sec.] 
NWM_soil_params.satpsi=0.355;   // [m]
NWM_soil_params.bb=4.05;        // [-]    -- called "bexp" in NWM
//NWM_soil_params.slop=1.0;    // calibration factor that varies from 0 (no flow) to 1 (free drainage) at soil bottom
NWM_soil_params.slop=0.01;    // calibration factor that varies from 0 (no flow) to 1 (free drainage) at soil bottom
NWM_soil_params.D=2.0;       // [m] soil thickness assumed in the NWM not from SOILPARM.TBL
NWM_soil_params.mult=1000.0;    // not from SOILPARM.TBL, This is actually calibration parameter: LKSATFAC

trigger_z_m=0.5;   // distance from the bottom of the soil column to the center of the lowest discretization

// calculate the activation storage ffor the secondary lateral flow outlet in the soil nonlinear reservoir.
// following the method in the NWM/t-shirt parameter equivalence document, assuming field capacity soil
// suction pressure = 1/3 atm= field_capacity_atm_press_fraction * atm_press_Pa.

field_capacity_atm_press_fraction=0.33;  //alpha in Eqn. 3.

// equation 3 from NWM/t-shirt parameter equivalence document
H_water_table_m=field_capacity_atm_press_fraction*atm_press_Pa/unit_weight_water_N_per_m3; 
soil_water_content_at_field_capacity=NWM_soil_params.smcmax*
                                     pow(H_water_table_m/NWM_soil_params.satpsi,(1.0/NWM_soil_params.bb));

// solve the integral given by Eqn. 5 in the parameter equivalence document.
// this equation calculates the amount of water stored in the 2 m thick soil column when the water content 
// at the center of the bottom discretization (trigger_z_m) is at field capacity
Omega=H_water_table_m-trigger_z_m;
lower_lim= pow(Omega,(1.0-1.0/NWM_soil_params.bb))/(1.0-1.0/NWM_soil_params.bb);
upper_lim= pow(Omega+NWM_soil_params.D,(1.0-1.0/NWM_soil_params.bb))/(1.0-1.0/NWM_soil_params.bb);
field_capacity_storage_threshold_m=NWM_soil_params.smcmax*pow(1.0/NWM_soil_params.satpsi,(-1.0/NWM_soil_params.bb))*
                                (upper_lim-lower_lim);
    
printf("field capacity storage threshold = %lf m\n", field_capacity_storage_threshold_m);


// initialize giuh parameters.  These would come from another source.
//------------------------------------------------------------------
num_giuh_ordinates=7;
if(num_giuh_ordinates>MAX_NUM_GIUH_ORDINATES)
  {
  fprintf(stderr,"Big problem, number of giuh ordinates greater than MAX_NUM_GIUH_ORDINATES.  Stopped.\n");
  exit(0);
  }
giuh_ordinates[0]=0.06;  // note these sum to 1.0.  If we have N ordinates, we need a queue sized N+1 to perform
giuh_ordinates[1]=0.51;  // the convolution.
giuh_ordinates[2]=0.28;
giuh_ordinates[3]=0.12;
giuh_ordinates[4]=0.03;

// initialize Schaake parameters
//--------------------------------
// in noah-mp refkdt=3.0. 
Schaake_adjusted_magic_constant_by_soil_type=refkdt*NWM_soil_params.satdk/2.0e-06; // 2.0e-06 is used in noah-mp


// initialize lateral flow function parameters
//---------------------------------------------
assumed_near_channel_water_table_slope=0.01; // [L/L]
lateral_flow_threshold_storage_m=field_capacity_storage_threshold_m;  // making them the same, but they don't have 2B

// Equation 10 in parameter equivalence document.
lateral_flow_linear_reservoir_constant=2.0*assumed_near_channel_water_table_slope*NWM_soil_params.mult*
                                       NWM_soil_params.satdk*NWM_soil_params.D*drainage_density_km_per_km2;   // m/s
lateral_flow_linear_reservoir_constant*=3600.0;  // convert to m/h

num_lateral_flow_nash_reservoirs=2;
if(num_lateral_flow_nash_reservoirs>MAX_NUM_NASH_CASCADE)
  {
  fprintf(stdout,"Number of Nash Cascade linear reservoirs greater than MAX_NUM_NASH_CASCADE.  Stopped.\n");
  return(-2);
  }
K_nash=0.03;                                          // lateral_flow_nash_cascade_reservoir_constant. TODO Calibrate.

d_alloc(&nash_storage,MAX_NUM_NASH_CASCADE);
for(i=0;i<MAX_NUM_NASH_CASCADE;i++) nash_storage[i]=0.0;  /* initialize */


//###################################################################################
//###########   INITIAL GW AND SOIL RESERVOIR PARAMS AND STATE CONDITIONS ###########
//###################################################################################


//  Populate the groundwater conceptual reservoir data structure
//-----------------------------------------------------------------------
// one outlet, 0.0 threshold, nonliner and exponential as in NWM
gw_reservoir.is_exponential=TRUE;         // set this true TRUE to use the exponential form of the discharge equation
gw_reservoir.storage_max_m=1.0;            // calibrated Sugar Creek WRF-Hydro value 16.0, I assume mm.
gw_reservoir.coeff_primary=0.01;           // per h
gw_reservoir.exponent_primary=6.0;              // linear iff 1.0, non-linear iff > 1.0
gw_reservoir.storage_threshold_primary_m=0.0;     // 0.0 means no threshold applied
gw_reservoir.storage_threshold_secondary_m=0.0;   // 0.0 means no threshold applied
gw_reservoir.coeff_secondary=0.0;                 // 0.0 means that secondary outlet is not applied
gw_reservoir.exponent_secondary=1.0;              // linear

gw_reservoir.storage_m=gw_reservoir.storage_max_m*0.5;  // INITIALIZE HALF FULL.
//-------------------------------------------------------------------------------------------------------------

volstart       += gw_reservoir.storage_m;    // initial mass balance checks in g.w. reservoir
vol_in_gw_start = gw_reservoir.storage_m;  

// Initialize the soil conceptual reservoir data structure.  Indented here to highlight different purposes
//-------------------------------------------------------------------------------------------------------------
// soil conceptual reservoir first, two outlets, two thresholds, linear (exponent=1.0).
soil_reservoir.is_exponential=FALSE;  // set this true TRUE to use the exponential form of the discharge equation
                                      // this should NEVER be set to true in the soil reservoir.
soil_reservoir.storage_max_m=NWM_soil_params.smcmax*NWM_soil_params.D;
  //  vertical percolation parameters------------------------------------------------
    soil_reservoir.coeff_primary=NWM_soil_params.satdk*NWM_soil_params.slop*3600.0; // m per h
    soil_reservoir.exponent_primary=1.0;      // 1.0=linear
    soil_reservoir.storage_threshold_primary_m=field_capacity_storage_threshold_m;
  // lateral flow parameters --------------------------------------------------------
    soil_reservoir.coeff_secondary=0.01;  // 0.0 to deactiv. else =lateral_flow_linear_reservoir_constant;   // m per h
    soil_reservoir.exponent_secondary=1.0;   // 1.0=linear
    soil_reservoir.storage_threshold_secondary_m=lateral_flow_threshold_storage_m;
    
soil_reservoir.storage_m=soil_reservoir.storage_max_m*0.667;  // INITIALIZE SOIL STORAGE
//-------------------------------------------------------------------------------------------------------------

volstart          += soil_reservoir.storage_m;    // initial mass balance checks in soil reservoir
vol_soil_start     = soil_reservoir.storage_m;


d_alloc(&rain_rate,MAX_NUM_RAIN_DATA);
l_alloc(&sim_year,MAX_NUM_RAIN_DATA);
l_alloc(&sim_month,MAX_NUM_RAIN_DATA);
l_alloc(&sim_day,MAX_NUM_RAIN_DATA);
l_alloc(&sim_hour,MAX_NUM_RAIN_DATA);
l_alloc(&sim_minute,MAX_NUM_RAIN_DATA);
l_alloc(&sim_second,MAX_NUM_RAIN_DATA);

//########################################################################################################
// read in the aorc forcing data, just save the rain rate in mm/h in an array, ignore other vars ffor now.
//########################################################################################################
if(yes_aorc==TRUE)  // reading a csv file containing all the AORC meteorological/radiation and rainfall
  {
  if((in_fptr=fopen("cat58_01Dec2015.csv","r"))==NULL)
    {printf("Can't open input file\n");exit(0);}

  num_rain_dat=720;  // hard coded number of lines in the forcing input file.
  fgets(theString,512,in_fptr);  // read in the header.
  for(i=0;i<num_rain_dat;i++)
    {
    fgets(theString,512,in_fptr);  // read in a line of AORC data.
    parse_aorc_line(theString,&year,&month,&day,&hour,&minute,&dsec,&aorc_data);
    sim_year[i]=year;
    sim_month[i]=month;
    sim_day[i]=day;
    sim_hour[i]=hour;
    sim_minute[i]=minute;
    sim_second[i]=(long)dsec;

    if(i==0) jdate_start=greg_2_jul(year,month,day,hour,minute,dsec);      // calc & save the starting julian date of the rainfall data
    // saving only the rainfall data ffor now.
    // rain_rate[i]=(double)aorc_data.precip_kg_per_m2;  // assumed 1000 kg/m3 density of water.  This result is mm/h;
    rain_rate[i]=(double)aorc_data.precip_kg_per_m2*3600.0;  // precip_kg_per_m2 is in m/s, convert to m/h
    //if ((i > 12) && (i < 24)) rain_rate[i] = 0.0;
    }
  }
else  // reading a single column txt file that contains only rainfall (mm/h) 
  {
  if((in_fptr=fopen("cat58_01Dec2015_rain_only.dat","r"))==NULL)
    {printf("Can't open input file\n");exit(0);}

  num_rain_dat=720;  // hard coded number of lines in the forcing input file.
  for(i=0;i<num_rain_dat;i++)
    {
    fscanf(in_fptr,"%lf",&rain_rate[i]);
    }
  }  
fclose(in_fptr);
//######################  END OF READ FORCING CODE  ######################################################


//fprintf(out_fptr,"#    ,           hourly,    direct,   giuh,   lateral,  base,   total\n");
//fprintf(out_fptr,"#Time,           rainfall,  runoff,   runoff, flow,     flow,   discharge\n");
//fprintf(out_fptr,"# (h),           (m/h),     (m/h),    (m/h),  (m/h),    (m/h),  (m/h)\n");
fprintf(out_fptr,"Time Step,Time,Rainfall,Direct Runoff,GIUH Runoff,Lateral Flow,Base Flow,Total Discharge,Flow\n");
//tstep,Schaake_output_runoff_m,giuh_runoff_m,
//                                         nash_lateral_runoff_m, flux_from_deep_gw_to_chan_m, Qout_m 

// run the model for num_timesteps+500
//---------------------------------------

double lateral_flux;      // flux from soil to lateral flow Nash cascade +to cascade
double percolation_flux;  // flux from soil to gw nonlinear researvoir, +downward
double oldval;

//###################################################################################
//############################      TIME LOOP       #################################
//###################################################################################

num_timesteps=num_rain_dat+279;  // run a little bit beyond the rain data to see what happens.
//num_timesteps=num_rain_dat;
for(tstep=0;tstep<num_timesteps;tstep++)
  {
  
  
  // if(tstep<num_rain_dat)  timestep_rainfall_input_m=rain_rate[tstep]/1000.0;  // convert from mm/h to m w/ 1h timestep
  if(tstep<num_rain_dat)  timestep_rainfall_input_m=rain_rate[tstep];  // in m/h already, no conversion needed
  else                    timestep_rainfall_input_m=0.0;
  volin+=timestep_rainfall_input_m;


  //##################################################
  // run the cfe code
  //##################################################
cfe(
        &soil_reservoir_storage_deficit_m,
        NWM_soil_params,
        soil_reservoir,
        timestep_h,
        Schaake_adjusted_magic_constant_by_soil_type,
        timestep_rainfall_input_m,
        &Schaake_output_runoff_m,
        &infiltration_depth_m,
        &flux_overland_m,
        &vol_sch_runoff,
        &vol_sch_infilt,
        &flux_perc_m,
        &vol_to_soil,
        &percolation_flux,
        &lateral_flux,
        &flux_lat_m,
        &gw_reservoir_storage_deficit_m,
        gw_reservoir,
        &vol_to_gw,
        &vol_soil_to_gw,
        &vol_soil_to_lat_flow,
        &volout,
        &flux_from_deep_gw_to_chan_m,
        &vol_from_gw,
        &giuh_runoff_m,
        num_giuh_ordinates,
        giuh_ordinates,
        runoff_queue_m_per_timestep,
        &vol_out_giuh,
        &nash_lateral_runoff_m,
        num_lateral_flow_nash_reservoirs,
        K_nash,
        nash_storage,
        &vol_in_nash,
        &vol_out_nash,
        &Qout_m
    );                      // #######################################################################


  // <<<<<<<<NOTE>>>>>>>>>
  // These fluxs are all in units of meters per time step.   Multiply them by the "catchment_area_km2" variable
  // to convert them into cubic meters per time step.
  /*
  timestep_rainfall_input_m *= catchment_area_km2*1000000.0;
  Schaake_output_runoff_m *= catchment_area_km2*1000000.0;
  giuh_runoff_m *= catchment_area_km2*1000000.0;
  nash_lateral_runoff_m *= catchment_area_km2*1000000.0;
  flux_from_deep_gw_to_chan_m *= catchment_area_km2*1000000.0;
  Qout_m *= catchment_area_km2*1000000.0;

  // WRITE OUTPUTS IN cubic meter/s
  fprintf(out_fptr,"%16.8lf %lf %lf %lf %lf %lf %lf\n",jdate_start+(double)tstep*1.0/24.0,
                           timestep_rainfall_input_m/3600.0,Schaake_output_runoff_m/3600.0,
                           giuh_runoff_m/3600.0,nash_lateral_runoff_m/3600.0,flux_from_deep_gw_to_chan_m/3600.0,
                           Qout_m/3600.0 );
  */
  
  /*
  // WRITE OUTPUTS IN mm/h ffor aiding in interpretation
  fprintf(out_fptr,"%16.8lf %lf %lf %lf %lf %lf %lf\n",jdate_start+(double)tstep*1.0/24.0,
                           timestep_rainfall_input_m*1000.0,Schaake_output_runoff_m*1000.0,
                           giuh_runoff_m*1000.0,nash_lateral_runoff_m*1000.0, flux_from_deep_gw_to_chan_m*1000.0,
                           Qout_m*1000.0 );
  */
  // WRITE OUTPUTS IN mm/h ffor aiding in interpretation
  //fprintf(out_fptr,"%d, %ld-%02ld-%02ld %02ld:%02ld:%02ld, %16.8lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
  double total_discharge = Qout_m*catchment_area_km2*1000000.0/3600.0;
  fprintf(out_fptr,"%d,%ld-%02ld-%02ld %02ld:%02ld:%02ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                           tstep,sim_year[tstep],sim_month[tstep],sim_day[tstep],
                           sim_hour[tstep],sim_minute[tstep],sim_second[tstep],
                           //jdate_start+(double)tstep*1.0/24.0,
                           timestep_rainfall_input_m,Schaake_output_runoff_m,
                           giuh_runoff_m, nash_lateral_runoff_m, flux_from_deep_gw_to_chan_m,
                           Qout_m, total_discharge);
  }  // END OF TIME LOOP


//
// PERFORM MASS BALANCE CHECKS AND WRITE RESULTS TO stderr.
//----------------------------------------------------------

volend= soil_reservoir.storage_m+gw_reservoir.storage_m;
vol_in_gw_end = gw_reservoir.storage_m;

#ifdef DEBUG
        fclose(out_debug_fptr);
#endif

// the GIUH queue might have water in it at the end of the simulation, so sum it up.
for(i=0;i<num_giuh_ordinates;i++) vol_end_giuh+=runoff_queue_m_per_timestep[i];

for(i=0;i<num_lateral_flow_nash_reservoirs;i++)  vol_in_nash_end+=nash_storage[i];


vol_soil_end=soil_reservoir.storage_m;

double global_residual;
double schaake_residual;
double giuh_residual;
double soil_residual;
double nash_residual;
double gw_residual;

global_residual=volstart+volin-volout-volend;
fprintf(stderr,"GLOBAL MASS BALANCE\n");
fprintf(stderr,"  initial volume: %8.4lf m\n",volstart);
fprintf(stderr,"    volume input: %8.4lf m\n",volin);
fprintf(stderr,"   volume output: %8.4lf m\n",volout);
fprintf(stderr,"    final volume: %8.4lf m\n",volend);
fprintf(stderr,"        residual: %6.4e m\n",global_residual);
if(volin>0.0) fprintf(stderr,"global pct. err: %6.4e percent of inputs\n",global_residual/volin*100.0);
else          fprintf(stderr,"global pct. err: %6.4e percent of initial\n",global_residual/volstart*100.0);
if(!is_fabs_less_than_epsilon(global_residual,1.0e-12)) 
              fprintf(stderr, "WARNING: GLOBAL MASS BALANCE CHECK FAILED\n");

schaake_residual=volin-vol_sch_runoff-vol_sch_infilt;
fprintf(stderr," SCHAAKE MASS BALANCE\n");
fprintf(stderr,"  surface runoff: %8.4lf m\n",vol_sch_runoff);
fprintf(stderr,"    infiltration: %8.4lf m\n",vol_sch_infilt);
fprintf(stderr,"schaake residual: %6.4e m\n",schaake_residual);  // should equal 0.0
if(!is_fabs_less_than_epsilon(schaake_residual,1.0e-12))
              fprintf(stderr,"WARNING: SCHAAKE PARTITIONING MASS BALANCE CHECK FAILED\n");

giuh_residual=vol_out_giuh-vol_sch_runoff-vol_end_giuh;
fprintf(stderr," GIUH MASS BALANCE\n");
fprintf(stderr,"  vol. into giuh: %8.4lf m\n",vol_sch_runoff);
fprintf(stderr,"   vol. out giuh: %8.4lf m\n",vol_out_giuh);
fprintf(stderr," vol. end giuh q: %8.4lf m\n",vol_end_giuh);
fprintf(stderr,"   giuh residual: %6.4e m\n",giuh_residual);  // should equal zero
if(!is_fabs_less_than_epsilon(giuh_residual,1.0e-12))
              fprintf(stderr,"WARNING: GIUH MASS BALANCE CHECK FAILED\n");

soil_residual=vol_soil_start+vol_sch_infilt-vol_soil_to_lat_flow-vol_soil_end-vol_to_gw;
fprintf(stderr," SOIL WATER CONCEPTUAL RESERVOIR MASS BALANCE\n");
fprintf(stderr,"   init soil vol: %8.4lf m\n",vol_soil_start);     
fprintf(stderr,"  vol. into soil: %8.4lf m\n",vol_sch_infilt);
fprintf(stderr,"vol.soil2latflow: %8.4lf m\n",vol_soil_to_lat_flow);
fprintf(stderr," vol. soil to gw: %8.4lf m\n",vol_soil_to_gw);
fprintf(stderr," final vol. soil: %8.4lf m\n",vol_soil_end);   
fprintf(stderr,"vol. soil resid.: %6.4e m\n",soil_residual);
if(!is_fabs_less_than_epsilon(soil_residual,1.0e-12))
               fprintf(stderr,"WARNING: SOIL CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");

nash_residual= vol_in_nash - vol_out_nash - vol_in_nash_end;
fprintf(stderr," NASH CASCADE CONCEPTUAL RESERVOIR MASS BALANCE\n");
fprintf(stderr,"    vol. to nash: %8.4lf m\n",vol_in_nash);
fprintf(stderr,"  vol. from nash: %8.4lf m\n",vol_out_nash);
fprintf(stderr," final vol. nash: %8.4lf m\n",vol_in_nash_end);
fprintf(stderr,"nash casc resid.: %6.4e m\n",nash_residual);
if(!is_fabs_less_than_epsilon(nash_residual,1.0e-12))
               fprintf(stderr,"WARNING: NASH CASCADE CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");


gw_residual=vol_in_gw_start+vol_to_gw-vol_from_gw-vol_in_gw_end;
fprintf(stderr," GROUNDWATER CONCEPTUAL RESERVOIR MASS BALANCE\n");
fprintf(stderr,"init gw. storage: %8.4lf m\n",vol_in_gw_start);
fprintf(stderr,"       vol to gw: %8.4lf m\n",vol_to_gw);
fprintf(stderr,"     vol from gw: %8.4lf m\n",vol_from_gw);
fprintf(stderr,"final gw.storage: %8.4lf m\n",vol_in_gw_end);
fprintf(stderr,"    gw. residual: %6.4e m\n",gw_residual);
if(!is_fabs_less_than_epsilon(gw_residual,1.0e-12))
               fprintf(stderr,"WARNING: GROUNDWATER CONCEPTUAL RESERVOIR MASS BALANCE CHECK FAILED\n");

fclose(out_fptr);

}  //<<<<<<<<<<<<  END OF MAIN PROGRAM



/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS JUST UTILITY MEMORY AND TIME FUNCTION CODE */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


/*####################################################################*/
/*########################### PARSE LINE #############################*/
/*####################################################################*/
void parse_aorc_line(char *theString,long *year,long *month, long *day,long *hour,long *minute, double *second,
                struct aorc_forcing_data *aorc)
{
char str[20];
long yr,mo,da,hr,mi;
double mm,julian,se;
float val;
int i,start,end,len;
int yes_pm,wordlen;
char theWord[150];

len=strlen(theString);

start=0; /* begin at the beginning of theString */
get_word(theString,&start,&end,theWord,&wordlen);
*year=atol(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
*month=atol(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
*day=atol(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
*hour=atol(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
*minute=atol(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
*second=(double)atof(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
aorc->precip_kg_per_m2=atof(theWord);
              
get_word(theString,&start,&end,theWord,&wordlen);
aorc->incoming_longwave_W_per_m2=atof(theWord);   

get_word(theString,&start,&end,theWord,&wordlen);
aorc->incoming_shortwave_W_per_m2=atof(theWord);   

get_word(theString,&start,&end,theWord,&wordlen);
aorc->surface_pressure_Pa=atof(theWord);           

get_word(theString,&start,&end,theWord,&wordlen);
aorc->specific_humidity_2m_kg_per_kg=atof(theWord);

get_word(theString,&start,&end,theWord,&wordlen);
aorc->air_temperature_2m_K=atof(theWord);          

get_word(theString,&start,&end,theWord,&wordlen);
aorc->u_wind_speed_10m_m_per_s=atof(theWord);      

get_word(theString,&start,&end,theWord,&wordlen);
aorc->v_wind_speed_10m_m_per_s=atof(theWord);      

get_word(theString,&start,&end,theWord,&wordlen);
aorc->precip_kg_per_m2=atof(theWord);

  
return;
}

/*####################################################################*/
/*############################## GET WORD ############################*/
/*####################################################################*/
void get_word(char *theString,int *start,int *end,char *theWord,int *wordlen)
{
int i,lenny,j;
lenny=strlen(theString);

while(theString[*start]==' ' || theString[*start]=='\t')
  {
  (*start)++;
  };
  
j=0;
for(i=(*start);i<lenny;i++)
  {
  if(theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' && theString[i]!='-')
    {
      theWord[j]=theString[i];
      j++;
    }
  else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
       && (theString[i]=='-' && theString[i-1]==',') )
    {
      theWord[j]=theString[i];
      j++;
    }
  else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
       && (theString[i]=='-' && theString[i-1]=='e') )
    {
      theWord[j]=theString[i];
      j++;
    }
  else
    {
    break;
    }
  }
theWord[j]='\0';
*start=i+1;
*wordlen=strlen(theWord);
return;
}

/****************************************/
void itwo_alloc(int ***array,int rows, int cols)
{
int  i,frows,fcols, numgood=0;
int error=0;

if ((rows==0)||(cols==0))
  {
  printf("Error: Attempting to allocate array of size 0\n");
  exit;
  }

frows=rows+1;  /* added one for FORTRAN numbering */
fcols=cols+1;  /* added one for FORTRAN numbering */

*array=(int **)malloc(frows*sizeof(int *));
if (*array) 
  {
  memset((*array), 0, frows*sizeof(int*));
  for (i=0; i<frows; i++)
    {
    (*array)[i] =(int *)malloc(fcols*sizeof(int ));
    if ((*array)[i] == NULL)
      {
      error = 1;
      numgood = i;
      i = frows;
      }
     else memset((*array)[i], 0, fcols*sizeof(int )); 
     }
   }
return;
}



void dtwo_alloc(double ***array,int rows, int cols)
{
int  i,frows,fcols, numgood=0;
int error=0;

if ((rows==0)||(cols==0))
  {
  printf("Error: Attempting to allocate array of size 0\n");
  exit;
  }

frows=rows+1;  /* added one for FORTRAN numbering */
fcols=cols+1;  /* added one for FORTRAN numbering */

*array=(double **)malloc(frows*sizeof(double *));
if (*array) 
  {
  memset((*array), 0, frows*sizeof(double *));
  for (i=0; i<frows; i++)
    {
    (*array)[i] =(double *)malloc(fcols*sizeof(double ));
    if ((*array)[i] == NULL)
      {
      error = 1;
      numgood = i;
      i = frows;
      }
     else memset((*array)[i], 0, fcols*sizeof(double )); 
     }
   }
return;
}



void d_alloc(double **var,int size)
{
  size++;  /* just for safety */

   *var = (double *)malloc(size * sizeof(double));
   if (*var == NULL)
      {
      printf("Problem allocating memory for array in d_alloc.\n");
      return;
      }
   else memset(*var,0,size*sizeof(double));
   return;
}

void i_alloc(int **var,int size)
{
   size++;  /* just for safety */

   *var = (int *)malloc(size * sizeof(int));
   if (*var == NULL)
      {
      printf("Problem allocating memory in i_alloc\n");
      return; 
      }
   else memset(*var,0,size*sizeof(int));
   return;
}

void l_alloc(long **var,int size)
{
   size++;  /* just for safety */

   *var = (long *)malloc(size * sizeof(long));
   if (*var == NULL)
      {
      printf("Problem allocating memory in l_alloc\n");
      return;
      }
   else memset(*var,0,size*sizeof(long));
   return;
}

/*
 * convert Gregorian days to Julian date
 *
 * Modify as needed for your application.
 *
 * The Julian day starts at noon of the Gregorian day and extends
 * to noon the next Gregorian day.
 *
 */
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
*/


double greg_2_jul(
long year, 
long mon, 
long day, 
long h, 
long mi, 
double se)
{
    long m = mon, d = day, y = year;
    long c, ya, j;
    double seconds = h * 3600.0 + mi * 60 + se;

    if (m > 2)
	m -= 3;
    else {
	m += 9;
	--y;
    }
    c = y / 100L;
    ya = y - (100L * c);
    j = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
    if (seconds < 12 * 3600.0) {
	j--;
	seconds += 12.0 * 3600.0;
    }
    else {
	seconds = seconds - 12.0 * 3600.0;
    }
    return (j + (seconds / 3600.0) / 24.0);
}

/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
**
** Modified by FLO 4/99 to account for nagging round off error 
**
*/
void calc_date(double jd, long *y, long *m, long *d, long *h, long *mi,
               double *sec)
{
    static int ret[4];

    long j;
    double tmp; 
    double frac;

    j=(long)jd;
    frac = jd - j;

    if (frac >= 0.5) {
	frac = frac - 0.5;
        j++;
    }
    else {
	frac = frac + 0.5;
    }

    ret[3] = (j + 1L) % 7L;
    j -= 1721119L;
    *y = (4L * j - 1L) / 146097L;
    j = 4L * j - 1L - 146097L * *y;
    *d = j / 4L;
    j = (4L * *d + 3L) / 1461L;
    *d = 4L * *d + 3L - 1461L * j;
    *d = (*d + 4L) / 4L;
    *m = (5L * *d - 3L) / 153L;
    *d = 5L * *d - 3 - 153L * *m;
    *d = (*d + 5L) / 5L;
    *y = 100L * *y + j;
    if (*m < 10)
	*m += 3;
    else {
	*m -= 9;
	*y=*y+1; /* Invalid use: *y++. Modified by Tony */
    }

    /* if (*m < 3) *y++; */
    /* incorrectly repeated the above if-else statement. Deleted by Tony.*/

    tmp = 3600.0 * (frac * 24.0);
    *h = (long) (tmp / 3600.0);
    tmp = tmp - *h * 3600.0;
    *mi = (long) (tmp / 60.0);
    *sec = tmp - *mi * 60.0;
}

int dayofweek(double j)
{
    j += 0.5;
    return (int) (j + 1) % 7;
}

