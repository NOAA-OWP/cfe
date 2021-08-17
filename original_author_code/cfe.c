#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
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

struct conceptual_reservoir
{
// this data structure describes a nonlinear reservoir having two outlets, one primary with an activation
// threshold that may be zero, and a secondary outlet with a threshold that may be zero
// this will also simulate a linear reservoir by setting the exponent parameter to 1.0 iff is_exponential==FALSE
// iff is_exponential==TRUE, then it uses the exponential discharge function from the NWM V2.0 forumulation
// as the primary discharge with a zero threshold, and does not calculate a secondary discharge.
//--------------------------------------------------------------------------------------------------
int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
double storage_max_m;   // maximum storage in this reservoir
double storage_m;       // state variable.
double coeff_primary;    // the primary outlet
double exponent_primary;
double storage_threshold_primary_m;
double storage_threshold_secondary_m;
double coeff_secondary;
double exponent_secondary;
};

struct NWM_soil_parameters
{
// using same variable names as used in NWM.  <sorry>
double smcmax;  // effective porosity [V/V]
double wltsmc;  // wilting point soil moisture content [V/V]
double satdk;   // saturated hydraulic conductivity [m s-1]
double satpsi;	// saturated capillary head [m]
double bb;      // beta exponent on Clapp-Hornberger (1978) soil water relations [-]
double mult;    // the multiplier applied to satdk to route water rapidly downslope
double slop;   // this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow.
double D;       // soil depth [m]
};



// function prototypes
// --------------------------------
extern void Schaake_partitioning_scheme(double dt, double magic_number, double deficit, double qinsur,
                                        double *runsrf, double *pddum);

extern void conceptual_reservoir_flux_calc(struct conceptual_reservoir *da_reservoir,
                                           double *primary_flux_m,double *secondary_flux_m);

extern double convolution_integral(double runoff_m, int num_giuh_ordinates, 
                                   double *giuh_ordinates, double *runoff_queue_m_per_timestep);
                                   
extern double nash_cascade(double flux_lat_m,int num_lateral_flow_nash_reservoirs,
                           double K_nash,double *nash_storage_arr);
extern int is_fabs_less_than_epsilon(double a,double epsilon);

extern void cfe(
        double *soil_reservoir_storage_deficit_m_ptr,
        struct NWM_soil_parameters NWM_soil_params_struct,
        struct conceptual_reservoir *soil_reservoir_struct,
        double timestep_h,
        double Schaake_adjusted_magic_constant_by_soil_type,
        double timestep_rainfall_input_m,
        double *Schaake_output_runoff_m_ptr,
        double *infiltration_depth_m_ptr,
        double *flux_overland_m_ptr,
        double *vol_sch_runoff_ptr,
        double *vol_sch_infilt_ptr,
        double *flux_perc_m_ptr,
        double *vol_to_soil_ptr,
        double *flux_lat_m_ptr,
        double *gw_reservoir_storage_deficit_m_ptr,
        struct conceptual_reservoir *gw_reservoir_struct,
        double *vol_to_gw_ptr,
        double *vol_soil_to_gw_ptr,
        double *vol_soil_to_lat_flow_ptr,
        double *volout_ptr,
        double *flux_from_deep_gw_to_chan_m_ptr,
        double *vol_from_gw_ptr,
        double *giuh_runoff_m_ptr,
        int num_giuh_ordinates,
        double *giuh_ordinates_arr,
        double *runoff_queue_m_per_timestep_arr,
        double *vol_out_giuh_ptr,
        double *nash_lateral_runoff_m_ptr,
        int num_lateral_flow_nash_reservoirs,
        double K_nash,
        double *nash_storage_arr,
        double *vol_in_nash_ptr,
        double *vol_out_nash_ptr,
        double *Qout_m_ptr
    );

  
//####################################################################################################
//####################################################################################################
//####################################################################################################
//####################################################################################################
// CFE STATE SPACE FUNCTION // #######################################################################
extern void cfe(
        double *soil_reservoir_storage_deficit_m_ptr,
        struct NWM_soil_parameters NWM_soil_params_struct,
        struct conceptual_reservoir *soil_reservoir_struct,
        double timestep_h,
        double Schaake_adjusted_magic_constant_by_soil_type,
        double timestep_rainfall_input_m,
        double *Schaake_output_runoff_m_ptr,
        double *infiltration_depth_m_ptr,
        double *flux_overland_m_ptr,
        double *vol_sch_runoff_ptr,
        double *vol_sch_infilt_ptr,
        double *flux_perc_m_ptr,
        double *vol_to_soil_ptr,
        double *flux_lat_m_ptr,
        double *gw_reservoir_storage_deficit_m_ptr,
        struct conceptual_reservoir *gw_reservoir_struct,
        double *vol_to_gw_ptr,
        double *vol_soil_to_gw_ptr,
        double *vol_soil_to_lat_flow_ptr,
        double *volout_ptr,
        double *flux_from_deep_gw_to_chan_m_ptr,
        double *vol_from_gw_ptr,
        double *giuh_runoff_m_ptr,
        int num_giuh_ordinates,
        double *giuh_ordinates_arr,
        double *runoff_queue_m_per_timestep_arr,
        double *vol_out_giuh_ptr,
        double *nash_lateral_runoff_m_ptr,
        int num_lateral_flow_nash_reservoirs,
        double K_nash,
        double *nash_storage_arr,
        double *vol_in_nash_ptr,
        double *vol_out_nash_ptr,
        double *Qout_m_ptr
    ){                      // #######################################################################
// CFE STATE SPACE FUNCTION // #######################################################################

    // ####    COPY THE MODEL FUNCTION STATE SPACE    ####
    double soil_reservoir_storage_deficit_m = *soil_reservoir_storage_deficit_m_ptr;
    double Schaake_output_runoff_m          = *Schaake_output_runoff_m_ptr;
    double infiltration_depth_m             = *infiltration_depth_m_ptr;
    double flux_overland_m                  = *flux_overland_m_ptr;
    double vol_sch_runoff                   = *vol_sch_runoff_ptr;
    double vol_sch_infilt                   = *vol_sch_infilt_ptr;
    double flux_perc_m                      = *flux_perc_m_ptr;
    double vol_to_soil                      = *vol_to_soil_ptr;
    double flux_lat_m                       = *flux_lat_m_ptr;
    double gw_reservoir_storage_deficit_m   = *gw_reservoir_storage_deficit_m_ptr;
    double vol_to_gw                        = *vol_to_gw_ptr;
    double vol_soil_to_gw                   = *vol_soil_to_gw_ptr;
    double vol_soil_to_lat_flow             = *vol_soil_to_lat_flow_ptr;
    double volout                           = *volout_ptr;
    double flux_from_deep_gw_to_chan_m      = *flux_from_deep_gw_to_chan_m_ptr;
    double vol_from_gw                      = *vol_from_gw_ptr;
    double giuh_runoff_m                    = *giuh_runoff_m_ptr;
    double vol_out_giuh                     = *vol_out_giuh_ptr;
    double nash_lateral_runoff_m            = *nash_lateral_runoff_m_ptr;
    double vol_in_nash                      = *vol_in_nash_ptr;
    double vol_out_nash                     = *vol_out_nash_ptr;
    double Qout_m                           = *Qout_m_ptr;

    // LOCAL PARAMETERS
    double diff=0;
    double primary_flux=0;
    double secondary_flux=0;
    double lateral_flux=0;      // flux from soil to lateral flow Nash cascade +to cascade
    double percolation_flux=0;  // flux from soil to gw nonlinear researvoir, +downward

  //##################################################
  // partition rainfall using Schaake function
  //##################################################

  soil_reservoir_storage_deficit_m=(NWM_soil_params_struct.smcmax*NWM_soil_params_struct.D-soil_reservoir_struct->storage_m);
  
  Schaake_partitioning_scheme(timestep_h,Schaake_adjusted_magic_constant_by_soil_type,soil_reservoir_storage_deficit_m,
                              timestep_rainfall_input_m,&Schaake_output_runoff_m,&infiltration_depth_m);
  
  // check to make sure that there is storage available in soil to hold the water that does not runoff
  //--------------------------------------------------------------------------------------------------
  if(soil_reservoir_storage_deficit_m<infiltration_depth_m)
    {
    Schaake_output_runoff_m+=(infiltration_depth_m-soil_reservoir_storage_deficit_m);  // put won't fit back into runoff
    infiltration_depth_m=soil_reservoir_storage_deficit_m;
    soil_reservoir_struct->storage_m=soil_reservoir_struct->storage_max_m;
    }
#ifdef DEBUG
  printf("After Schaake function: rain:%8.5lf mm  runoff:%8.5lf mm  infiltration:%8.5lf mm  residual:%e m\n",
                                 timestep_rainfall_input_m,Schaake_output_runoff_m*1000.0,infiltration_depth_m*1000.0,
                                 timestep_rainfall_input_m-Schaake_output_runoff_m-infiltration_depth_m);
#endif

  flux_overland_m=Schaake_output_runoff_m;

  vol_sch_runoff   += flux_overland_m;
  vol_sch_infilt   += infiltration_depth_m;
  
  // put infiltration flux into soil conceptual reservoir.  If not enough room
  // limit amount transferred to deficit
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  REDUNDANT soil_reservoir_storage_deficit_m=soil_reservoir.storage_max_m-soil_reservoir.storage_m;  <- commented out by FLO based on comments from Bartel
  
  if(flux_perc_m>soil_reservoir_storage_deficit_m)
    {
    diff=flux_perc_m-soil_reservoir_storage_deficit_m;  // the amount that there is not capacity ffor
    infiltration_depth_m=soil_reservoir_storage_deficit_m;  
    vol_sch_runoff+=diff;  // send excess water back to GIUH runoff
    vol_sch_infilt-=diff;  // correct overprediction of infilt.
    flux_overland_m+=diff; // bug found by Nels.  This was missing and fixes it.
    }

  vol_to_soil              += infiltration_depth_m; 
  soil_reservoir_struct->storage_m += infiltration_depth_m;  // put the infiltrated water in the soil.

  
  // calculate fluxes from the soil storage into the deep groundwater (percolation) and to lateral subsurface flow
  //--------------------------------------------------------------------------------------------------------------
  conceptual_reservoir_flux_calc(soil_reservoir_struct,&percolation_flux,&lateral_flux);
  flux_perc_m=percolation_flux;  // m/h   <<<<<<<<<<<  flux of percolation from soil to g.w. reservoir >>>>>>>>>
  
  flux_lat_m=lateral_flux;  // m/h        <<<<<<<<<<<  flux into the lateral flow Nash cascade >>>>>>>>
  

  // calculate flux of base flow from deep groundwater reservoir to channel
  //--------------------------------------------------------------------------
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  gw_reservoir_storage_deficit_m= gw_reservoir_struct->storage_max_m-gw_reservoir_struct->storage_m;
  
  // limit amount transferred to deficit iff there is insuffienct avail. storage
  if(flux_perc_m>gw_reservoir_storage_deficit_m)
    {
    diff=flux_perc_m-gw_reservoir_storage_deficit_m;
    flux_perc_m=gw_reservoir_storage_deficit_m;
    vol_sch_runoff+=diff;  // send excess water back to GIUH runoff
    vol_sch_infilt-=diff;  // correct overprediction of infilt.
    }
    
  vol_to_gw                +=flux_perc_m;
  vol_soil_to_gw           +=flux_perc_m;
  
  gw_reservoir_struct->storage_m   += flux_perc_m;
  soil_reservoir_struct->storage_m -= flux_perc_m;
  soil_reservoir_struct->storage_m -= flux_lat_m;
  vol_soil_to_lat_flow     += flux_lat_m;  //TODO add this to nash cascade as input
  volout=volout+flux_lat_m;
  
  conceptual_reservoir_flux_calc(gw_reservoir_struct,&primary_flux,&secondary_flux);
  
  flux_from_deep_gw_to_chan_m=primary_flux;  // m/h   <<<<<<<<<< BASE FLOW FLUX >>>>>>>>>
  vol_from_gw+=flux_from_deep_gw_to_chan_m;
  
  // in the instance of calling the gw reservoir the secondary flux should be zero- verify
  if(is_fabs_less_than_epsilon(secondary_flux,1.0e-09)==FALSE) printf("problem with nonzero flux point 1\n");

  
  // adjust state of deep groundwater conceptual nonlinear reservoir
  //-----------------------------------------------------------------
  
  gw_reservoir_struct->storage_m -= flux_from_deep_gw_to_chan_m;

  
  // Solve the convolution integral ffor this time step 

  giuh_runoff_m = convolution_integral(Schaake_output_runoff_m,num_giuh_ordinates,
                                              giuh_ordinates_arr,runoff_queue_m_per_timestep_arr);
  vol_out_giuh+=giuh_runoff_m;

  volout+=giuh_runoff_m;
  volout+=flux_from_deep_gw_to_chan_m;
  
  // Route lateral flow through the Nash cascade.
  nash_lateral_runoff_m = nash_cascade(flux_lat_m,num_lateral_flow_nash_reservoirs,
                                       K_nash,nash_storage_arr);
  vol_in_nash   += flux_lat_m;
  vol_out_nash  += nash_lateral_runoff_m;

#ifdef DEBUG
        fprintf(out_debug_fptr,"%d %lf %lf\n",tstep,flux_lat_m,nash_lateral_runoff_m);
#endif

  Qout_m = giuh_runoff_m + nash_lateral_runoff_m + flux_from_deep_gw_to_chan_m;
    
    // #### COPY BACK FUNCTION VALUES    ####    
    *soil_reservoir_storage_deficit_m_ptr = soil_reservoir_storage_deficit_m;
    *Schaake_output_runoff_m_ptr          = Schaake_output_runoff_m;
    *infiltration_depth_m_ptr             = infiltration_depth_m;
    *flux_overland_m_ptr                  = flux_overland_m;
    *vol_sch_runoff_ptr                   = vol_sch_runoff;
    *vol_sch_infilt_ptr                   = vol_sch_infilt;
    *flux_perc_m_ptr                      = flux_perc_m;
    *vol_to_soil_ptr                      = vol_to_soil;
    *flux_lat_m_ptr                       = flux_lat_m;
    *gw_reservoir_storage_deficit_m_ptr   = gw_reservoir_storage_deficit_m;
    *vol_to_gw_ptr                        = vol_to_gw;
    *vol_soil_to_gw_ptr                   = vol_soil_to_gw;
    *vol_soil_to_lat_flow_ptr             = vol_soil_to_lat_flow;
    *volout_ptr                           = volout;
    *flux_from_deep_gw_to_chan_m_ptr      = flux_from_deep_gw_to_chan_m;
    *vol_from_gw_ptr                      = vol_from_gw;
    *giuh_runoff_m_ptr                    = giuh_runoff_m;
    *vol_out_giuh_ptr                     = vol_out_giuh;
    *nash_lateral_runoff_m_ptr            = nash_lateral_runoff_m;
    *vol_in_nash_ptr                      = vol_in_nash;
    *vol_out_nash_ptr                     = vol_out_nash;
    *Qout_m_ptr                           = Qout_m;



} // END CFE STATE SPACE FUNCTIONS
  //####################################################################################################
  //####################################################################################################



//##############################################################
//###################   NASH CASCADE   #########################
//##############################################################
extern double nash_cascade(double flux_lat_m,int num_lateral_flow_nash_reservoirs,
                           double K_nash,double *nash_storage_arr)
{
//##############################################################
// Solve for the flow through the Nash cascade to delay the 
// arrival of the lateral flow into the channel
//##############################################################
// local vars
int i;
double outflow_m;
static double Q[MAX_NUM_NASH_CASCADE];

//Loop through reservoirs
for(i = 0; i < num_lateral_flow_nash_reservoirs; i++)
  {
  Q[i] = K_nash*nash_storage_arr[i];
  nash_storage_arr[i]  -= Q[i];

  if (i==0) nash_storage_arr[i] += flux_lat_m; 
  else      nash_storage_arr[i] +=  Q[i-1];
  }

/*  Get Qout */
outflow_m = Q[num_lateral_flow_nash_reservoirs-1];

//Return the flow output
return (outflow_m);

}

//##############################################################
//############### GIUH CONVOLUTION INTEGRAL   ##################
//##############################################################
extern double convolution_integral(double runoff_m,int num_giuh_ordinates, 
                                   double *giuh_ordinates, double *runoff_queue_m_per_timestep)
{
//##############################################################
// This function solves the convolution integral involving N
//  GIUH ordinates.
//##############################################################
double runoff_m_now;
int N,i;

N=num_giuh_ordinates;
runoff_queue_m_per_timestep[N]=0.0;

for(i=0;i<N;i++)
  {
  runoff_queue_m_per_timestep[i]+=giuh_ordinates[i]*runoff_m;
  }
runoff_m_now=runoff_queue_m_per_timestep[0];

for(i=1;i<N;i++)  // shift all the entries in preperation ffor the next timestep
  {
  runoff_queue_m_per_timestep[i-1]=runoff_queue_m_per_timestep[i];
  }
runoff_queue_m_per_timestep[N-1]=0.0;

return(runoff_m_now);
}

//##############################################################
//########## SINGLE OUTLET EXPONENTIAL RESERVOIR ###############
//##########                -or-                 ###############
//##########    TWO OUTLET NONLINEAR RESERVOIR   ###############
//################################################################
// This function calculates the flux from a linear, or nonlinear 
// conceptual reservoir with one or two outlets, or from an
// exponential nonlinear conceptual reservoir with only one outlet.
// In the non-exponential instance, each outlet can have its own
// activation storage threshold.  Flow from the second outlet is 
// turned off by setting the discharge coeff. to 0.0.
//################################################################
extern void conceptual_reservoir_flux_calc(struct conceptual_reservoir *da_reservoir,
                                           double *primary_flux_m,double *secondary_flux_m)
{
//struct conceptual_reservoir  <<<<INCLUDED HERE FOR REFERENCE.>>>>
//{
// int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
// double storage_max_m;
// double storage_m;
// double coeff_primary;
// double exponent_secondary;
// double storage_threshold_primary_m;
// double storage_threshold_secondary_m;
// double coeff_secondary;
// double exponent_secondary;
// };
// THIS FUNCTION CALCULATES THE FLUXES FROM A CONCEPTUAL NON-LINEAR (OR LINEAR) RESERVOIR WITH TWO OUTLETS
// all fluxes calculated by this routine are instantaneous with units of the coefficient.

//local variables
double storage_above_threshold_m;

if(da_reservoir->is_exponential==TRUE)  // single outlet reservoir like the NWM V1.2 exponential conceptual gw reservoir
  {
  // calculate the one flux and return.
  *primary_flux_m=da_reservoir->coeff_primary*
                    (exp(da_reservoir->exponent_primary*da_reservoir->storage_m/da_reservoir->storage_max_m)-1.0);
  *secondary_flux_m=0.0;
  return;
  }
// code goes past here iff it is not a single outlet exponential deep groundwater reservoir of the NWM variety
// The vertical outlet is assumed to be primary and satisfied first.

*primary_flux_m=0.0;
storage_above_threshold_m=da_reservoir->storage_m-da_reservoir->storage_threshold_primary_m;
if(storage_above_threshold_m>0.0)
  {
  // flow is possible from the primary outlet
  *primary_flux_m=da_reservoir->coeff_primary*
                pow(storage_above_threshold_m/(da_reservoir->storage_max_m-da_reservoir->storage_threshold_primary_m),
                    da_reservoir->exponent_primary);
  if(*primary_flux_m > storage_above_threshold_m) 
                    *primary_flux_m=storage_above_threshold_m;  // limit to max. available
  }
*secondary_flux_m=0.0;
storage_above_threshold_m=da_reservoir->storage_m-da_reservoir->storage_threshold_secondary_m;
if(storage_above_threshold_m>0.0)
  {
  // flow is possible from the secondary outlet
  *secondary_flux_m=da_reservoir->coeff_secondary*
                  pow(storage_above_threshold_m/(da_reservoir->storage_max_m-da_reservoir->storage_threshold_secondary_m),
                      da_reservoir->exponent_secondary);
  if(*secondary_flux_m > (storage_above_threshold_m-(*primary_flux_m))) 
                    *secondary_flux_m=storage_above_threshold_m-(*primary_flux_m);  // limit to max. available
  }
return;
}


//##############################################################
//#########   SCHAAKE RUNOFF PARTITIONING SCHEME   #############
//##############################################################
void Schaake_partitioning_scheme(double timestep_h, double Schaake_adjusted_magic_constant_by_soil_type, 
           double column_total_soil_moisture_deficit_m,
           double water_input_depth_m,double *surface_runoff_depth_m,double *infiltration_depth_m)
{


/*! ===============================================================================
  This subtroutine takes water_input_depth_m and partitions it into surface_runoff_depth_m and
  infiltration_depth_m using the scheme from Schaake et al. 1996. 
! --------------------------------------------------------------------------------
! ! modified by FLO April 2020 to eliminate reference to ice processes, 
! ! and to de-obfuscate and use descriptive and dimensionally consistent variable names.
! --------------------------------------------------------------------------------
    IMPLICIT NONE
! --------------------------------------------------------------------------------
! inputs
  double timestep_h
  double Schaake_adjusted_magic_constant_by_soil_type = C*Ks(soiltype)/Ks_ref, where C=3, and Ks_ref=2.0E-06 m/s
  double column_total_soil_moisture_deficit_m
  double water_input_depth_m  amount of water input to soil surface this time step [m]

! outputs
  double surface_runoff_depth_m      amount of water partitioned to surface water this time step [m]


--------------------------------------------------------------------------------*/
int k;
double timestep_d,Schaake_parenthetical_term,Ic,Px,infilt_dep_m;


if(0.0 < water_input_depth_m) 
  {
  if (0.0 > column_total_soil_moisture_deficit_m)
    {
    *surface_runoff_depth_m=water_input_depth_m;
    *infiltration_depth_m=0.0;
    }
  else
    {
    // partition time-step total applied water as per Schaake et al. 1996.
                                         // change from dt in [s] to dt1 in [d] because kdt has units of [d^(-1)] 
    timestep_d = timestep_h /24.0;    // timestep_d is the time step in days.
      
    // calculate the parenthetical part of Eqn. 34 from Schaake et al. Note the magic constant has units of [d^(-1)]
    
    Schaake_parenthetical_term = (1.0 - exp ( - Schaake_adjusted_magic_constant_by_soil_type * timestep_d));      
    
    // From Schaake et al. Eqn. 2., using the column total moisture deficit 
    // BUT the way it is used here, it is the cumulative soil moisture deficit in the entire soil profile. 
    // "Layer" info not used in this subroutine in noah-mp, except to sum up the total soil moisture storage.
    // NOTE: when column_total_soil_moisture_deficit_m becomes zero, which occurs when the soil column is saturated, 
    // then Ic=0, where Ic in the Schaake paper is called the "spatially averaged infiltration capacity", 
    // and is defined in Eqn. 12. 

    Ic = column_total_soil_moisture_deficit_m * Schaake_parenthetical_term; 
                                     
    Px=water_input_depth_m;   // Total water input to partitioning scheme this time step [m]
  
    // This is eqn 24 from Schaake et al.  NOTE: this is 0 in the case of a saturated soil column, when Ic=0.  
    // Physically happens only if soil has no-flow lower b.c.
    
    *infiltration_depth_m = (Px * (Ic / (Px + Ic)));  


    if( 0.0 < (water_input_depth_m-(*infiltration_depth_m)) )
      {
      *surface_runoff_depth_m = water_input_depth_m - (*infiltration_depth_m);
      }
    else  *surface_runoff_depth_m=0.0;
    *infiltration_depth_m =  water_input_depth_m - (*surface_runoff_depth_m);
    }
  }
else
  {
  *surface_runoff_depth_m = 0.0;
  *infiltration_depth_m = 0.0;
  }
return;
}

extern int is_fabs_less_than_epsilon(double a,double epsilon)  // returns true if fabs(a)<epsilon
{
if(fabs(a)<epsilon) return(TRUE);
else                return(FALSE);
}

