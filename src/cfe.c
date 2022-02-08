#include "../include/cfe.h"
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a < _b ? _a : _b; })

// CFE STATE SPACE FUNCTION // #######################################################################
// Adapted version of Conceptual Functional Equivalent model re-written in state-space form July, 2021
//####################################################################################################
extern void cfe(
        double *soil_reservoir_storage_deficit_m_ptr,
        struct NWM_soil_parameters NWM_soil_params_struct,
        struct conceptual_reservoir *soil_reservoir_struct,
        double timestep_h,

        /* xinanjiang_dev: since we are doing the option for Schaake and XinJiang, 
                           instead of passing in the constants
                           pass in a structure with the constants for both subroutines.
        //double Schaake_adjusted_magic_constant_by_soil_type,*/
        struct direct_runoff_parameters_structure direct_runoff_params_struct,

        double timestep_rainfall_input_m,

        /* xinanjiang_dev: rename to the general "direct runoff"
        double *Schaake_output_runoff_m_ptr,*/
        double *flux_output_direct_runoff_m,

        double *infiltration_depth_m_ptr,
        double *flux_perc_m_ptr,
        double *flux_lat_m_ptr,
        double *gw_reservoir_storage_deficit_m_ptr,
        struct conceptual_reservoir *gw_reservoir_struct,
        double *flux_from_deep_gw_to_chan_m_ptr,
        double *giuh_runoff_m_ptr,
        int num_giuh_ordinates,
        double *giuh_ordinates_arr,
        double *runoff_queue_m_per_timestep_arr,
        double *nash_lateral_runoff_m_ptr,
        int num_lateral_flow_nash_reservoirs,
        double K_nash,
        double *nash_storage_arr,
        struct evapotranspiration_structure *evap_struct,
        double *Qout_m_ptr,
        struct massbal *massbal_struct,
        double time_step_size
    ){                      // #######################################################################
// CFE STATE SPACE FUNCTION // #######################################################################

// ####    COPY THE MODEL FUNCTION STATE SPACE TO LOCAL VARIABLES    ####
// ####    Reason: so we don't have to re-write domain science code to de-reference a whole bunch of pointers
// ####        Note: all of thes variables are storages in [m] or fluxes in [m/timestep]    
    double soil_reservoir_storage_deficit_m = *soil_reservoir_storage_deficit_m_ptr;   // storage [m]

    /* xinanjiang_dev: rename to the general "direct runoff"
    double Schaake_output_runoff_m          = *Schaake_output_runoff_m_ptr;            // Schaake partitioned runoff this timestep [m]*/
    double direct_output_runoff_m          = *flux_output_direct_runoff_m;            // Schaake partitioned runoff this timestep [m]*/

    double infiltration_depth_m             = *infiltration_depth_m_ptr;               // Schaake partitioned infiltration this timestep [m]
    double flux_perc_m                      = *flux_perc_m_ptr;                        // water moved from soil reservoir to gw reservoir this timestep [m]
    double flux_lat_m                       = *flux_lat_m_ptr;                         // water moved from soil reservoir to lateral flow Nash cascad this timestep [m]
    double gw_reservoir_storage_deficit_m   = *gw_reservoir_storage_deficit_m_ptr;     // deficit in gw reservoir storage [m]
    double flux_from_deep_gw_to_chan_m      = *flux_from_deep_gw_to_chan_m_ptr;        // water moved from gw reservoir to catchment outlet nexus this timestep [m]
    double giuh_runoff_m                    = *giuh_runoff_m_ptr;                      // water leaving GIUH to outlet this timestep [m]
    double nash_lateral_runoff_m            = *nash_lateral_runoff_m_ptr;              // water leaving lateral subsurface flow Nash cascade this timestep [m]
    double Qout_m                           = *Qout_m_ptr;                             // the total runoff this timestep (GIUH+Nash+GW) [m]

    // LOCAL VARIABLES, the values of which are not important to describe the model state.  They are like notes on scrap paper.
 
    double diff=0.0;
    double primary_flux=0.0;      // pointers to these variables passed to conceptual nonlinear reservoir which has two outlets, primary & secondary
    double secondary_flux=0.0;    // pointers to these variables passed to conceptual nonlinear reservoir which has two outlets, primary & secondary
    double lateral_flux=0.0;      // flux from soil to lateral flow Nash cascade +to cascade  [m/timestep]
    double percolation_flux=0.0;  // flux from soil to gw nonlinear researvoir, +downward  [m/timestep]

  //##################################################
  // partition rainfall using Schaake function
  //##################################################

  soil_reservoir_storage_deficit_m=(NWM_soil_params_struct.smcmax*NWM_soil_params_struct.D-soil_reservoir_struct->storage_m);

  evap_struct->potential_et_m_per_timestep = evap_struct->potential_et_m_per_s * time_step_size;

  et_from_rainfall(&timestep_rainfall_input_m,evap_struct);
  massbal_struct->vol_et_from_rain = massbal_struct->vol_et_from_rain + evap_struct->actual_et_m_per_timestep;
  massbal_struct->vol_et_to_atm = massbal_struct->vol_et_to_atm + evap_struct->actual_et_m_per_timestep;
  massbal_struct->volout=massbal_struct->volout+evap_struct->actual_et_m_per_timestep;
  // NEW FLO
  if(0.0 < timestep_rainfall_input_m) 
    {
    if (direct_runoff_params_struct.surface_partitioning_scheme == Schaake)
      {
      Schaake_partitioning_scheme(timestep_h,direct_runoff_params_struct.Schaake_adjusted_magic_constant_by_soil_type,soil_reservoir_storage_deficit_m,
                                  timestep_rainfall_input_m,&direct_output_runoff_m,&infiltration_depth_m);
      }
    else if (direct_runoff_params_struct.surface_partitioning_scheme == Xinanjiang)
      {
      Xinanjiang_partitioning_scheme(timestep_rainfall_input_m, soil_reservoir_struct->storage_threshold_primary_m,
                                     soil_reservoir_struct->storage_max_m, soil_reservoir_struct->storage_m,
                                     &direct_runoff_params_struct, 
                                     &direct_output_runoff_m, &infiltration_depth_m);
      }
    else
      {
      fprintf(stderr,"Problem, must specify one of Schaake of Xinanjiang partitioning scheme.\n");
      fprintf(stderr,"Program terminating.\n");
      exit(-1);   // note -1 is arbitrary   #############BOMB################ NEW FLO
      }
    }
  else  // NEW FLO no need to call Schaake or Xinanjiang if it's not raining.
    {
    direct_output_runoff_m = 0.0;
    infiltration_depth_m = 0.0;
    }
  et_from_soil(soil_reservoir_struct, evap_struct, &NWM_soil_params_struct);
  massbal_struct->vol_et_from_soil = massbal_struct->vol_et_from_soil + evap_struct->actual_et_m_per_timestep;
  massbal_struct->vol_et_to_atm = massbal_struct->vol_et_to_atm + evap_struct->actual_et_m_per_timestep;
  massbal_struct->volout=massbal_struct->volout+evap_struct->actual_et_m_per_timestep;
  // check to make sure that there is storage available in soil to hold the water that does not runoff
  //--------------------------------------------------------------------------------------------------
  if(soil_reservoir_storage_deficit_m<infiltration_depth_m)
    {
    direct_output_runoff_m+=(infiltration_depth_m-soil_reservoir_storage_deficit_m);  // put infiltration that won't fit back into runoff
    infiltration_depth_m=soil_reservoir_storage_deficit_m;
    soil_reservoir_struct->storage_m=soil_reservoir_struct->storage_max_m;
    }
#ifdef DEBUG
  /* xinanjiang_dev
  printf("After Schaake function: rain:%8.5lf mm  runoff:%8.5lf mm  infiltration:%8.5lf mm  residual:%e m\n",
                                 timestep_rainfall_input_m,Schaake_output_runoff_m*1000.0,infiltration_depth_m*1000.0,
                                 timestep_rainfall_input_m-Schaake_output_runoff_m-infiltration_depth_m);                   */
  printf("After direct runoff function: rain:%8.5lf mm  runoff:%8.5lf mm  infiltration:%8.5lf mm  residual:%e m\n",
                                 timestep_rainfall_input_m,direct_output_runoff_m*1000.0,infiltration_depth_m*1000.0,
                                 timestep_rainfall_input_m-direct_output_runoff_m-infiltration_depth_m);
#endif  

  massbal_struct->vol_runoff   += direct_output_runoff_m;       // edit FLO
  massbal_struct->vol_infilt   += infiltration_depth_m;  // edit FLO
  
  // put infiltration flux into soil conceptual reservoir.  If not enough room
  // limit amount transferred to deficit
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if(flux_perc_m>soil_reservoir_storage_deficit_m)
    {
    diff=flux_perc_m-soil_reservoir_storage_deficit_m;  // the amount that there is not capacity ffor
    infiltration_depth_m=soil_reservoir_storage_deficit_m;  
    massbal_struct->vol_runoff+=diff;  // send excess water back to GIUH runoff  edit FLO
    massbal_struct->vol_infilt-=diff;  // correct overprediction of infilt.      edit FLO
    direct_output_runoff_m+=diff; // bug found by Nels.  This was missing and fixes it. 
    }

  massbal_struct->vol_to_soil              += infiltration_depth_m; 
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
    massbal_struct->vol_runoff+=diff;  // send excess water back to GIUH runoff
    massbal_struct->vol_infilt-=diff;  // correct overprediction of infilt.
    }
  
    
  massbal_struct->vol_to_gw                +=flux_perc_m;
  massbal_struct->vol_soil_to_gw           +=flux_perc_m;
  
  gw_reservoir_struct->storage_m   += flux_perc_m;
  soil_reservoir_struct->storage_m -= flux_perc_m;
  soil_reservoir_struct->storage_m -= flux_lat_m;
  massbal_struct->vol_soil_to_lat_flow     += flux_lat_m;  //TODO add this to nash cascade as input
  massbal_struct->volout=massbal_struct->volout+flux_lat_m;
  
  conceptual_reservoir_flux_calc(gw_reservoir_struct,&primary_flux,&secondary_flux);
  
  
  flux_from_deep_gw_to_chan_m=primary_flux;  // m/h   <<<<<<<<<< BASE FLOW FLUX >>>>>>>>>
  if(flux_from_deep_gw_to_chan_m >  gw_reservoir_struct->storage_m)  {
  flux_from_deep_gw_to_chan_m=gw_reservoir_struct->storage_m;
  // TODO: set a flag when flux larger than storage
  printf("WARNING: Groundwater flux larger than storage \n");
  }
  	
 
  massbal_struct->vol_from_gw+=flux_from_deep_gw_to_chan_m;
  
  // in the instance of calling the gw reservoir the secondary flux should be zero- verify
  if(is_fabs_less_than_epsilon(secondary_flux,1.0e-09)==FALSE) printf("problem with nonzero flux point 1\n");

  
  // adjust state of deep groundwater conceptual nonlinear reservoir
  //-----------------------------------------------------------------
 
  gw_reservoir_struct->storage_m -= flux_from_deep_gw_to_chan_m;

  
  // Solve the convolution integral ffor this time step 

  /* xinanjiang_dev
  giuh_runoff_m = convolution_integral(Schaake_output_runoff_m,num_giuh_ordinates,    */
  giuh_runoff_m = convolution_integral(direct_output_runoff_m,num_giuh_ordinates,
                                       giuh_ordinates_arr,runoff_queue_m_per_timestep_arr);
  massbal_struct->vol_out_giuh+=giuh_runoff_m;

  massbal_struct->volout+=giuh_runoff_m;
  massbal_struct->volout+=flux_from_deep_gw_to_chan_m;
  
  // Route lateral flow through the Nash cascade.
  nash_lateral_runoff_m = nash_cascade(flux_lat_m,num_lateral_flow_nash_reservoirs,
                                       K_nash,nash_storage_arr);
  massbal_struct->vol_in_nash   += flux_lat_m;
  massbal_struct->vol_out_nash  += nash_lateral_runoff_m;

#ifdef DEBUG
        fprintf(out_debug_fptr,"%d %lf %lf  %lf\n",tstep,flux_lat_m,nash_lateral_runoff_m,flux_from_deep_gw_to_chan_m);
#endif

  Qout_m = giuh_runoff_m + nash_lateral_runoff_m + flux_from_deep_gw_to_chan_m;
    
    // #### COPY BACK STATE VALUES BY POINTER REFERENCE SO VISIBLE TO FRAMEWORK    ####    
    *soil_reservoir_storage_deficit_m_ptr = soil_reservoir_storage_deficit_m;

    /* xinanjiang_dev
    *Schaake_output_runoff_m_ptr          = Schaake_output_runoff_m;    */
    *flux_output_direct_runoff_m          = direct_output_runoff_m;

    *infiltration_depth_m_ptr             = infiltration_depth_m;
    *flux_perc_m_ptr                      = flux_perc_m;
    *flux_lat_m_ptr                       = flux_lat_m;
    *gw_reservoir_storage_deficit_m_ptr   = gw_reservoir_storage_deficit_m;
    *flux_from_deep_gw_to_chan_m_ptr      = flux_from_deep_gw_to_chan_m;
    *giuh_runoff_m_ptr                    = giuh_runoff_m;
    *nash_lateral_runoff_m_ptr            = nash_lateral_runoff_m;
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
  /*printf("primary_flux_m %lf  \n",*primary_flux_m);
  printf("da_reservoir->coeff_primary %lf \n",da_reservoir->coeff_primary);
  printf("da_reservoir->exponent_primary %lf \n",da_reservoir->exponent_primary);
  printf("da_reservoir->storage_m %lf \n",da_reservoir->storage_m);
  printf("da_reservoir->storage_max_m %lf \n",da_reservoir->storage_max_m);*/
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
double timestep_d,Schaake_parenthetical_term,Ic,Px;

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


//##############################################################
//########   XINANJIANG RUNOFF PARTITIONING SCHEME   ###########
//##############################################################

void Xinanjiang_partitioning_scheme(double water_input_depth_m, double field_capacity_m,
                                    double max_soil_moisture_storage_m, double column_total_soil_water_m,
                                    struct direct_runoff_parameters_structure *parms, 
                                    double *surface_runoff_depth_m, double *infiltration_depth_m)
{
  //------------------------------------------------------------------------
  //  This module takes the water_input_depth_m and separates it into surface_runoff_depth_m
  //  and infiltration_depth_m by calculating the saturated area and runoff based on a scheme developed
  //  for the Xinanjiang model by Jaywardena and Zhou (2000). According to Knoben et al.
  //  (2019) "the model uses a variable contributing area to simulate runoff.  [It] uses
  //  a double parabolic curve to simulate tension water capacities within the catchment, 
  //  instead of the original single parabolic curve" which is also used as the standard 
  //  VIC fomulation.  This runoff scheme was selected for implementation into NWM v3.0.
  //  REFERENCES:
  //  1. Jaywardena, A.W. and M.C. Zhou, 2000. A modified spatial soil moisture storage 
  //     capacity distribution curve for the Xinanjiang model. Journal of Hydrology 227: 93-113
  //  2. Knoben, W.J.M. et al., 2019. Supplement of Modular Assessment of Rainfall-Runoff Models
  //     Toolbox (MARRMoT) v1.2: an open-source, extendable framework providing implementations
  //     of 46 conceptual hydrologic models as continuous state-space formulations. Supplement of 
  //     Geosci. Model Dev. 12: 2463-2480.
  //-------------------------------------------------------------------------
  //  Written by RLM May 2021
  //  Adapted by JMFrame September 2021 for new version of CFE
  //-------------------------------------------------------------------------
  // Inputs
  //   double  water_input_depth_m           amount of water input to soil surface this time step [m]
  //   double  field_capacity_m              amount of water stored in soil reservoir when at field capacity [m]
  //   double  max_soil_moisture_storage_m   total storage of the soil moisture reservoir (porosity*soil thickness) [m]
  //   double  column_total_soil_water_m     current storage of the soil moisture reservoir [m]
  //   double  a_inflection_point_parameter  a parameter
  //   double  b_shape_parameter             b parameter
  //   double  x_shape_parameter             x parameter
  //
  // Outputs
  //   double  surface_runoff_depth_m        amount of water partitioned to surface water this time step [m]
  //   double  infiltration_depth_m          amount of water partitioned as infiltration (soil water input) this time step [m]
  //------------------------------------------------------------------------- 

  double tension_water_m, free_water_m, max_tension_water_m, max_free_water_m, pervious_runoff_m;

  //could move this if statement outside of both the Schaake and Xinanjiang subroutines  edit FLO- moved to main().

  // partition the total soil water in the column between free water and tension water
  free_water_m = column_total_soil_water_m - field_capacity_m;

  if(0.0 < free_water_m) {      //edit FLO
    tension_water_m = field_capacity_m;
  } else {
    free_water_m = 0.0;
    tension_water_m = column_total_soil_water_m;
  }

  // estimate the maximum free water and tension water available in the soil column
  max_free_water_m = max_soil_moisture_storage_m - field_capacity_m;
  max_tension_water_m = field_capacity_m;

  // check that the free_water_m and tension_water_m do not exceed the maximum and if so, change to the max value
  if(max_free_water_m < free_water_m) free_water_m = max_free_water_m;
  if(max_tension_water_m < tension_water_m) tension_water_m = max_tension_water_m;

  // NOTE: the impervious surface runoff assumptions due to frozen soil used in NWM 3.0 have not been included.
  // We are assuming an impervious area due to frozen soils equal to 0 (see eq. 309 from Knoben et al).

  // The total (pervious) runoff is first estimated before partitioning into surface and subsurface components.
  // See Knoben et al eq 310 for total runoff and eqs 313-315 for partitioning between surface and subsurface
  // components.

  // Calculate total estimated pervious runoff. 
  // NOTE: If the impervious surface runoff due to frozen soils is added,
  // the pervious_runoff_m equation will need to be adjusted by the fraction of pervious area.
  if ((tension_water_m/max_tension_water_m) <= (0.5 - parms->a_Xinanjiang_inflection_point_parameter)) {
      pervious_runoff_m = water_input_depth_m * (pow((0.5 - parms->a_Xinanjiang_inflection_point_parameter), 
                                                     (1.0 - parms->b_Xinanjiang_shape_parameter)) *
                                                 pow((1.0 - (tension_water_m/max_tension_water_m)),
                                                     parms->b_Xinanjiang_shape_parameter));

  } else {
      pervious_runoff_m = water_input_depth_m * (1.0 - pow((0.5 + parms->a_Xinanjiang_inflection_point_parameter), 
                                                         (1.0 - parms->b_Xinanjiang_shape_parameter)) * 
                                                     pow((1.0 - (tension_water_m/max_tension_water_m)),
                                                         (parms->b_Xinanjiang_shape_parameter)));
  }
  // Separate the surface water from the pervious runoff 
  // NOTE: If impervious runoff is added to this subroutine, impervious runoff should be added to
  // the surface_runoff_depth_m.
  *surface_runoff_depth_m = pervious_runoff_m * (1.0 - pow((1.0 - (free_water_m/max_free_water_m)),parms->x_Xinanjiang_shape_parameter));

  // The surface runoff depth is bounded by a minimum of 0 and a maximum of the water input depth.
  // Check that the estimated surface runoff is not less than 0.0 and if so, change the value to 0.0.
  if(*surface_runoff_depth_m < 0.0) *surface_runoff_depth_m = 0.0;
  // Check that the estimated surface runoff does not exceed the amount of water input to the soil surface.  If it does,
  // change the surface water runoff value to the water input depth.
  if(*surface_runoff_depth_m > water_input_depth_m) *surface_runoff_depth_m = water_input_depth_m;
  // Separate the infiltration from the total water input depth to the soil surface.
  *infiltration_depth_m = water_input_depth_m - *surface_runoff_depth_m;    

return;
}


//##############################################################
//####################   ET FROM RAINFALL   ####################
//##############################################################
void et_from_rainfall(double *timestep_rainfall_input_m, struct evapotranspiration_structure *et_struct)
{
    /*
        iff it is raining, take PET from rainfall first.  Wet veg. is efficient evaporator.
    */

    if (*timestep_rainfall_input_m >0.0){

        if (*timestep_rainfall_input_m > et_struct->potential_et_m_per_timestep){
    
            et_struct->actual_et_m_per_timestep = et_struct->potential_et_m_per_timestep;
            *timestep_rainfall_input_m -= et_struct->actual_et_m_per_timestep;
        }

        else{
            // LKC: This was incorrectly set to potential instead of actual
	    et_struct->actual_et_m_per_timestep = *timestep_rainfall_input_m;          
            *timestep_rainfall_input_m = 0.0;

        }
    // Move this out of the loop since EVPT needs to be corrected wheter R > EVPT or not
    et_struct->potential_et_m_per_timestep = et_struct->potential_et_m_per_timestep-et_struct->actual_et_m_per_timestep;
    }
}

//##############################################################
//####################   ET FROM SOIL   ########################
//##############################################################
void et_from_soil(struct conceptual_reservoir *soil_res, struct evapotranspiration_structure *et_struct, struct NWM_soil_parameters *soil_parms)
{
    /*
        take AET from soil moisture storage, 
        using Budyko type function to limit PET if wilting<soilmoist<field_capacity
    */
    double Budyko_numerator;
    double Budyko_denominator;
    double Budyko;
    
    if (et_struct->potential_et_m_per_timestep > 0){
        
        if (soil_res->storage_m >= soil_res->storage_threshold_primary_m){        
            et_struct->actual_et_m_per_timestep = min(et_struct->potential_et_m_per_timestep, soil_res->storage_m);
            soil_res->storage_m -= et_struct->actual_et_m_per_timestep;
	    // Modified - no reason to set potential et to zero
            et_struct->potential_et_m_per_timestep = et_struct->potential_et_m_per_timestep - et_struct->actual_et_m_per_timestep;
        }                 
        else if (soil_res->storage_m > soil_parms->wilting_point_m && soil_res->storage_m < soil_res->storage_threshold_primary_m){
        
            Budyko_numerator = soil_res->storage_m - soil_parms->wilting_point_m;
            Budyko_denominator = soil_res->storage_threshold_primary_m - soil_parms->wilting_point_m;
            Budyko = Budyko_numerator / Budyko_denominator;
                           
            et_struct->actual_et_m_per_timestep = Budyko * (et_struct->potential_et_m_per_timestep);
                           
            soil_res->storage_m -= et_struct->actual_et_m_per_timestep;
            et_struct->potential_et_m_per_timestep = et_struct->potential_et_m_per_timestep - et_struct->actual_et_m_per_timestep;
        }
    }
}

extern int is_fabs_less_than_epsilon(double a,double epsilon)  // returns true if fabs(a)<epsilon
{
if(fabs(a)<epsilon) return(TRUE);
else                return(FALSE);
}

