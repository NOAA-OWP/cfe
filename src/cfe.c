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

#include "cfe.h"
#include "discrete_soil_moisture.h"

#define max(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a > _b ? _a : _b; })
#define min(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a < _b ? _a : _b; })

// mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
//  CFE STATE SPACE FUNCTION
//   Conceptual Functional Equivalent (CFE) to WRF-Hydro Based National Water Model (Version 3.1 and earlier).
//   Conceived by Fred Ogden, May 2020.
//   - Adapted to BMI by Jonathan Frame and Jessica Garrett, May 2021.
//   - Re-written in state-space form July, 2021.
//   - Modified to optionally replace surface routing GIUH with Nash Cascade, adding the ability to simulate
//     retention depth and runon/channel infiltration by Ahmad Jan Khattak, May 2024.
//   - Refactored to use the same Nash Cascade routing ffor both surface and subsurface routing, FLO 6/25
//   - Retention depth and runon/channel infiltration added more parameters and no more model skill.
//     Nash cascade surface routing scheme, retention, and infiltration losses from surface Nash cascade
//     removed by FLO 5/26. (the sake of parsimony)
//
//  Description:
//    The CFE model is a simplified conceptual model designed to emulate as closely as possible the
//    WRF-Hydro based NOAA National Water Model (NWM).  This version of the CFE model is a stateless
//    BMI (Basic Model Interface) implementation. The CFE model was designed to run at scales less than about
//    10 sq. km. using an hourly timestep.  The user is cautioned not to use it at larger scales or time
//    steps significantly different from 1 hour.  Assumptions include spatially uniform precipitation over
//    the entire catchment at the hourly time scale.  At larger space scales (>>10 km) and longer (or shorter)
//    time scales, this assumption will likely be invalid, particularly shorter time scales.
//
//    The CFE model structure does the following:
//
//    - partition rainfall into soil moisture or surface water using either:
//        - Schaake (1995) method
//              -or-
//        - Xinanjiang method as modified by Jayawardena and Zhou (2012)
//    - manage soil moisture storage as a single conceptual reservoir,
//    - pass soil water to deep groundwater and to lateral subsurface flow using a field capacity threshold
//    - route surface water to the catchment outlet using:
//        - a geomorphic instantaneous unit hydrograph, or other unit hydrograph
//    - route shallow lateral subsurface flow to the catchment outlet using 2-reservoir Nash cascade
//    - manage deep groundwater storage and simulate base flow
//
//  The CFE model uses methods identical to the WRF-Hydro based NWM (versions 3.1 and earlier) to simulate
//  the following processes:
//    - rainfall partitioning
//    - percolation to groundwater and activation of lateral subsurface flow
//    - deep groundwater storage accounting and base flow simulation
//  ----------------------------------------------------------------
//  Differences between CFE and WRF-Hydro based National Water Model
//    - Noah-MP solves the soil moisture form of the Richardson/Richards' equation using a coarse
//      discretization (0.1, 0.3, 0.6, 1.0 m) of a homogeneous, (non-layered) soil.  The CFE
//      model replaces this with a single conceptual linear reservoir having two outlets, both using a
//      field capacity threshold, as was applied in the WRF-Hydro based NWM, to activate percolation to
//      groundwater and lateral subsurface flow.
//    - CFE eliminates the quasi-2D (steepest descent) overland flow and lateral subsurface flow routing on
//      a 250m grid with appropriate small catchment conceptualizations (unit hydrograph, Nash cascade).
//      This was done because: (1) a 250 m grid does not accurately represent the pdf of land slope (the
//      primary driver of lateral hydrologic fluxes), (2) the 250 m routing grid creates a signficant scale
//      mismatch between parameter identifiability and representativeness, and results in large area-average
//      "effective" parameter values and breaks the tie with physics, and (3) the 250 m grid routing
//      requires a lot of unnecessary spatial parameter estimation and massive computational power.
//
//   Version history:
//   2.1  July, 2025, Modified to read/parse keyword based input config file, prototyped model definition file.
//   3.0  August, 2025, Modified to include discretized stateless soil moisture balance model that solves
//        soil moisture at Noah-MP discretizations (top-down) of 0.1, 0.3, 0.6, and 1.0 m
// mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
extern void cfe(
    double *soil_reservoir_storage_deficit_m_ptr,
    struct NWM_SOIL_PARAMETERS_STRUCTURE NWM_soil_params_struct,
    struct CONCEPTUAL_RESERVOIR_STRUCTURE *soil_reservoir_struct,
    double timestep_h,
    struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE infiltration_excess_params_struct,
    double timestep_rainfall_input_m,
    double *surface_runoff_m_ptr,
    double *infiltration_depth_m_ptr,
    double *flux_perc_m_ptr,
    double *flux_lat_m_ptr,
    double *gw_reservoir_storage_deficit_m_ptr,
    struct CONCEPTUAL_RESERVOIR_STRUCTURE *gw_reservoir_struct,
    double *flux_from_deep_gw_to_chan_m_ptr,
    double *flux_direct_runoff_to_channel_m_ptr,
    int num_giuh_ordinates,
    double *giuh_ordinates_arr,
    double *giuh_runoff_queue_m_per_timestep_arr,
    double *flux_nash_subsurface_lateral_runoff_m_ptr,
    struct NASH_CASCADE_PARAMETERS_STRUCTURE *nash_subsurface_params,
    struct EVAPOTRANSPIRATION_STRUCTURE *evap_struct,
    double *Qout_m_ptr,
    cfe_volbal_struct *volbal_struct,
    double time_step_seconds,
    int surface_runoff_scheme,
    int yes_simulate_discrete_soil_moisture,
    SoilControl *soil_control,
    SoilGeometry *soil_geometry,
    SoilParameters *soil_parameters,
    SoilStateIn *soil_state_in,
    SoilStateOut *soil_state_out,
    SoilFluxes *soil_fluxes,
    SoilLookupTables *ch_lookup_tables,
    double *soil_discrete_storage_theta,
    int verbosity)
{
  // #######################################################################
  // CFE STATE SPACE FUNCTION
  // #######################################################################

  // ####    COPY THE MODEL FUNCTION STATE SPACE TO LOCAL VARIABLES    ####
  // ####    Reason: so we don't have to re-write domain science code to de-reference a whole bunch of pointers
  // ####    Note: all of thes variables are storages in [m] or fluxes in [m/timestep]

  double soil_reservoir_storage_deficit_m = *soil_reservoir_storage_deficit_m_ptr;           // storage [m]
  double flux_surface_runoff_input_to_surface_routing_m = *surface_runoff_m_ptr;             // Schaake/Xinanjiang partitioned excess water this timestep [m]*/
  double infiltration_depth_m = *infiltration_depth_m_ptr;                                   // Schaake partitioned infiltration this timestep [m]
  double impervious_runoff_m = 0.0;                                                          // runoff from mapped impervious area [m]
  double flux_perc_soil_to_gw_m = *flux_perc_m_ptr;                                          // water moved from soil reservoir to gw reservoir this timestep [m]
  double flux_soil_to_subsurface_lat_m = *flux_lat_m_ptr;                                    // water moved from soil reservoir to lateral flow Nash cascad this timestep [m]
  double gw_reservoir_storage_deficit_m = *gw_reservoir_storage_deficit_m_ptr;               // deficit in gw reservoir storage [m]
  double flux_from_deep_gw_to_chan_m = *flux_from_deep_gw_to_chan_m_ptr;                     // water moved from gw reservoir to catchment outlet nexus this timestep [m]
  double flux_direct_runoff_to_channel_m = *flux_direct_runoff_to_channel_m_ptr;             // water leaving GIUH or NASH cascade reservoir to outlet this timestep [m]
  double flux_nash_subsurface_lateral_runoff_m = *flux_nash_subsurface_lateral_runoff_m_ptr; // water leaving lateral subsurface flow Nash cascade this timestep [m]
  double Qout_m = *Qout_m_ptr;                                                               // the total runoff this timestep (Surface+Nash+GW) [m]

  // LOCAL VARIABLES, the values of which are not important to describe the model state.  They are like notes on scrap paper.

  double diff = 0.0;
  double primary_flux = 0.0;   // pointers to these variables passed to conceptual nonlinear reservoir which has two outlets, primary & secondary
  double secondary_flux = 0.0; // pointers to these variables passed to conceptual nonlinear reservoir which has two outlets, primary & secondary
  double lateral_flux = 0.0;   // flux from soil to lateral flow Nash cascade +to cascade  [m/timestep]

  double percolation_flux = 0.0; // flux from soil to gw nonlinear researvoir, +downward  [m/timestep]
  double soil_storage_temp_m = 0.0;

  // IMPORTANT NOTE: The DSMB option accounts ffor the infiltration from either
  // the Schaake or Xinanjiang partitioning
  // function that won't fit into the upper disc of the soil.  In general
  // it should almost always be zero, because the Discrete Soil Balance Model
  // is strictly NOT an infiltration excess code.  The discretization is too coarse
  // plus it's not accurately solving Richardson/Richards or SMVE.

  // store current soil storage_m in a temp. variable
  if (yes_simulate_discrete_soil_moisture)
  { // simulating discretized soil
    for (int i = 0; i < NDISC; i++)
      soil_storage_temp_m += soil_discrete_storage_theta[i] * soil_geometry->dz_m[i];
  }
  else
  { // conceptual soil reservoir
    soil_storage_temp_m = soil_reservoir_struct->storage_m;
  }

#if CFE_DEBUG > 1
  printf("CFE KERNEL DEBUG:\n");
  printf("  GIUH num ordinates: %d\n", num_giuh_ordinates);
  if (num_giuh_ordinates > 0)
    printf("  GIUH ordinate[0]: %.6f\n", giuh_ordinates_arr[0]);
  printf("  Input to surface routing: %.6f\n", flux_surface_runoff_input_to_surface_routing_m);
#endif

  // ET Demand
  evap_struct->potential_et_m_per_timestep = evap_struct->potential_et_m_per_s * time_step_seconds;
  evap_struct->reduced_potential_et_m_per_timestep = evap_struct->potential_et_m_per_s * time_step_seconds;

  if (isnan(timestep_rainfall_input_m)) {
    fprintf(stderr, "WARNING: NaN rainfall input — treating as 0.0\n");
    timestep_rainfall_input_m = 0.0;
  }

  evap_struct->actual_et_from_rain_m_per_timestep = 0.0;
  if (timestep_rainfall_input_m > 0)
  {
    // calculate evaporation from rainfall
    et_from_rainfall(&timestep_rainfall_input_m, evap_struct);
  }

  /* vol_et_from_rain: tracks ET by source (rain interception)
   * vol_et_to_atm:    tracks total ET leaving domain (all sources)
   * These are intentionally parallel accounting — not double-counting. */
  volbal_struct->vol_et_from_rain += evap_struct->actual_et_from_rain_m_per_timestep;
  volbal_struct->vol_et_to_atm    += evap_struct->actual_et_from_rain_m_per_timestep;
  volbal_struct->volout           += evap_struct->actual_et_from_rain_m_per_timestep;

  evap_struct->actual_et_from_soil_m_per_timestep = 0.0;

  if (!yes_simulate_discrete_soil_moisture)
  {
    // CONCEPTUAL RESERVOIR: Take remaining ET demand from soil, with Budyko-style limiter
    evap_struct->actual_et_from_soil_m_per_timestep = 0.0;
    if (soil_reservoir_struct->storage_m > NWM_soil_params_struct.wilting_point_m)
    {
      et_from_soil(soil_reservoir_struct, evap_struct, &NWM_soil_params_struct);
    }
    soil_reservoir_storage_deficit_m = (NWM_soil_params_struct.smcmax * NWM_soil_params_struct.D - soil_reservoir_struct->storage_m);
  }
// are we missing an else condition here, or in the case of DSBM is ET demand taken from soil elsewhere in the code? -FLO

  if (0.0 < timestep_rainfall_input_m) // it's raining
  {
    if (infiltration_excess_params_struct.surface_water_partitioning_scheme == PARTITION_SCHAAKE)
    {

      if (yes_simulate_discrete_soil_moisture)
      {

        double column_total_soil_moisture_deficit_m = 0.0;
        double column_total_soil_moisture_m = 0.0;
        double field_capacity_storage_m = soil_parameters->theta_fc * NWM_soil_params_struct.D;
        for (int i = 0; i < NDISC; i++)
        {
          column_total_soil_moisture_m += soil_state_in->theta_in[i] * soil_geometry->dz_m[i];
        }
        column_total_soil_moisture_deficit_m = NWM_soil_params_struct.smcmax * NWM_soil_params_struct.D -
                                               column_total_soil_moisture_m;

        if (!soil_control->is_sft_coupled)
          soil_state_in->ice_fraction = 0.0;

        // ##################################################
        //  partition rainfall using Schaake scheme
        // ##################################################
        Schaake_partitioning_scheme(
            timestep_h,
            field_capacity_storage_m,
            infiltration_excess_params_struct.Schaake_adjusted_magic_constant_by_soil_type,
            column_total_soil_moisture_deficit_m,
            timestep_rainfall_input_m,
            NWM_soil_params_struct.smcmax,
            NWM_soil_params_struct.D,
            infiltration_excess_params_struct.urban_decimal_fraction,
            &flux_surface_runoff_input_to_surface_routing_m,
            &infiltration_depth_m,
            &impervious_runoff_m,
            soil_state_in->ice_fraction,
            infiltration_excess_params_struct.ice_content_threshold);
      }

      if (!yes_simulate_discrete_soil_moisture)
      {

        if (!soil_reservoir_struct->is_sft_coupled)
          soil_reservoir_struct->ice_fraction_schaake = 0.0; // to ensure that ice_fraction_schaake is set to 0.0 ffor uncoupled SFT

        Schaake_partitioning_scheme(
            timestep_h,
            soil_reservoir_struct->storage_threshold_primary_m,
            infiltration_excess_params_struct.Schaake_adjusted_magic_constant_by_soil_type,
            soil_reservoir_storage_deficit_m,
            timestep_rainfall_input_m,
            NWM_soil_params_struct.smcmax,
            NWM_soil_params_struct.D,
            infiltration_excess_params_struct.urban_decimal_fraction,
            &flux_surface_runoff_input_to_surface_routing_m,
            &infiltration_depth_m,
            &impervious_runoff_m,
            soil_reservoir_struct->ice_fraction_schaake,
            infiltration_excess_params_struct.ice_content_threshold);
      }
    }
    else if (infiltration_excess_params_struct.surface_water_partitioning_scheme == PARTITION_XINANJIANG)
    {
      // ##################################################
      //  partition rainfall using Xinanjiang scheme
      // ##################################################

      if (yes_simulate_discrete_soil_moisture)
      {

        if (!soil_control->is_sft_coupled)
          soil_state_in->ice_fraction = 0.0;

        double thickness_of_upper_two_discs_m          = soil_geometry->dz_m[0] + soil_geometry->dz_m[1];
        double upper_two_disc_field_capacity_storage_m = soil_parameters->theta_fc * thickness_of_upper_two_discs_m;
        double upper_two_disc_max_soil_moisture_m      = NWM_soil_params_struct.smcmax * thickness_of_upper_two_discs_m;
        double upper_two_disc_total_soil_moisture_m    = 0.0;
        for (int i = 0; i < 2; i++)
        {
          upper_two_disc_total_soil_moisture_m += soil_state_in->theta_in[i] * soil_geometry->dz_m[i];
        }

        Xinanjiang_partitioning_scheme(
            timestep_rainfall_input_m,
            upper_two_disc_field_capacity_storage_m,
            upper_two_disc_max_soil_moisture_m,
            upper_two_disc_total_soil_moisture_m,
            &infiltration_excess_params_struct,
            &flux_surface_runoff_input_to_surface_routing_m,
            &infiltration_depth_m,
            &impervious_runoff_m,
            soil_reservoir_struct->ice_fraction_xinanjiang);
      }

      if (!yes_simulate_discrete_soil_moisture)
      {

        // to ensure that ice_fraction_xinanjiang is set to 0.0 ffor uncoupled SFT
        if (!soil_reservoir_struct->is_sft_coupled)
          soil_reservoir_struct->ice_fraction_xinanjiang = 0.0;

        Xinanjiang_partitioning_scheme(
            timestep_rainfall_input_m,
            soil_reservoir_struct->storage_threshold_primary_m,
            soil_reservoir_struct->storage_max_m,
            soil_reservoir_struct->storage_m,
            &infiltration_excess_params_struct,
            &flux_surface_runoff_input_to_surface_routing_m,
            &infiltration_depth_m,
            &impervious_runoff_m,
            soil_reservoir_struct->ice_fraction_xinanjiang);
      }
    }
    else
    {
      fprintf(stderr, "Problem, must specify one of Schaake of Xinanjiang partitioning scheme.\n");
      fprintf(stderr, "Program terminating.\n");
      exit(-1); // note -1 is arbitrary
    }
  }
  else
  { //  it's not raining
    flux_surface_runoff_input_to_surface_routing_m = 0.0;
    infiltration_depth_m = 0.0;
  }

  // Iff running a single conceptual soil balance reservoir model- check to make
  // sure that there is storage available in soil to hold the water that the Schaake
  // or Xinanjiang function says infiltrated...
  //--------------------------------------------------------------------------------------------------
  if (!yes_simulate_discrete_soil_moisture)
  {
    // CONCEPTUAL SOIL RESERVOIR
    if (soil_reservoir_storage_deficit_m < infiltration_depth_m)
    {
      // OVERFLOW CASE: infiltration calculated by Schaake or Xinanjiang doesn't all fit
      double excess = infiltration_depth_m - soil_reservoir_storage_deficit_m;
      flux_surface_runoff_input_to_surface_routing_m += excess;
      volbal_struct->vol_direct_runoff += excess;
      infiltration_depth_m = soil_reservoir_storage_deficit_m;
      soil_reservoir_struct->storage_m = soil_reservoir_struct->storage_max_m;
      soil_reservoir_storage_deficit_m = 0.0;
    }
    else
    {
      // NORMAL CASE: all infiltration fits
      soil_reservoir_struct->storage_m += infiltration_depth_m;
    }
  }

  // Don't doo the above check in the discrete soil moisture situation.  It handles that internally.
  if (yes_simulate_discrete_soil_moisture)
  {
    // Use the discrete soil balance module (DSBM)
    SoilForcing soil_forcing = {0};
    soil_forcing.rain_mm_per_h = infiltration_depth_m * 1000.0 / timestep_h; // Convert m/timestep to mm/h
    soil_forcing.pet_mm_per_h = evap_struct->reduced_potential_et_m_per_timestep * 1000.0 / timestep_h;

    TimestepSoilVolbal soil_volbal = {0}; // initializes all elements to 0

    // soil_state_in->total_storage_m  is initialized at run time in initialize()

    soil_state_in->storage_deficit_m = (NWM_soil_params_struct.smcmax * NWM_soil_params_struct.D) -
                                       soil_state_in->total_storage_m;

    // call discrete soil balance module to solve soil moistures in 0.1, 0.3, 0.6 and 1.0 m discs.
    DSBM_step_one_hour_stateless(
        soil_control,
        soil_geometry,
        soil_parameters,
        ch_lookup_tables, // Clapp-Hornberger look-up-table, may be NULL
        soil_state_in,
        evap_struct, // Pass the actual ET struct
        &soil_forcing,
        soil_state_out,
        soil_fluxes,
        &soil_volbal,
        NULL);

    if (verbosity > 2 && yes_simulate_discrete_soil_moisture)
    { // THIS PRODUCES A TON OF OUTPUT TO stdout
      check_dsbm_local_volume_balance(
          &soil_volbal,
          soil_fluxes,
          soil_state_in,
          soil_state_out,
          soil_geometry,
          infiltration_depth_m,
          evap_struct->actual_et_from_soil_m_per_timestep +
              evap_struct->actual_bare_soil_evaporation_m_per_timestep,
          1.0e-9 // tolerance in meters  1.0e-9 is a good value.  Set to 1.0e-19 to print all the time
      );
    }

    // EXTRACT RESULTS AND CONNECT TO VOLUME BALANCE
    flux_perc_soil_to_gw_m = soil_volbal.perc_m;
    flux_soil_to_subsurface_lat_m = soil_volbal.lateral_m;

    // Extract infiltration excess from discrete soil model
    double DSBM_infiltration_excess_runoff_m = soil_fluxes->rain_excess_m;
    volbal_struct->vol_direct_runoff += soil_fluxes->rain_excess_m; // water that won't fit into disc0 of DSBM this time step

    // Add this excess to surface routing input
    flux_surface_runoff_input_to_surface_routing_m += DSBM_infiltration_excess_runoff_m;

    // The full pre-DSBM infiltration_depth_m estimate gets added to vol_infilt below
    // even though part of it never entered the soil column - it was just counted a
    // second time as runoff via DSBM_infiltration_excess_runoff_m above.
    infiltration_depth_m -= DSBM_infiltration_excess_runoff_m;
    if (infiltration_depth_m < 0.0) infiltration_depth_m = 0.0;

    // Update discrete soil storage array ffor persistence
    for (int i = 0; i < NDISC; i++)
    {
      soil_discrete_storage_theta[i] = soil_state_out->theta_out[i];
    }

    // Update soil state input ffor next timestep
    for (int i = 0; i < NDISC; i++)
    {
      soil_state_in->theta_in[i] = soil_state_out->theta_out[i];
    }

    // Volume balance accounting
    volbal_struct->vol_soil_start = soil_state_in->total_storage_m; // Initial from input struct
    volbal_struct->vol_soil_end = soil_state_out->total_storage_m;  // Final from output struct

    // Use output deficit ffor any surface routing needs
    soil_reservoir_storage_deficit_m = soil_state_out->storage_deficit_m; // If needed elsewhere
  }
  // MOVED THIS TO AFTER DSBM CALL______________________________________
  {
    double total_soil_et_m =
        evap_struct->actual_et_from_soil_m_per_timestep +
        evap_struct->actual_bare_soil_evaporation_m_per_timestep;
    volbal_struct->vol_forest_aet +=
        evap_struct->actual_et_from_soil_m_per_timestep;
    volbal_struct->vol_bare_soil_evaporation +=
        evap_struct->actual_bare_soil_evaporation_m_per_timestep;
    volbal_struct->vol_et_from_soil += total_soil_et_m;
    volbal_struct->vol_et_to_atm   += total_soil_et_m;
    volbal_struct->volout          += total_soil_et_m;
  }

  evap_struct->actual_et_m_per_timestep = evap_struct->actual_et_from_rain_m_per_timestep +
                                          evap_struct->actual_et_from_soil_m_per_timestep +
                                          evap_struct->actual_bare_soil_evaporation_m_per_timestep;

  //-- NEW DSBM
  // CONSOLIDATED VOLUME BALANCE - SAME LOGIC FOR BOTH MODELS
  volbal_struct->vol_to_soil += infiltration_depth_m;
  volbal_struct->vol_impervious_runoff += impervious_runoff_m;
  volbal_struct->vol_pervious_runoff +=
      flux_surface_runoff_input_to_surface_routing_m - impervious_runoff_m;
  volbal_struct->vol_runoff += flux_surface_runoff_input_to_surface_routing_m;
  volbal_struct->vol_infilt += infiltration_depth_m;

  // Soil fluxes (calculated differently but accounted the same way)
  if (yes_simulate_discrete_soil_moisture)
  {
    // DSBM fluxes already calculated and stored in flux_perc_soil_to_gw_m and flux_soil_to_subsurface_lat_m
    volbal_struct->vol_to_gw += flux_perc_soil_to_gw_m;
    volbal_struct->vol_soil_to_gw += flux_perc_soil_to_gw_m;
    volbal_struct->vol_soil_to_lat_flow += flux_soil_to_subsurface_lat_m;
    volbal_struct->vol_in_subsurf_nash += flux_soil_to_subsurface_lat_m;
    // printf("DEBUG flux_soil_to_subsurface_lat_m %.5f\n", flux_soil_to_subsurface_lat_m);
  }

#ifdef DEBUG

  printf("After direct runoff function: rain:%8.5lf mm  runoff:%8.5lf mm  infiltration:%8.5lf mm  residual:%e m\n",
         timestep_rainfall_input_m * 1000.0, flux_surface_runoff_input_to_surface_routing_m * 1000.0, infiltration_depth_m * 1000.0,
         timestep_rainfall_input_m - flux_surface_runoff_input_to_surface_routing_m - infiltration_depth_m);
#endif

  // calculate fluxes from the soil storage into the deep groundwater (percolation) and to lateral subsurface flow
  //--------------------------------------------------------------------------------------------------------------

  if (!yes_simulate_discrete_soil_moisture)
  {
    // CONCEPTUAL RESERVOIR: existing flux calculation
    conceptual_reservoir_flux_calc(soil_reservoir_struct, &percolation_flux, &lateral_flux);
    flux_perc_soil_to_gw_m = percolation_flux;    // m/h <----  flux of percolation from soil to g.w. reservoir
    flux_soil_to_subsurface_lat_m = lateral_flux; // m/h   <----  flux from soil reservoir into the lateral subsurface flow Nash cascade
  }

  gw_reservoir_storage_deficit_m =
      gw_reservoir_struct->storage_max_m - gw_reservoir_struct->storage_m;

  if (!yes_simulate_discrete_soil_moisture)
  {
    // Cap recharge iff GW storage is nearly full
    if (flux_perc_soil_to_gw_m > gw_reservoir_storage_deficit_m)
    {
      diff = flux_perc_soil_to_gw_m - gw_reservoir_storage_deficit_m;
      flux_perc_soil_to_gw_m = gw_reservoir_storage_deficit_m;
      volbal_struct->vol_runoff += diff;
      volbal_struct->vol_infilt -= diff;
    }
  }

  //--------------------------------------------
  // Always add recharge to groundwater storage
  //--------------------------------------------
  if (!yes_simulate_discrete_soil_moisture)
  {
    volbal_struct->vol_to_gw += flux_perc_soil_to_gw_m;
    volbal_struct->vol_soil_to_gw += flux_perc_soil_to_gw_m;
  }
  gw_reservoir_struct->storage_m += flux_perc_soil_to_gw_m;

  if (!yes_simulate_discrete_soil_moisture)
  {
    soil_reservoir_struct->storage_m -= flux_perc_soil_to_gw_m;
    soil_reservoir_struct->storage_m -= flux_soil_to_subsurface_lat_m;
  }

  //--------------------------------------------
  // Compute baseflow AFTER recharge is applied
  //--------------------------------------------
  conceptual_reservoir_flux_calc(gw_reservoir_struct, &primary_flux, &secondary_flux);

  flux_from_deep_gw_to_chan_m = primary_flux; // baseflow flux
  if (flux_from_deep_gw_to_chan_m > gw_reservoir_struct->storage_m)
  {
    flux_from_deep_gw_to_chan_m = gw_reservoir_struct->storage_m;
    printf("WARNING: Groundwater flux larger than storage\n");
  }

  // === UPDATE STORAGES FIRST ===
  gw_reservoir_struct->storage_m -= flux_from_deep_gw_to_chan_m;

  // === THEN UPDATE VOLUME BALANCES ===
  volbal_struct->vol_from_gw += flux_from_deep_gw_to_chan_m;

  //--------------------------------------------
  // Soil bookkeeping ffor conceptual reservoir
  //--------------------------------------------
  if (!yes_simulate_discrete_soil_moisture)
  {
    volbal_struct->vol_soil_to_lat_flow += flux_soil_to_subsurface_lat_m;
    soil_reservoir_struct->storage_change_m =
        soil_reservoir_struct->storage_m - soil_storage_temp_m;
  }

  if (is_fabs_less_than_epsilon(secondary_flux, 1.0e-09) == FALSE)
    printf("problem with nonzero flux point 1\n");


  if (surface_runoff_scheme == SURF_ROUTE_GIUH)
  { // Solve the convolution integral ffor this time step
    flux_direct_runoff_to_channel_m = giuh_convolution_integral(flux_surface_runoff_input_to_surface_routing_m,
                                                                num_giuh_ordinates,
                                                                giuh_ordinates_arr,
                                                                giuh_runoff_queue_m_per_timestep_arr);
  }

  // Route lateral flow through the Nash cascade (new use of this function ffor lateral subsurface flow routing -FLO 6/25)

  flux_nash_subsurface_lateral_runoff_m = nash_cascade_routing(flux_soil_to_subsurface_lat_m, 0.0, nash_subsurface_params); //<- 0.0 here means no losses.

  // moved here after the lateral subsurface routing call is done.
  volbal_struct->vol_out_surface += flux_direct_runoff_to_channel_m;
  volbal_struct->volout += flux_direct_runoff_to_channel_m;
  volbal_struct->volout += flux_nash_subsurface_lateral_runoff_m; // THIS WAS MISSING.
  volbal_struct->volout += flux_from_deep_gw_to_chan_m;

  if (!yes_simulate_discrete_soil_moisture)
  {
    volbal_struct->vol_in_subsurf_nash += flux_soil_to_subsurface_lat_m;
  }
  volbal_struct->vol_out_subsurf_nash += flux_nash_subsurface_lateral_runoff_m;

  Qout_m = flux_direct_runoff_to_channel_m + flux_nash_subsurface_lateral_runoff_m + flux_from_deep_gw_to_chan_m;

  //----DEBUG
  // printf("CFE KERNEL DEBUG: surface_out=%lf, nash_out=%lf, gw_out=%lf\n",
  //       flux_direct_runoff_to_channel_m, flux_nash_subsurface_lateral_runoff_m, flux_from_deep_gw_to_chan_m);
  // printf("CFE KERNEL DEBUG: Qout_total=%.4f, volout_running=%.4f\n",
  //       Qout_m, volbal_struct->volout);
  //-----

  // #### COPY BACK STATE VALUES BY POINTER REFERENCE SO VISIBLE TO FRAMEWORK    ####
  *soil_reservoir_storage_deficit_m_ptr = soil_reservoir_storage_deficit_m;
  *surface_runoff_m_ptr = flux_surface_runoff_input_to_surface_routing_m;
  *infiltration_depth_m_ptr = infiltration_depth_m;
  *flux_perc_m_ptr = flux_perc_soil_to_gw_m;
  *flux_lat_m_ptr = flux_soil_to_subsurface_lat_m;
  *gw_reservoir_storage_deficit_m_ptr = gw_reservoir_storage_deficit_m;
  *flux_from_deep_gw_to_chan_m_ptr = flux_from_deep_gw_to_chan_m;
  *flux_direct_runoff_to_channel_m_ptr = flux_direct_runoff_to_channel_m;
  *flux_nash_subsurface_lateral_runoff_m_ptr = flux_nash_subsurface_lateral_runoff_m;
  *Qout_m_ptr = Qout_m;

} // END CFE STATE SPACE FUNCTIONS
  // ####################################################################################################
  // ####################################################################################################

// ##############################################################
// #########   SCHAAKE RUNOFF PARTITIONING SCHEME   #############
// ##############################################################
void Schaake_partitioning_scheme(
    double timestep_h,
    double field_capacity_m,
    double Schaake_adjusted_magic_constant_by_soil_type,
    double column_total_soil_moisture_deficit_m,
    double water_input_depth_m,
    double smcmax,
    double soil_depth,
    double catchment_impervious_fraction,
    double *flux_surface_runoff_input_to_surface_routing_m,
    double *infiltration_depth_m,
    double *impervious_runoff_m,
    double ice_fraction_schaake, double ice_content_threshold)
{

  /*! ===============================================================================
    This subtroutine takes water_input_depth_m and partitions it into flux_surface_runoff_input_to_surface_routing_m and
    infiltration_depth_m using the scheme from Schaake et al. 1996.
  ! --------------------------------------------------------------------------------
  ! ! modified by FLO April 2020 to eliminate reference to ice processes,
  ! ! and to de-obfuscate and use descriptive and dimensionally consistent variable names.
  ! ! Frozen soil effects added from Noah-MP by AJK
  ! --------------------------------------------------------------------------------
      IMPLICIT NONE
  ! --------------------------------------------------------------------------------
  ! inputs
    double timestep_h
    double Schaake_adjusted_magic_constant_by_soil_type = C*Ks(soiltype)/Ks_ref, where C=3, and Ks_ref=2.0E-06 m/s
    double column_total_soil_moisture_deficit_m
    double water_input_depth_m  amount of water input to soil surface this time step [m]

  ! outputs
    double flux_surface_runoff_input_to_surface_routing_m      amount of water partitioned to surface water this time step [m]


  --------------------------------------------------------------------------------*/

  assert(ice_fraction_schaake >= 0.0);
  assert(catchment_impervious_fraction >= 0.0);
  assert(catchment_impervious_fraction <= 1.0);

  double timestep_d, Schaake_parenthetical_term, Ic, Px;
  double pervious_water_input_m;
  double pervious_runoff_m;

  *impervious_runoff_m = water_input_depth_m *
                         catchment_impervious_fraction * 0.95;
  pervious_water_input_m = water_input_depth_m - *impervious_runoff_m;
  pervious_runoff_m = 0.0;

  if (0.0 < pervious_water_input_m)
  {
    if (0.0 >= column_total_soil_moisture_deficit_m)
    {
      pervious_runoff_m = pervious_water_input_m;
      *infiltration_depth_m = 0.0;
    }
    else
    {
      // partition time-step total applied water as per Schaake et al. 1996.
      // change from dt in [s] to dt1 in [d] because kdt has units of [d^(-1)]
      timestep_d = timestep_h / 24.0; // timestep_d is the time step in days.

      // calculate the parenthetical part of Eqn. 34 from Schaake et al. Note the magic constant has units of [d^(-1)]

      Schaake_parenthetical_term = (1.0 - exp(-Schaake_adjusted_magic_constant_by_soil_type * timestep_d));

      // From Schaake et al. Eqn. 2., using the column total moisture deficit
      // BUT the way it is used here, it is the cumulative soil moisture deficit in the entire soil profile.
      // "Layer" info not used in this subroutine in noah-mp, except to sum up the total soil moisture storage.
      // NOTE: when column_total_soil_moisture_deficit_m becomes zero, which occurs when the soil column is saturated,
      // then Ic=0, where Ic in the Schaake paper is called the "spatially averaged infiltration capacity",
      // and is defined in Eqn. 12.

      Ic = column_total_soil_moisture_deficit_m * Schaake_parenthetical_term;

      Px = pervious_water_input_m; // Pervious water input to partitioning scheme this time step [m]

      // This is eqn 24 from Schaake et al.  NOTE: this is 0 in the case of a saturated soil column, when Ic=0.
      // Physically happens only if soil has no-flow lower b.c.

      *infiltration_depth_m = (Px * (Ic / (Px + Ic)));

      if (0.0 < (pervious_water_input_m - (*infiltration_depth_m)))
      {
        pervious_runoff_m = pervious_water_input_m - (*infiltration_depth_m);
      }
      else
      {
        pervious_runoff_m = 0.0;
      }
      *infiltration_depth_m = pervious_water_input_m - pervious_runoff_m;
    }
  }
  else
  {
    pervious_runoff_m = 0.0;
    *infiltration_depth_m = 0.0;
  }

  // Impermeable fraction due to frozen soil, taken directly from Noah-MP
  double factor = 1.0;
  // field_capacity_m = SMCREF in NOAH_MP
  //  ice_content_threshold = frzk in NOAH_MP
  // double frzk = 0.15; // Ice content above which soil is impermeable

  if (ice_fraction_schaake > 1.0E-2)
  {
    int cv_frz = 3;                                              // This is a constant parameter in Noah-MP
    double field_capacity = field_capacity_m / soil_depth;       // divide by the reservior depth to get SMCREF [m/m] (unitless)
    double frz_fact = smcmax / field_capacity * (0.412 / 0.468); // Directly from Noah-MP
    double frzx = ice_content_threshold * frz_fact;

    double acrt = cv_frz * frzx / ice_fraction_schaake;
    double sum1 = 1;

    for (int i1 = 1; i1 < cv_frz; i1++)
    {
      int k = 1;
      for (int i2 = i1 + 1; i2 < cv_frz; i2++)
      {
        k *= i2;
      }
      sum1 += pow(acrt, (cv_frz - i1)) / (double)k;
    }

    factor = 1. - exp(-acrt) * sum1;
  }

  *infiltration_depth_m = factor * (*infiltration_depth_m);

  pervious_runoff_m = pervious_water_input_m - *infiltration_depth_m;

  *flux_surface_runoff_input_to_surface_routing_m =
      *impervious_runoff_m + pervious_runoff_m;

  return;
}

// ##############################################################
// ########   XINANJIANG RUNOFF PARTITIONING SCHEME   ###########
// ##############################################################

void Xinanjiang_partitioning_scheme(
    double water_input_depth_m,
    double field_capacity_m,
    double max_soil_moisture_storage_m,
    double column_total_soil_water_m,
    struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE *parms,
    double *flux_surface_runoff_input_to_surface_routing_m,
    double *infiltration_depth_m,
    double *impervious_runoff_m,
    double ice_fraction_xinanjiang)
{
  //------------------------------------------------------------------------
  //  This module takes the water_input_depth_m and separates it into flux_surface_runoff_input_to_surface_routing_m
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
  //  Reviewed by FLO Feb. 2025, compared against Jayawardena & Zhou (2000) and NWM 3.0 code.  Fixed bug and refactored.
  //        Compared against NWM 3.0 and found identical performance within approx. 1e-05 on calculated outputs
  //-------------------------------------------------------------------------
  // Inputs
  //   double  water_input_depth_m           amount of water input to soil surface this time step [m]
  //   double  field_capacity_m              amount of water stored in soil reservoir when at field capacity [m]
  //   double  max_soil_moisture_storage_m   total storage of the soil moisture reservoir (porosity*soil thickness) [m]
  //   double  column_total_soil_water_m     current storage of the soil moisture reservoir [m]
  //   double  a_inflection_point_parameter  a parameter
  //   double  b_shape_parameter             b parameter
  //   double  x_shape_parameter             x parameter
  //   double  urban_decimal_fraction        fraction of land cover in the modeled area that is classified as urban [unitless decimal]
  //   double  ice_fraction_xinanjiang       fraction of top soil discretization that is frozen [unitless decimal]
  //
  // Outputs
  //   double  flux_surface_runoff_input_to_surface_routing_m  amount of water partitioned to surface water this time step [m]
  //   double  infiltration_depth_m          amount of water partitioned as infiltration (soil water input) this time step [m]
  //   double  impervious_runoff_m           runoff from mapped impervious area [m]
  //-------------------------------------------------------------------------

  // local variables
  double tension_water_m;
  double free_water_m;
  double max_tension_water_m;
  double max_free_water_m;
  double water_input_pervious_fraction_m;
  double frozen_soil_runoff_m;
  double xinanjiang_runoff_m;
  double f_over_F; // notation from Jayawardena and Zhou (2000) see Fig 2.

  // first thing, check whether we can just rreturn without calculating anything to save compute
  //----------------------                  NWM variable name
  *flux_surface_runoff_input_to_surface_routing_m = 0.0; // RUNSRF
  *infiltration_depth_m = 0.0;                           // PDDUM
  *impervious_runoff_m = 0.0;

  if (water_input_depth_m < 1.0e-08)
  { // zero or really close to zero so calculations not needed
    return;
  }

  // calculations required

  // initialize                          // NWM variable name
  tension_water_m = 0.0;     // WM
  max_tension_water_m = 0.0; // WM_MAX
  free_water_m = 0.0;        // SM
  max_free_water_m = 0.0;    // SM_MAX
  frozen_soil_runoff_m = 0.0;
  xinanjiang_runoff_m = 0.0;

  // Partition the total soil water in the column between free water and tension water assuming that
  // total pore space in the soil Vtot, given by (porosity * soil_thickness) also equals max_tension_water_m + max_free_water_m.
  //
  // Water input to soil fills up the tension water storage first.  Once it is full, any additiona
  // input of water goes into free water storage.  This means that iff there is _any_ free water then
  // the tension water storage is full.

  if (column_total_soil_water_m - field_capacity_m > 0.0)
  { // soil moisture greater than field capacity
    free_water_m = column_total_soil_water_m - field_capacity_m;
    tension_water_m = field_capacity_m;
  }
  else
  {
    tension_water_m = column_total_soil_water_m;
  }
  max_tension_water_m = field_capacity_m;
  max_free_water_m = max_soil_moisture_storage_m - field_capacity_m;

  if (tension_water_m > max_tension_water_m)
    tension_water_m = max_tension_water_m; // as done in NWM - could cause volbal err.
  if (free_water_m > max_free_water_m)
    free_water_m = max_free_water_m; // as done in NWM - could cause volbal err.

  // Route 95 percent of rainfall on the mapped impervious fraction directly
  // to runoff.  The remaining rainfall is passed to the Xinanjiang scheme.
  *impervious_runoff_m = water_input_depth_m *
                         parms->urban_decimal_fraction * 0.95;
  water_input_pervious_fraction_m =
      water_input_depth_m - *impervious_runoff_m;

  // Frozen-soil runoff is generated over the remaining area and is tracked as
  // pervious runoff, not mapped-impervious runoff.
  frozen_soil_runoff_m =
      water_input_pervious_fraction_m * ice_fraction_xinanjiang;
  water_input_pervious_fraction_m -= frozen_soil_runoff_m;

  // edited by RLM; added logic block to handle what happens when porosity or field capacity = 0
  // FLO changed from 0.95 to 1.0 because what happens to the other 5% ?
  if (max_free_water_m <= 0.0 || max_tension_water_m <= 0.0)
  {
    *flux_surface_runoff_input_to_surface_routing_m =
        *impervious_runoff_m + frozen_soil_runoff_m +
        water_input_pervious_fraction_m;
    *infiltration_depth_m = 0.0;
    return;
  }

  // if code gets to here, then field capacity and soil porosity were both nonzero, so there is some space in
  // the soil ffor more water, AND impervious runoff has already been abstracted from precipitation as impervious_runoff_m.

  // solve pervious surface runoff (m) based on Eq. 310
  // Use notation from Knoben et al. (2019) p. 71

  double a = parms->a_Xinanjiang_inflection_point_parameter; // tension water contributing area curve inflection point (called c in Jayawardena & Zhue)
  double b = parms->b_Xinanjiang_shape_parameter;            // tension water contributing area curve shape parameter (also b in Jayawardena & Zhue)
  double Ex = parms->x_Xinanjiang_shape_parameter;           // free water contributing area curve shape parameter ffor direct runoff
  double W = tension_water_m;
  double Wmax = max_tension_water_m;

  if ((W / Wmax) <= (0.5 - a))
  {
    // THIS LINE IN ORIGINAL CFE XINANJIANG CODE IS A BUG
    //  f_over_F = (pow((0.5 - a),(1.0 - b)) * pow((1.0 - (tension_water_m/max_tension_water_m)) , b));

    // FLO correct comparing against Eqn. 2a in Jayawardena & Zhou (2000)
    f_over_F = (pow((0.5 - a), (1.0 - b)) * pow(W / Wmax, b)); // Eqn. 2a in Jayawardena & Zhue
  }
  else
  {                                                                        // if ( (0.5-a) < W / Wmax )
    f_over_F = 1.0 - pow((0.5 + a), (1.0 - b)) * pow((1.0 - W / Wmax), b); // Eqn. 2b in Jayawardena & Zhue
  }

  double R = water_input_pervious_fraction_m * f_over_F; // water moved from tension to free water storage

  double S = free_water_m;
  double Smax = max_free_water_m;

  xinanjiang_runoff_m = R * (1.0 - pow((1.0 - (S / Smax)), Ex));
  *flux_surface_runoff_input_to_surface_routing_m =
      *impervious_runoff_m + frozen_soil_runoff_m + xinanjiang_runoff_m;

  // Infiltration is the portion of the remaining pervious-area input that was
  // not converted to frozen-soil or Xinanjiang runoff.
  *infiltration_depth_m =
      water_input_pervious_fraction_m - xinanjiang_runoff_m;

#ifdef DEBUG
  if (fabs(water_input_depth_m - (*infiltration_depth_m) - (*flux_surface_runoff_input_to_surface_routing_m)) > 1.0e-06)
    printf("volball err. warning: %f\n",
           water_input_depth_m - (*infiltration_depth_m) - (*flux_surface_runoff_input_to_surface_routing_m));
#endif
  return;
}

// ##############################################################
// ####################   ET FROM RAINFALL   ####################
// ##############################################################
void et_from_rainfall(double *timestep_rainfall_input_m, struct EVAPOTRANSPIRATION_STRUCTURE *et_struct)
{
  /*
      iff it is raining, take PET from rainfall first.  Wet veg. is efficient evaporator.
  */

  if (*timestep_rainfall_input_m > 0.0)
  {

    if (*timestep_rainfall_input_m > et_struct->potential_et_m_per_timestep)
    {

      et_struct->actual_et_from_rain_m_per_timestep = et_struct->potential_et_m_per_timestep;
      *timestep_rainfall_input_m -= et_struct->actual_et_from_rain_m_per_timestep;
    }

    else
    {
      // LKC: This was incorrectly set to potential instead of actual
      et_struct->actual_et_from_rain_m_per_timestep = *timestep_rainfall_input_m;
      *timestep_rainfall_input_m = 0.0;
    }
    // Move this out of the loop since EVPT needs to be corrected wheter R > EVPT or not
    et_struct->reduced_potential_et_m_per_timestep = et_struct->potential_et_m_per_timestep - et_struct->actual_et_from_rain_m_per_timestep;
  }
}

void et_from_soil(struct CONCEPTUAL_RESERVOIR_STRUCTURE *soil_res,
                  struct EVAPOTRANSPIRATION_STRUCTURE *et_struct,
                  struct NWM_SOIL_PARAMETERS_STRUCTURE *soil_parms)
{
  /*
      take AET from soil moisture storage,
      using Budyko type function to limit PET if wilting<soilmoist<field_capacity
  */
  double Budyko_numerator;
  double Budyko_denominator;
  double Budyko;

  /*-------------------- Root zone adjusted AET development -rlm --------------------------*/

  et_struct->actual_et_from_soil_m_per_timestep = 0;

  // if rootzone-based AET turned ON
  if (soil_res->is_aet_rootzone)
  {
    // Assuming the disc with the most moisture is the bottom disc of the root zone (max_rootzone_disc)
    // Convert volumetric soil moisture from the max root zone disc to moisture content [m] (disc_storage_m)
    double disc_storage_m = soil_res->smc_profile[soil_res->max_rootzone_disc] *
                            soil_res->delta_soil_disc_depth_m[soil_res->max_rootzone_disc];

    // If the moisture content from the disc with the most moisture is less than the
    // wilting point, actual et from this timestep is 0 and no water is removed.
    if (soil_res->smc_profile[soil_res->max_rootzone_disc] <= soil_parms->wltsmc)
    {
      et_struct->actual_et_from_soil_m_per_timestep = 0;
    }
    // calculate the amount of moisture removed by evapotranspiration for the bottom disc of the root zone
    else if (soil_res->smc_profile[soil_res->max_rootzone_disc] >= soil_res->soil_water_content_field_capacity)
    {
      et_struct->actual_et_from_soil_m_per_timestep = min(et_struct->reduced_potential_et_m_per_timestep, disc_storage_m);
    }
    else
    {
      Budyko_numerator = soil_res->smc_profile[soil_res->max_rootzone_disc] - soil_parms->wltsmc;
      Budyko_denominator = soil_res->soil_water_content_field_capacity - soil_parms->wltsmc;
      Budyko = Budyko_numerator / Budyko_denominator;

      et_struct->actual_et_from_soil_m_per_timestep = min(Budyko * et_struct->reduced_potential_et_m_per_timestep, disc_storage_m);
    }

    // Reduce remaining PET and remove moisture from soil profile equal to the calculated AET (actual_et_from_soil_m_per_timestep
    et_struct->reduced_potential_et_m_per_timestep -= et_struct->actual_et_from_soil_m_per_timestep;
    soil_res->smc_profile[soil_res->max_rootzone_disc] -= (et_struct->actual_et_from_soil_m_per_timestep /
                                                           soil_res->delta_soil_disc_depth_m[soil_res->max_rootzone_disc]);
    soil_res->storage_m -= et_struct->actual_et_from_soil_m_per_timestep;
  }

  else if (et_struct->reduced_potential_et_m_per_timestep > 0)
  {

    if (soil_res->storage_m >= soil_res->storage_threshold_primary_m)
    {
      et_struct->actual_et_from_soil_m_per_timestep = min(et_struct->reduced_potential_et_m_per_timestep, soil_res->storage_m);
    }
    else if (soil_res->storage_m > soil_parms->wilting_point_m && soil_res->storage_m < soil_res->storage_threshold_primary_m)
    {

      Budyko_numerator = soil_res->storage_m - soil_parms->wilting_point_m;
      Budyko_denominator = soil_res->storage_threshold_primary_m - soil_parms->wilting_point_m;
      Budyko = Budyko_numerator / Budyko_denominator;
      // LKC: Include check to guarantee EAT is not larger than soil storage
      et_struct->actual_et_from_soil_m_per_timestep = min(Budyko * et_struct->reduced_potential_et_m_per_timestep, soil_res->storage_m);
    }
    soil_res->storage_m -= et_struct->actual_et_from_soil_m_per_timestep;
    et_struct->reduced_potential_et_m_per_timestep = et_struct->reduced_potential_et_m_per_timestep - et_struct->actual_et_from_soil_m_per_timestep;
  }
}
