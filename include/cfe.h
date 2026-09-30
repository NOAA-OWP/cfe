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
#ifndef CFE_CFE_H
#define CFE_CFE_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "cfe_types.h"
#include "cfe_helpers.h"
#include "parser_helpers.h"
#include "cfe.h"
#include "conceptual_reservoir.h"
#include "giuh.h"
#include "nash_cascade.h"
#include "version.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif 
//define DEBUG

#define MAX_NUM_RAIN_DATA 720

#define WATER_LIQUID_DENSITY_kg_per_m3 998.0
#define GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2 9.80655
#define STANDARD_ATM_PRESS_Pa 101325.0

// t-shirt approximation of the hydrologic routing funtionality of the National Water Model V. 3.1 and earlier.
// This code was developed to test the hypothesis that the National Water Model runoff generation, vadose zone
// dynamics, and conceptual groundwater model can be greatly simplified by acknowledging that it is truly a
// conceptual model.  The hypothesis is supported by a number of observations made during a 2017-2018 deep dive
// into the NWM code.  Thesed are:
//
// 1. Rainfall/throughfall/melt partitioning in the NWM is based on a simple curve-number like approach that
//    was developed by Schaake et al. (1996) and which is very similar to the Probability Distributed Moisture (PDM)
//    function by Moore, 1985.   The Schaake function is a single valued function of soil moisture deficit,
//    predicts 100% runoff when the soil is saturated, like the curve-number method, and is fundamentally simple.
// 2. Run-on infiltration is strictly not calculated.  In the WRF-Hydro based NWM, the overland flow routing 
//    scheme applies the Schaake function repeatedly to predict this phenomenon, which violates the 
//    underlying assumption of the PDM method that only rainfall inputs affect soil moisture.  In reality, runon
//    infiltration happens in rill/gully flow, not over entire areas (grids in the NWM situation).
// 3. The water-content based Richards' equation, applied using a coarse-discretization, can be replaced with a simple
//    conceptual reservoir because it never allows saturation or infiltration-excess runoff unless deactivated by
//    assuming no-flow lower boundary condition.  Since this form of Richards' equation cannot simulate heterogeneous
//    soil discretizations, it can be replaced with a conceptual reservoir.
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
//       is identical at both outlets, and is based on an integral solution of the storage in the soil assuming
//       Clapp-Hornberger parameters equal to those used in the NWM to determine that storage corresponding to a
//       soil water content 0.5 m above the soil column bottom that produces a soil suction head equal to -1/3 atm,
//       which is a commonly applied assumption used to estimate the field capacity water content.
//       The first outlet calculates vertical percolation of water to deep groundwater using the saturated hydraulic
//       conductivity of the soil multiplied by the NWM "slope" parameter, which when 1.0 indicates free drainage and
//       when 0.0 indicates a no-flow lower boundary condition.   The second outlet is used to calculate the flux to
//       the soil lateral flow path, using a conceptual LKSATFAC-like calibration parameter.
//    4. The lateral flow is routed to the catchment outlet using a Nash-cascade of reservoirs to produce a 
//       volume conserving delayed response, and elminates the need to apply 250 m lateral flow routing on a grid.
//    5. The groundwater contribution to base flow is modeled using either (a) an exponential nonlinear reservoir
//       identical to the one in the NWM formulation, or (b) a nonlinear reservoir forumulation, which can also be
//       made linear by assuming an exponent value equal to 1.0.
//
//    This model was conceived and first written by Fred L. Ogden, May, 2020, in the service of the 
//    NOAA-NWS Office of Water Prediction, in Tuscaloosa, Alabama.
//
//  Version history:
//  2.1  July, 2025, Modified to read/parse keyword based input config file, prototyped model definition file.
//  3.0  August, 2025, Mofified to include discretized stateless soil moisture balance model that solves
//       soil moisture at Noah-MP discretizations (top-down) of 0.1, 0.3, 0.6, and 1.0 m
//########################################################################################################

// define data structures
//--------------------------

struct NWM_SOIL_PARAMETERS_STRUCTURE {
    // using same nondescriptive variable names as used in NWM.  <sorry>
    double smcmax;  // effective porosity [V/V]
    double wltsmc;  // wilting point soil moisture content [V/V]
    double satdk;   // saturated hydraulic conductivity [m s-1]
    double satpsi;	// saturated capillary head [m]
    double bb;      // beta exponent on Clapp-Hornberger (1978) soil water relations [-]
    // NOT USED: double mult;    // the multiplier applied to satdk to route water rapidly downslope
    double slop;   // this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow.
    double D;       // soil depth [m]
    double wilting_point_m;
    double alpha_fc;
    double refkdt;
    double soil_storage;
    };

struct EVAPOTRANSPIRATION_STRUCTURE {
    double potential_et_m_per_s;
    double potential_et_m_per_timestep;
    double reduced_potential_et_m_per_timestep;
    double actual_et_from_rain_m_per_timestep;
    double actual_et_from_soil_m_per_timestep;
    double actual_bare_soil_evaporation_m_per_timestep;
    double actual_et_m_per_timestep;
};
typedef struct EVAPOTRANSPIRATION_STRUCTURE evapotranspiration_structure;


/* Xinanjiang*/
struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE {
    cfe_partition_scheme_t surface_water_partitioning_scheme;
    double Schaake_adjusted_magic_constant_by_soil_type;
    double a_Xinanjiang_inflection_point_parameter;
    double b_Xinanjiang_shape_parameter;
    double x_Xinanjiang_shape_parameter;
    double urban_decimal_fraction;
    double ice_content_threshold; // ice content above which soil is impermeable
};
typedef struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE rainfall_partitioning_parameters_structure;


// function prototypes
// --------------------------------
void Schaake_partitioning_scheme
      (
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
      double ice_fraction_schaake, double ice_content_threshold
      );

// xinanjiang_dev: XinJiang function written by Rachel adapted by Jmframe with review and bug fixes by FLO producing v. 2.1,
extern void Xinanjiang_partitioning_scheme
        (
        double water_input_depth_m,
        double field_capacity_m,
        double max_soil_moisture_storage_m,
        double column_total_soil_water_m,
        struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE *parms,
        double *flux_surface_runoff_input_to_surface_routing_m,
        double *infiltration_depth_m,
        double *impervious_runoff_m,
        double ice_fraction_xinanjiang
        );

extern void et_from_rainfall(double *timestep_rainfall_input_m, evapotranspiration_structure *et_struct);


extern void et_from_soil
    (
    struct CONCEPTUAL_RESERVOIR_STRUCTURE *soil_res, 
    struct EVAPOTRANSPIRATION_STRUCTURE *et_struct,
    struct NWM_SOIL_PARAMETERS_STRUCTURE *soil_parms
    );



extern int is_fabs_less_than_epsilon(double a,double epsilon);

extern void cfe(
        double *soil_reservoir_storage_deficit_m_ptr,
        struct NWM_SOIL_PARAMETERS_STRUCTURE NWM_soil_params_struct,
        struct CONCEPTUAL_RESERVOIR_STRUCTURE *soil_reservoir_struct,
        double timestep_h,
        struct RAINFALL_PARTITIONING_PARAMETERS_STRUCTURE infiltration_excess_params_struct,
        double timestep_rainfall_input_m,
        double *infiltration_excess_m_ptr,
        double *infiltration_depth_m_ptr,
        double *flux_perc_m_ptr,
        double *flux_lat_m_ptr,
        double *gw_reservoir_storage_deficit_m_ptr,
        struct CONCEPTUAL_RESERVOIR_STRUCTURE *gw_reservoir_struct,
        double *flux_from_deep_gw_to_chan_m_ptr,
        double *giuh_runoff_m_ptr,
        int num_giuh_ordinates,
        double *giuh_ordinates_arr,
        double *giuh_runoff_queue_m_per_timestep_arr,
        double *nash_lateral_runoff_m_ptr,
        struct NASH_CASCADE_PARAMETERS_STRUCTURE *nash_subsurface_params,
        struct EVAPOTRANSPIRATION_STRUCTURE *evap_struct,
        double *Qout_m_ptr,
        cfe_volbal_struct *volbal_struct,
        double time_step_seconds,
        int surface_runoff_scheme,
        int yes_simulate_discrete_soil_moisture,
        SoilControl* soil_control,
        SoilGeometry* soil_geometry,
        SoilParameters* soil_parameters,
        SoilStateIn* soil_state_in,
        SoilStateOut* soil_state_out,
        SoilFluxes* soil_fluxes,
        SoilLookupTables* ch_lookup_tables,
        double* soil_discrete_storage_theta,
        int verbosity      
    );

#endif //CFE_CFE_H
