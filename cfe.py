import time
import numpy as np
import pandas as pd
import sys

class CFE():
    def __init__(self):
        super(CFE, self).__init__()
        
    # __________________________________________________________________________________________________________
    # MAIN MODEL FUNCTION
    def run_cfe(self, cfe_state):
        
        # ________________________________________________
        cfe_state.volin += cfe_state.timestep_rainfall_input_m
        
        # ________________________________________________
        cfe_state.potential_et_m_per_timestep = cfe_state.potential_et_m_per_s * cfe_state.time_step_size
        cfe_state.reduced_potential_et_m_per_timestep = cfe_state.potential_et_m_per_s * cfe_state.time_step_size

        # ________________________________________________
        # SUBROUTINE
        # timestep_rainfall_input_m = f(timestep_rainfall_input_m, potential_et_m_per_timestep)
        cfe_state.actual_et_from_rain_m_per_timestep = 0
        if(cfe_state.timestep_rainfall_input_m > 0):
            self.et_from_rainfall(cfe_state)
        
        # ________________________________________________
        cfe_state.vol_et_from_rain += cfe_state.actual_et_from_rain_m_per_timestep
        #cfe_state.vol_et_to_atm += cfe_state.actual_et_from_rain_m_per_timestep
        cfe_state.volout += cfe_state.actual_et_from_rain_m_per_timestep

        # ---------------------- SUBROUTINE ---------------------- #
        # ET from soil
        cfe_state.actual_et_from_soil_m_per_timestep = 0
        if(cfe_state.soil_reservoir['storage_m'] > cfe_state.soil_reservoir['wilting_point_m']): 
            self.et_from_soil(cfe_state)

        cfe_state.vol_et_from_soil += cfe_state.actual_et_from_soil_m_per_timestep
        #cfe_state.vol_et_to_atm += cfe_state.actual_et_from_soil_m_per_timestep;
        cfe_state.volout += cfe_state.actual_et_from_soil_m_per_timestep 

        cfe_state.actual_et_m_per_timestep= cfe_state.actual_et_from_rain_m_per_timestep + cfe_state.actual_et_from_soil_m_per_timestep
  
        # ________________________________________________        
        cfe_state.soil_reservoir_storage_deficit_m = (cfe_state.soil_params['smcmax'] * \
                                                 cfe_state.soil_params['D'] - \
                                                 cfe_state.soil_reservoir['storage_m'])
        
        # ________________________________________________
        # SUBROUTINE
        # Calculates the value for surface_runoff_depth_m
        if (cfe_state.timestep_rainfall_input_m > 0.0): 
            if cfe_state.surface_partitioning_scheme == "Schaake": 
                self.Schaake_partitioning_scheme(cfe_state)
            elif cfe_state.surface_partitioning_scheme == "Xinanjiang": 
                self.Xinanjiang_partitioning_scheme(cfe_state)
            else: 
                print("Problem: must specify one of Schaake of Xinanjiang partitioning scheme.\n")
                print("Program terminating.:( \n");
                sys.exit(1)
        else: 
            cfe_state.surface_runoff_depth_m = 0.0
            cfe_state.infiltration_depth_m = 0.0

        # ________________________________________________
        if cfe_state.soil_reservoir_storage_deficit_m < cfe_state.infiltration_depth_m:
            # put won't fit back into runoff
            cfe_state.surface_runoff_depth_m += (cfe_state.infiltration_depth_m - cfe_state.soil_reservoir_storage_deficit_m)
            cfe_state.infiltration_depth_m = cfe_state.soil_reservoir_storage_deficit_m
            cfe_state.soil_reservoir['storage_m'] = cfe_state.soil_reservoir['storage_max_m']
            cfe_state.soil_reservoir_storage_deficit_m = 0

        # ________________________________________________
        cfe_state.vol_sch_runoff += cfe_state.surface_runoff_depth_m
        cfe_state.vol_sch_infilt += cfe_state.infiltration_depth_m

        # ________________________________________________
        if cfe_state.current_time_step == 0:
            cfe_state.previous_flux_perc_m = cfe_state.flux_perc_m
            
        # ________________________________________________
        if cfe_state.previous_flux_perc_m > cfe_state.soil_reservoir_storage_deficit_m:
            diff = cfe_state.previous_flux_perc_m - cfe_state.soil_reservoir_storage_deficit_m
            cfe_state.infiltration_depth_m = cfe_state.soil_reservoir_storage_deficit_m
            cfe_state.vol_sch_runoff += diff
            cfe_state.vol_sch_infilt -= diff
            cfe_state.surface_runoff_depth_m += diff
            
            cfe_state.soil_reservoir_storage_deficit_m = 0
        # ________________________________________________
        cfe_state.vol_to_soil += cfe_state.infiltration_depth_m
        cfe_state.soil_reservoir['storage_m'] += cfe_state.infiltration_depth_m

        # ________________________________________________
        # SUBROUTINE
        # primary_flux, secondary_flux = f(reservoir)
        self.conceptual_reservoir_flux_calc(cfe_state, cfe_state.soil_reservoir)

        # ________________________________________________
        cfe_state.flux_perc_m = cfe_state.primary_flux_m
        cfe_state.flux_lat_m = cfe_state.secondary_flux_m

        # ________________________________________________
        cfe_state.gw_reservoir_storage_deficit_m = cfe_state.gw_reservoir['storage_max_m'] - cfe_state.gw_reservoir['storage_m']
        
        # ________________________________________________
        if cfe_state.flux_perc_m > cfe_state.gw_reservoir_storage_deficit_m:
            diff = cfe_state.flux_perc_m - cfe_state.gw_reservoir_storage_deficit_m
            cfe_state.flux_perc_m = cfe_state.gw_reservoir_storage_deficit_m
            cfe_state.vol_sch_runoff+=diff 
            cfe_state.vol_sch_infilt-=diff 
            
        # ________________________________________________
        cfe_state.vol_to_gw                += cfe_state.flux_perc_m
        cfe_state.vol_soil_to_gw           += cfe_state.flux_perc_m

        cfe_state.gw_reservoir['storage_m']   += cfe_state.flux_perc_m
        cfe_state.soil_reservoir['storage_m'] -= cfe_state.flux_perc_m
        cfe_state.soil_reservoir['storage_m'] -= cfe_state.flux_lat_m
        cfe_state.vol_soil_to_lat_flow        += cfe_state.flux_lat_m  #TODO add this to nash cascade as input
        cfe_state.volout                       = cfe_state.volout + cfe_state.flux_lat_m

            
        # ________________________________________________
        # SUBROUTINE
        # primary_flux, secondary_flux = f(reservoir)
        self.conceptual_reservoir_flux_calc(cfe_state, cfe_state.gw_reservoir) 
            
        # ________________________________________________
        cfe_state.flux_from_deep_gw_to_chan_m = cfe_state.primary_flux_m
        if (cfe_state.flux_from_deep_gw_to_chan_m > cfe_state.gw_reservoir['storage_m']): 
            cfe_state.flux_from_deep_gw_to_chan_m = cfe_state.gw_reservoir['storage_m']
            print("WARNING: Groundwater flux larger than storage. \n")
        
        cfe_state.vol_from_gw += cfe_state.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________        
        if not self.is_fabs_less_than_epsilon(cfe_state.secondary_flux, 1.0e-09):
            print("problem with nonzero flux point 1\n")
                        
        # ________________________________________________                               
        cfe_state.gw_reservoir['storage_m'] -= cfe_state.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________
        # SUBROUTINE
        # giuh_runoff_m = f(Schaake_output, giuh_ordinates, runoff_queue_m_per_timestep)
        self.convolution_integral(cfe_state)
        
        # ________________________________________________
        cfe_state.vol_out_giuh += cfe_state.flux_giuh_runoff_m
        cfe_state.volout += cfe_state.flux_giuh_runoff_m + cfe_state.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________
        # SUBROUTINE
        self.nash_cascade(cfe_state)

        # ________________________________________________
        cfe_state.vol_in_nash += cfe_state.flux_lat_m
        cfe_state.vol_out_nash += cfe_state.flux_nash_lateral_runoff_m
        
        # ________________________________________________
        cfe_state.flux_Qout_m = cfe_state.flux_giuh_runoff_m + cfe_state.flux_nash_lateral_runoff_m + cfe_state.flux_from_deep_gw_to_chan_m
        cfe_state.total_discharge = cfe_state.flux_Qout_m * cfe_state.catchment_area_km2 * 1000000.0 / 3600.0
        
        # ________________________________________________
        cfe_state.current_time_step += 1
        cfe_state.current_time      += pd.Timedelta(value=cfe_state.time_step_size, unit='s')

        return
    
    
    # __________________________________________________________________________________________________________
    def nash_cascade(self,cfe_state):
        """
            Solve for the flow through the Nash cascade to delay the 
            arrival of the lateral flow into the channel
        """
        Q = np.zeros(cfe_state.num_lateral_flow_nash_reservoirs)
        
        for i in range(cfe_state.num_lateral_flow_nash_reservoirs):
            
            Q[i] = cfe_state.K_nash * cfe_state.nash_storage[i]
            
            cfe_state.nash_storage[i] -= Q[i]
            
            if i == 0:
                
                cfe_state.nash_storage[i] += cfe_state.flux_lat_m
                
            else:
                
                cfe_state.nash_storage[i] += Q[i-1]
        
        cfe_state.flux_nash_lateral_runoff_m = Q[cfe_state.num_lateral_flow_nash_reservoirs - 1]
        
        return
    
                               
    # __________________________________________________________________________________________________________
    def convolution_integral(self,cfe_state):
        """
            This function solves the convolution integral involving N GIUH ordinates.
            
            Inputs:
                Schaake_output_runoff_m
                num_giuh_ordinates
                giuh_ordinates
            Outputs:
                runoff_queue_m_per_timestep
        """

#        cfe_state.runoff_queue_m_per_timestep[-1] = 0
        
        for i in range(cfe_state.num_giuh_ordinates): 

            cfe_state.runoff_queue_m_per_timestep[i] += cfe_state.giuh_ordinates[i] * cfe_state.surface_runoff_depth_m
            
        cfe_state.flux_giuh_runoff_m = cfe_state.runoff_queue_m_per_timestep[0]
        
        # __________________________________________________________________
        # shift all the entries in preperation for the next timestep
        
        for i in range(1, cfe_state.num_giuh_ordinates):  
            
            cfe_state.runoff_queue_m_per_timestep[i-1] = cfe_state.runoff_queue_m_per_timestep[i]

        cfe_state.runoff_queue_m_per_timestep[-1] = 0

        return
    

    # __________________________________________________________________________________________________________
    def et_from_rainfall(self,cfe_state):
        
        """
            iff it is raining, take PET from rainfall first.  Wet veg. is efficient evaporator.
        """
        
        if cfe_state.timestep_rainfall_input_m >0.0:

            if cfe_state.timestep_rainfall_input_m > cfe_state.potential_et_m_per_timestep:
        
                cfe_state.actual_et_from_rain_m_per_timestep = cfe_state.potential_et_m_per_timestep
                cfe_state.timestep_rainfall_input_m -= cfe_state.actual_et_from_rain_m_per_timestep

            else: 

                cfe_state.actual_et_from_rain_m_per_timestep = cfe_state.timestep_rainfall_input_m
                cfe_state.timestep_rainfall_input_m=0.0
        
        cfe_state.reduced_potential_et_m_per_timestep = cfe_state.potential_et_m_per_timestep-cfe_state.actual_et_from_rain_m_per_timestep

        return
                
                
    # __________________________________________________________________________________________________________
    ########## SINGLE OUTLET EXPONENTIAL RESERVOIR ###############
    ##########                -or-                 ###############
    ##########    TWO OUTLET NONLINEAR RESERVOIR   ###############                        
    def conceptual_reservoir_flux_calc(self,cfe_state,reservoir):
        """
            This function calculates the flux from a linear, or nonlinear 
            conceptual reservoir with one or two outlets, or from an
            exponential nonlinear conceptual reservoir with only one outlet.
            In the non-exponential instance, each outlet can have its own
            activation storage threshold.  Flow from the second outlet is 
            turned off by setting the discharge coeff. to 0.0.
        """

        if reservoir['is_exponential'] == True: 
            flux_exponential = np.exp(reservoir['exponent_primary'] * \
                                      reservoir['storage_m'] / \
                                      reservoir['storage_max_m']) - 1.0
            cfe_state.primary_flux_m = reservoir['coeff_primary'] * flux_exponential
            cfe_state.secondary_flux_m=0.0
            return
    
        cfe_state.primary_flux_m=0.0
        
        storage_above_threshold_m = reservoir['storage_m'] - reservoir['storage_threshold_primary_m']
        
        if storage_above_threshold_m > 0.0:
                               
            storage_diff = reservoir['storage_max_m'] - reservoir['storage_threshold_primary_m']
            storage_ratio = storage_above_threshold_m / storage_diff
            storage_power = np.power(storage_ratio, reservoir['exponent_primary'])
            
            cfe_state.primary_flux_m = reservoir['coeff_primary'] * storage_power

            if cfe_state.primary_flux_m > storage_above_threshold_m:
                cfe_state.primary_flux_m = storage_above_threshold_m
                
        cfe_state.secondary_flux_m = 0.0
            
        storage_above_threshold_m = reservoir['storage_m'] - reservoir['storage_threshold_secondary_m']
        
        if storage_above_threshold_m > 0.0:
            
            storage_diff = reservoir['storage_max_m'] - reservoir['storage_threshold_secondary_m']
            storage_ratio = storage_above_threshold_m / storage_diff
            storage_power = np.power(storage_ratio, reservoir['exponent_secondary'])
            
            cfe_state.secondary_flux_m = reservoir['coeff_secondary'] * storage_power
            
            if cfe_state.secondary_flux_m > (storage_above_threshold_m - cfe_state.primary_flux_m):
                cfe_state.secondary_flux_m = storage_above_threshold_m - cfe_state.primary_flux_m
                
        return
    
    
    # __________________________________________________________________________________________________________
    #  SCHAAKE RUNOFF PARTITIONING SCHEME
    def Schaake_partitioning_scheme(self,cfe_state):
        """
            This subtroutine takes water_input_depth_m and partitions it into surface_runoff_depth_m and
            infiltration_depth_m using the scheme from Schaake et al. 1996. 
            !--------------------------------------------------------------------------------
            modified by FLO April 2020 to eliminate reference to ice processes, 
            and to de-obfuscate and use descriptive and dimensionally consistent variable names.
            
            inputs:
              timestep_d
              Schaake_adjusted_magic_constant_by_soil_type = C*Ks(soiltype)/Ks_ref, where C=3, and Ks_ref=2.0E-06 m/s
              column_total_soil_moisture_deficit_m (soil_reservoir_storage_deficit_m)
              water_input_depth_m (timestep_rainfall_input_m) amount of water input to soil surface this time step [m]
            outputs:
              surface_runoff_depth_m      amount of water partitioned to surface water this time step [m]
              infiltration_depth_m
        """
        
        if 0 < cfe_state.timestep_rainfall_input_m:
            
            if 0 > cfe_state.soil_reservoir_storage_deficit_m:
                
                cfe_state.surface_runoff_depth_m = cfe_state.timestep_rainfall_input_m
                
                cfe_state.infiltration_depth_m = 0.0
                
            else:
                
                schaake_exp_term = np.exp( - cfe_state.Schaake_adjusted_magic_constant_by_soil_type * cfe_state.timestep_d)
                
                Schaake_parenthetical_term = (1.0 - schaake_exp_term)
                
                Ic = cfe_state.soil_reservoir_storage_deficit_m * Schaake_parenthetical_term
                
                Px = cfe_state.timestep_rainfall_input_m
                
                cfe_state.infiltration_depth_m = (Px * (Ic / (Px + Ic)))
                
                if 0.0 < (cfe_state.timestep_rainfall_input_m - cfe_state.infiltration_depth_m):
                    
                    cfe_state.surface_runoff_depth_m = cfe_state.timestep_rainfall_input_m - cfe_state.infiltration_depth_m
                    
                else:
                    
                    cfe_state.surface_runoff_depth_m = 0.0
                    
                    cfe_state.infiltration_depth_m =  cfe_state.timestep_rainfall_input_m - cfe_state.surface_runoff_depth_m
                    
        else:
            
            cfe_state.surface_runoff_depth_m = 0.0
            
            cfe_state.infiltration_depth_m = 0.0
            
        return

    # __________________________________________________________________________________________________________
    def Xinanjiang_partitioning_scheme(self,cfe_state): 
        """
            This module takes the water_input_depth_m and separates it into surface_runoff_depth_m
            and infiltration_depth_m by calculating the saturated area and runoff based on a scheme developed
            for the Xinanjiang model by Jaywardena and Zhou (2000). According to Knoben et al.
            (2019) "the model uses a variable contributing area to simulate runoff.  [It] uses
            a double parabolic curve to simulate tension water capacities within the catchment, 
            instead of the original single parabolic curve" which is also used as the standard 
            VIC fomulation.  This runoff scheme was selected for implementation into NWM v3.0.
            REFERENCES:
            1. Jaywardena, A.W. and M.C. Zhou, 2000. A modified spatial soil moisture storage 
                capacity distribution curve for the Xinanjiang model. Journal of Hydrology 227: 93-113
            2. Knoben, W.J.M. et al., 2019. Supplement of Modular Assessment of Rainfall-Runoff Models
                Toolbox (MARRMoT) v1.2: an open-source, extendable framework providing implementations
                of 46 conceptual hydrologic models as continuous state-space formulations. Supplement of 
                Geosci. Model Dev. 12: 2463-2480.
            -------------------------------------------------------------------------
            Written by RLM May 2021
            Adapted by JMFrame September 2021 for new version of CFE
            Further adapted by QiyueL August 2022 for python version of CFE
            ------------------------------------------------------------------------
            Inputs
            double  time_step_rainfall_input_m           amount of water input to soil surface this time step [m]
            double  field_capacity_m                     amount of water stored in soil reservoir when at field capacity [m]
            double  max_soil_moisture_storage_m          total storage of the soil moisture reservoir (porosity*soil thickness) [m]
            double  column_total_soil_water_m     current storage of the soil moisture reservoir [m]
            double  a_inflection_point_parameter  a parameter
            double  b_shape_parameter             b parameter
            double  x_shape_parameter             x parameter
                //
            Outputs
            double  surface_runoff_depth_m        amount of water partitioned to surface water this time step [m]
            double  infiltration_depth_m          amount of water partitioned as infiltration (soil water input) this time step [m]
            ------------------------------------------------------------------------- 
        """

        # partition the total soil water in the column between free water and tension water
        free_water_m = cfe_state.soil_reservoir['storage_m']- cfe_state.soil_reservoir['storage_threshold_primary_m'];

        if (0.0 < free_water_m):

            tension_water_m = cfe_state.soil_reservoir['storage_threshold_primary_m']

        else: 

            free_water_m = 0.0
            tension_water_m = cfe_state.soil_reservoir['storage_m']
        
        # estimate the maximum free water and tension water available in the soil column
        max_free_water_m = cfe_state.soil_reservoir['storage_max_m'] - cfe_state.soil_reservoir['storage_threshold_primary_m']
        max_tension_water_m = cfe_state.soil_reservoir['storage_threshold_primary_m']

        # check that the free_water_m and tension_water_m do not exceed the maximum and if so, change to the max value
        if(max_free_water_m < free_water_m): 
            free_water_m = max_free_water_m

        if(max_tension_water_m < tension_water_m): 
            tension_water_m = max_tension_water_m

        """
            NOTE: the impervious surface runoff assumptions due to frozen soil used in NWM 3.0 have not been included.
            We are assuming an impervious area due to frozen soils equal to 0 (see eq. 309 from Knoben et al).

            The total (pervious) runoff is first estimated before partitioning into surface and subsurface components.
            See Knoben et al eq 310 for total runoff and eqs 313-315 for partitioning between surface and subsurface
            components.

            Calculate total estimated pervious runoff. 
            NOTE: If the impervious surface runoff due to frozen soils is added,
            the pervious_runoff_m equation will need to be adjusted by the fraction of pervious area.
        """
        a_Xinanjiang_inflection_point_parameter = 1
        b_Xinanjiang_shape_parameter = 1
        x_Xinanjiang_shape_parameter = 1

        if ((tension_water_m/max_tension_water_m) <= (0.5 - a_Xinanjiang_inflection_point_parameter)): 
            pervious_runoff_m = cfe_state.timestep_rainfall_input_m * \
                (np.power((0.5 - a_Xinanjiang_inflection_point_parameter),\
                    (1.0 - b_Xinanjiang_shape_parameter)) * \
                        np.power((1.0 - (tension_water_m/max_tension_water_m)),\
                            b_Xinanjiang_shape_parameter))

        else: 
            pervious_runoff_m = cfe_state.timestep_rainfall_input_m* \
                (1.0 - np.power((0.5 + a_Xinanjiang_inflection_point_parameter), \
                    (1.0 - b_Xinanjiang_shape_parameter)) * \
                        np.power((1.0 - (tension_water_m/max_tension_water_m)),\
                            (b_Xinanjiang_shape_parameter)))
    
        # Separate the surface water from the pervious runoff 
        ## NOTE: If impervious runoff is added to this subroutine, impervious runoff should be added to
        ## the surface_runoff_depth_m.
        
        cfe_state.surface_runoff_depth_m = pervious_runoff_m * \
             (1.0 - np.power((1.0 - (free_water_m/max_free_water_m)),x_Xinanjiang_shape_parameter))

        # The surface runoff depth is bounded by a minimum of 0 and a maximum of the water input depth.
        # Check that the estimated surface runoff is not less than 0.0 and if so, change the value to 0.0.
        if(cfe_state.surface_runoff_depth_m < 0.0): 
            cfe_state.surface_runoff_depth_m = 0.0;
    
        # Check that the estimated surface runoff does not exceed the amount of water input to the soil surface.  If it does,
        # change the surface water runoff value to the water input depth.
        if(cfe_state.surface_runoff_depth_m > cfe_state.timestep_rainfall_input_m): 
             cfe_state.surface_runoff_depth_m = cfe_state.timestep_rainfall_input_m
        
        # Separate the infiltration from the total water input depth to the soil surface.
        cfe_state.infiltration_depth_m = cfe_state.timestep_rainfall_input_m - cfe_state.surface_runoff_depth_m;    

        return        
                               
    # __________________________________________________________________________________________________________
    def et_from_soil(self,cfe_state):
        """
            take AET from soil moisture storage, 
            using Budyko type curve to limit PET if wilting<soilmoist<field_capacity
        """
        
        if cfe_state.reduced_potential_et_m_per_timestep > 0:
            
            # print("this should not happen yet. Still debugging the other functions.")
            
            if cfe_state.soil_reservoir['storage_m'] >= cfe_state.soil_reservoir['storage_threshold_primary_m']:
            
                cfe_state.actual_et_from_soil_m_per_timestep = np.min(cfe_state.reduced_potential_et_m_per_timestep, 
                                                       cfe_state.soil_reservoir['storage_m'])

            elif (cfe_state.soil_reservoir['storage_m'] > cfe_state.soil_reservoir['wilting_point_m'] and 
                  cfe_state.soil_reservoir['storage_m'] < cfe_state.soil_reservoir['storage_threshold_primary_m']):
            
                Budyko_numerator = cfe_state.soil_reservoir['storage_m'] - cfe_state.soil_reservoir['wilting_point_m']
                Budyko_denominator = cfe_state.soil_reservoir['storage_threshold_primary_m'] - \
                                     cfe_state.soil_reservoir['wilting_point_m']
                Budyko = Budyko_numerator / Budyko_denominator
                               
                cfe_state.actual_et_from_soil_m_per_timestep = np.min(Budyko * cfe_state.reduced_potential_et_m_per_timestep,cfe_state.soil_reservoir['storage_m'])
                               
            cfe_state.soil_reservoir['storage_m'] -= cfe_state.actual_et_from_soil_m_per_timestep
            cfe_state.reduced_potential_et_m_per_timestep = cfe_state.reduced_potential_et_m_per_timestep - cfe_state.actual_et_from_soil_m_per_timestep
        return
            
            
    # __________________________________________________________________________________________________________
    def is_fabs_less_than_epsilon(self,a,epsilon):
        
        if np.abs(a) < epsilon:
            
            return True
        
        else:
            
            return False 
    
