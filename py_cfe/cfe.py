import time
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

class CFE():
    def __init__(self, cfe_state):
        super(CFE, self).__init__()
        
        self.timestep_rainfall_input_m = cfe_state.timestep_rainfall_input_m
gw_reservoir
soil_reservoir

        self.volin

    # __________________________________________________________________________________________________________
    # MAIN MODEL FUNCTION
    def run_cfe(self):
        
        # ________________________________________________
        self.volin += self.timestep_rainfall_input_m
        
        # ________________________________________________
        self.potential_et_m_per_timestep = self.potential_et_m_per_s * self.time_step_size
        
        # ________________________________________________
        # timestep_rainfall_input_m = f(timestep_rainfall_input_m, potential_et_m_per_timestep)
        self.et_from_rainfall()
        
        # ________________________________________________        
        self.soil_reservoir_storage_deficit_m = (self.soil_params['smcmax'] * \
                                                 self.soil_params['D'] - \
                                                 self.soil_reservoir['storage_m'])
        
        # ________________________________________________
        # Calculates the value for surface_runoff_depth_m
        self.Schaake_partitioning_scheme()

        # ________________________________________________
        self.et_from_soil()
        
        # ________________________________________________
        if self.soil_reservoir_storage_deficit_m < self.infiltration_depth_m:
            # put won't fit back into runoff
            self.surface_runoff_depth_m += (self.infiltration_depth_m - soil_reservoir_storage_deficit_m)
            self.infiltration_depth_m = self.soil_reservoir_storage_deficit_m
            self.soil_reservoir['storage_m'] = self.soil_reservoir['storage_max_m']

        # ________________________________________________
        self.vol_sch_runoff += self.surface_runoff_depth_m
        self.vol_sch_infilt += self.infiltration_depth_m

        # ________________________________________________
        if self.current_time_step == 0:
            self.previous_flux_perc_m = self.flux_perc_m
            
        # ________________________________________________
        if self.previous_flux_perc_m > self.soil_reservoir_storage_deficit_m:
            diff = self.previous_flux_perc_m - self.soil_reservoir_storage_deficit
            self.infiltration_depth_m = self.soil_reservoir_storage_deficit_m
            self.vol_sch_runoff += diff
            self.vol_sch_infilt -= diff
            self.surface_runoff_depth_m += diff
            
        # ________________________________________________
        self.vol_to_soil += self.infiltration_depth_m
        self.soil_reservoir['storage_m'] += self.infiltration_depth_m

        # ________________________________________________
        # primary_flux, secondary_flux = f(reservoir)
        self.conceptual_reservoir_flux_calc( self.soil_reservoir )

        # ________________________________________________
        self.flux_perc_m = self.primary_flux
        self.flux_lat_m = self.secondary_flux

        # ________________________________________________
        self.gw_reservoir_storage_deficit_m = self.gw_reservoir['storage_max_m'] - self.gw_reservoir['storage_m']
        
        # ________________________________________________
        if self.flux_perc_m > self.gw_reservoir_storage_deficit_m:
            diff = self.flux_perc_m - self.gw_reservoir_storage_deficit_m
            self.flux_perc_m = self.gw_reservoir_storage_deficit_m
            self.vol_sch_runoff+=diff 
            self.vol_sch_infilt-=diff 
            
        # ________________________________________________
        self.vol_to_gw                += self.flux_perc_m
        self.vol_soil_to_gw           += self.flux_perc_m

        self.gw_reservoir['storage_m']   += self.flux_perc_m
        self.soil_reservoir['storage_m'] -= self.flux_perc_m
        self.soil_reservoir['storage_m'] -= self.flux_lat_m
        self.vol_soil_to_lat_flow        += self.flux_lat_m  #TODO add this to nash cascade as input
        self.volout                       = self.volout + self.flux_lat_m;

            
        # ________________________________________________
        # primary_flux, secondary_flux = f(reservoir)
        self.conceptual_reservoir_flux_calc( self.gw_reservoir )
            
        # ________________________________________________
        self.flux_from_deep_gw_to_chan_m = self.primary_flux
        self.vol_from_gw += self.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________        
        if not self.is_fabs_less_than_epsilon(self.secondary_flux, 1.0e-09):
            print("problem with nonzero flux point 1\n")
                        
        # ________________________________________________                               
        self.gw_reservoir['storage_m'] -= self.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________
        # giuh_runoff_m = f(Schaake_output, giuh_ordinates, runoff_queue_m_per_timestep)
        self.convolution_integral()
        
        # ________________________________________________
        self.vol_out_giuh += self.flux_giuh_runoff_m
        self.volout += self.flux_giuh_runoff_m + self.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________
        self.nash_cascade()

        # ________________________________________________
        self.vol_in_nash += self.flux_lat_m
        self.vol_out_nash += self.flux_nash_lateral_runoff_m
        
        # ________________________________________________
        self.flux_Qout_m = self.flux_giuh_runoff_m + self.flux_nash_lateral_runoff_m + self.flux_from_deep_gw_to_chan_m
        self.total_discharge = self.flux_Qout_m * self.catchment_area_km2 * 1000000.0 / 3600.0
        
        # ________________________________________________
        self.current_time_step += 1
        self.current_time      += pd.Timedelta(value=self.time_step_size, unit='s')

        return
    
    
    # __________________________________________________________________________________________________________
    def nash_cascade(self):
        """
            Solve for the flow through the Nash cascade to delay the 
            arrival of the lateral flow into the channel
        """
        Q = np.zeros(self.num_lateral_flow_nash_reservoirs)
        
        for i in range(self.num_lateral_flow_nash_reservoirs):
            
            Q[i] = self.K_nash * self.nash_storage[i]
            
            self.nash_storage[i] -= Q[i]
            
            if i == 0:
                
                self.nash_storage[i] += self.flux_lat_m
                
            else:
                
                self.nash_storage[i] += Q[i-1]
        
        self.flux_nash_lateral_runoff_m = Q[self.num_lateral_flow_nash_reservoirs - 1]
        
        return
    
                               
    # __________________________________________________________________________________________________________
    def convolution_integral(self):
        """
            This function solves the convolution integral involving N GIUH ordinates.
            
            Inputs:
                Schaake_output_runoff_m
                num_giuh_ordinates
                giuh_ordinates
            Outputs:
                runoff_queue_m_per_timestep
        """

#        self.runoff_queue_m_per_timestep[-1] = 0
        
        for i in range(self.num_giuh_ordinates): 

            self.runoff_queue_m_per_timestep[i] += self.giuh_ordinates[i] * self.surface_runoff_depth_m
            
        self.flux_giuh_runoff_m = self.runoff_queue_m_per_timestep[0]
        
        # __________________________________________________________________
        # shift all the entries in preperation for the next timestep
        
        for i in range(1, self.num_giuh_ordinates):  
            
            self.runoff_queue_m_per_timestep[i-1] = self.runoff_queue_m_per_timestep[i]

        self.runoff_queue_m_per_timestep[-1] = 0

        return
    

    # __________________________________________________________________________________________________________
    def et_from_rainfall(self):
        
        """
            iff it is raining, take PET from rainfall first.  Wet veg. is efficient evaporator.
        """
        
        if self.timestep_rainfall_input_m >0.0:

            if self.timestep_rainfall_input_m > self.potential_et_m_per_timestep:
        
                self.actual_et_m_per_timestep = self.potential_et_m_per_timestep
                self.timestep_rainfall_input_m -= self.actual_et_m_per_timestep

            else: 

                self.potential_et_m_per_timestep -= self.timestep_rainfall_input_m
                self.timestep_rainfall_input_m=0.0
        return
                
                
    # __________________________________________________________________________________________________________
    ########## SINGLE OUTLET EXPONENTIAL RESERVOIR ###############
    ##########                -or-                 ###############
    ##########    TWO OUTLET NONLINEAR RESERVOIR   ###############                        
    def conceptual_reservoir_flux_calc(self, reservoir):
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
            self.primary_flux_m = reservoir['coeff_primary'] * flux_exponential
            self.secondary_flux_m=0.0
            return
    
        self.primary_flux_m=0.0
        
        storage_above_threshold_m = reservoir['storage_m'] - reservoir['storage_threshold_primary_m']
        
        if storage_above_threshold_m > 0.0:
                               
            storage_diff = reservoir['storage_max_m'] - reservoir['storage_threshold_primary_m']
            storage_ratio = storage_above_threshold_m / storage_diff
            storage_power = np.power(storage_ratio, reservoir['exponent_primary'])
            
            self.primary_flux_m = reservoir['coeff_primary'] * storage_power

            if self.primary_flux_m > storage_above_threshold_m:
                self.primary_flux_m = storage_above_threshold_m
                
        self.secondary_flux_m = 0.0;
            
        storage_above_threshold_m = reservoir['storage_m'] - reservoir['storage_threshold_secondary_m']
        
        if storage_above_threshold_m > 0.0:
            
            storage_diff = reservoir['storage_max_m'] - reservoir['storage_threshold_secondary_m']
            storage_ratio = storage_above_threshold_m / storage_diff
            storage_power = np.power(storage_ratio, reservoir['exponent_secondary'])
            
            self.secondary_flux_m = reservoir['coeff_secondary'] * storage_power
            
            if self.secondary_flux_m > (storage_above_threshold_m - self.primary_flux_m):
                self.secondary_flux_m = storage_above_threshold_m - self.primary_flux_m
                
        return
    
    
    # __________________________________________________________________________________________________________
    #  SCHAAKE RUNOFF PARTITIONING SCHEME
    def Schaake_partitioning_scheme(self):
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
        
        if 0 < self.timestep_rainfall_input_m:
            
            if 0 > self.soil_reservoir_storage_deficit_m:
                
                self.surface_runoff_depth_m = self.timestep_rainfall_input_m
                
                self.infiltration_depth_m = 0.0
                
            else:
                
                schaake_exp_term = np.exp( - self.Schaake_adjusted_magic_constant_by_soil_type * self.timestep_d)
                
                Schaake_parenthetical_term = (1.0 - schaake_exp_term)
                
                Ic = self.soil_reservoir_storage_deficit_m * Schaake_parenthetical_term
                
                Px = self.timestep_rainfall_input_m
                
                self.infiltration_depth_m = (Px * (Ic / (Px + Ic)))
                
                if 0.0 < (self.timestep_rainfall_input_m - self.infiltration_depth_m):
                    
                    self.surface_runoff_depth_m = self.timestep_rainfall_input_m - self.infiltration_depth_m
                    
                else:
                    
                    self.surface_runoff_depth_m = 0.0
                    
                    self.infiltration_depth_m =  self.timestep_rainfall_input_m - self.surface_runoff_depth_m
                    
        else:
            
            self.surface_runoff_depth_m = 0.0
            
            self.infiltration_depth_m = 0.0
            
        return
            
                               
    # __________________________________________________________________________________________________________
    def et_from_soil(self):
        """
            take AET from soil moisture storage, 
            using Budyko type curve to limit PET if wilting<soilmoist<field_capacity
        """
        
        if self.potential_et_m_per_timestep > 0:
            
            print("this should not happen yet. Still debugging the other functions.")
            
            if self.soil_reservoir['storage_m'] >= self.soil_reservoir['storage_threshold_primary_m']:
            
                self.actual_et_m_per_timestep = np.min(self.potential_et_m_per_timestep, 
                                                       self.soil_reservoir['storage_m'])

                self.soil_reservoir['storage_m'] -= self.actual_et_m_per_timestep

                self.et_struct['potential_et_m_per_timestep'] = 0.0
                               
            elif (self.soil_reservoir['storage_m'] > self.soil_reservoir['wilting_point_m'] and 
                  self.soil_reservoir['storage_m'] < self.soil_reservoir['storage_threshold_primary_m']):
            
                Budyko_numerator = self.soil_reservoir['storage_m'] - self.soil_reservoir['wilting_point_m']
                Budyko_denominator = self.soil_reservoir['storage_threshold_primary_m'] - \
                                     self.soil_reservoir['wilting_point_m']
                Budyki = Budyko_numerator / Budyko_denominator
                               
                self.actual_et_m_per_timestep = Budyko * self.potential_et_m_per_timestep
                               
                self.soil_reservoir['storage_m'] -= self.actual_et_m_per_timestep
        return
            
            
    # __________________________________________________________________________________________________________
    def is_fabs_less_than_epsilon(self,a,epsilon):
        
        if np.abs(a) < epsilon:
            
            return True
        
        else:
            
            return False 
    