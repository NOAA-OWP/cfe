import time
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import math

def conceptual_reservoir_flux_calc(t, S, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc):
    # ODE of soil moisture reservoir mass balance
    """
    Case 1: S (Soil moisture storage ) > storage_threshold_primary_m
        Interpretation: When the soil moisture is plenty, AET(=PET), percolation, and lateral flow are all active. Thus,
        Equation: dS/dt = Infiltration - PET - (Klf+Kperc) * (S - storage_threshold_primary_m)/(storage_max_m - storage_threshold_primary_m)
    Case 2: storage_threshold_primary_m > S (Soil moisture storage) > storage_threshold_primary_m - wltsmc
        Interpretation: When the soil moisture is in the middle range, AET is active and proportional to the soil moisture storage ratio
        Equation: dS/dt = Infiltration - PET * (S - wltsmc)/(storage_threshold_primary_m - wltsmc)
    Case 3: wltsmc > S (Soil moisture storage)
        Interpretation: When the soil moisture is depleted, no outflux is active
        Equation: dS/dt = Infitlation
    :param t:
    :param S:
    :param storage_threshold_primary_m:
    :param storage_max_m: maximum soil moisture storage, i.e., porosity
    :param coeff_primary: K_perc, percolation coefficient
    :param coeff_secondary: K_lf, lateral flow coefficient
    :param PET: potential evapotranspiration
    :param infilt: infiltration
    :param wltsmc: wilting point (in meter)
    :return: dS
    """

    storage_above_threshold_m = S - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply(S - wltsmc > 0, 1)

    storage_above_threshold_m_paw = S - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = np.minimum(storage_above_threshold_m_paw/storage_diff_paw, 1) # Equation 11 (Ogden's document).

    dS = infilt -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * storage_ratio - ET_switch * PET * storage_ratio_paw
    return dS

def jac(t, S, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc):
    # Jacobian matrix for the conceptual_reservoir_flux_calc
    storage_diff = storage_max_m - storage_threshold_primary_m

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply((S - wltsmc > 0) and (S - storage_threshold_primary_m < 0), 1)

    storage_diff_paw = storage_threshold_primary_m - wltsmc

    dfdS = -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * 1/storage_diff - ET_switch * PET * 1/storage_diff_paw
    return [dfdS]

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
        cfe_state.reduced_potential_et_m_per_timestep =  cfe_state.potential_et_m_per_s * cfe_state.time_step_size
        cfe_state.vol_PET += cfe_state.potential_et_m_per_timestep

        # ________________________________________________
        # SUBROUTINE
        # timestep_rainfall_input_m = f(timestep_rainfall_input_m, potential_et_m_per_timestep)
        cfe_state.actual_et_from_rain_m_per_timestep = 0
        self.et_from_rainfall(cfe_state)

        cfe_state.vol_et_from_rain += cfe_state.actual_et_from_rain_m_per_timestep
        cfe_state.vol_et_to_atm += cfe_state.actual_et_from_rain_m_per_timestep
        cfe_state.volout += cfe_state.actual_et_from_rain_m_per_timestep
        
        # ________________________________________________        
        cfe_state.soil_reservoir_storage_deficit_m = (cfe_state.soil_params['smcmax'] * \
                                                 cfe_state.soil_params['D'] - \
                                                 cfe_state.soil_reservoir['storage_m'])
        
        # ________________________________________________
        # SUBROUTINE
        # Calculates the value for surface_runoff_depth_m (infiltration excess overland flow)
        self.Schaake_partitioning_scheme(cfe_state)
        cfe_state.vol_sch_runoff += cfe_state.surface_runoff_depth_m
        
        # ________________________________________________
        if cfe_state.soil_reservoir_storage_deficit_m < cfe_state.infiltration_depth_m:
            # Calculates saturation excess overland flow
            cfe_state.surface_runoff_depth_m += (cfe_state.infiltration_depth_m - cfe_state.soil_reservoir_storage_deficit_m)
            cfe_state.infiltration_depth_m = cfe_state.soil_reservoir_storage_deficit_m
        # ________________________________________________
        cfe_state.vol_sch_infilt += cfe_state.infiltration_depth_m

        # ________________________________________________
        # Solve ODE for the soil reservoir
        # Calculates primary_flux, secondary_flux, AET, and infiltration simultaneously
        cfe_state.actual_et_from_soil_m_per_timestep = 0
        self.soil_reservoir_flux_calc(cfe_state, cfe_state.soil_reservoir)

        cfe_state.flux_perc_m = cfe_state.primary_flux_m
        cfe_state.flux_lat_m = cfe_state.secondary_flux_m
        cfe_state.vol_et_from_soil += cfe_state.actual_et_from_soil_m_per_timestep
        cfe_state.vol_et_to_atm += cfe_state.actual_et_from_soil_m_per_timestep
        cfe_state.volout += cfe_state.actual_et_from_soil_m_per_timestep

        # ________________________________________________
        cfe_state.gw_reservoir_storage_deficit_m = cfe_state.gw_reservoir['storage_max_m'] - cfe_state.gw_reservoir['storage_m']

        # ________________________________________________
        if cfe_state.flux_perc_m > cfe_state.gw_reservoir_storage_deficit_m:
            # (?) When the groundwater storage is full, the overflowing amount goes to direct runoff?
            self.diff_perc = cfe_state.flux_perc_m - cfe_state.gw_reservoir_storage_deficit_m
            cfe_state.flux_perc_m = cfe_state.gw_reservoir_storage_deficit_m
            cfe_state.vol_sch_runoff += self.diff_perc
            cfe_state.vol_sch_infilt -= self.diff_perc
            
        # ________________________________________________
        cfe_state.gw_reservoir['storage_m']     += cfe_state.flux_perc_m
        cfe_state.vol_to_gw                     += cfe_state.flux_perc_m
        cfe_state.vol_soil_to_gw                += cfe_state.flux_perc_m
        cfe_state.vol_soil_to_lat_flow          += cfe_state.flux_lat_m  #TODO add this to nash cascade as input
        cfe_state.volout                        += cfe_state.flux_lat_m

        # ________________________________________________
        # SUBROUTINE
        # primary_flux, secondary_flux = f(reservoir)
        self.conceptual_reservoir_flux_calc(cfe_state, cfe_state.gw_reservoir) 
            
        # ________________________________________________
        cfe_state.flux_from_deep_gw_to_chan_m = cfe_state.primary_flux_m
        cfe_state.vol_from_gw += cfe_state.flux_from_deep_gw_to_chan_m
        
        # ________________________________________________        
        if not self.is_fabs_less_than_epsilon(cfe_state.secondary_flux_m, 1.0e-09):
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

        if cfe_state.timestep_rainfall_input_m > 0.0:
            if cfe_state.timestep_rainfall_input_m > cfe_state.potential_et_m_per_timestep:

                cfe_state.actual_et_from_rain_m_per_timestep = cfe_state.potential_et_m_per_timestep
                cfe_state.timestep_rainfall_input_m -= cfe_state.actual_et_from_rain_m_per_timestep

            else:
                cfe_state.actual_et_from_rain_m_per_timestep = cfe_state.timestep_rainfall_input_m
                cfe_state.timestep_rainfall_input_m = 0.0

        cfe_state.reduced_potential_et_m_per_timestep = cfe_state.potential_et_m_per_timestep - cfe_state.actual_et_from_rain_m_per_timestep

        return
                

    def soil_reservoir_flux_calc(self, cfe_state, reservoir):
        """
        This function solves ODE for a soil reservoir
        :param cfe_state:
        :param reservoir:
        :return: reservoir['storage_m']
        """

        # Initialization
        y0 = [reservoir['storage_m']]

        t = np.array([0, 0.05, 0.15, 0.3, 0.6, 1.0])

        # Solve ODE
        sol = odeint(
            conceptual_reservoir_flux_calc,
            y0,
            t,
            args=(
                reservoir['storage_threshold_primary_m'],
                reservoir['storage_max_m'],
                reservoir['coeff_primary'],
                reservoir['coeff_secondary'],
                cfe_state.reduced_potential_et_m_per_timestep,
                cfe_state.infiltration_depth_m,
                cfe_state.soil_params['wltsmc'] * cfe_state.soil_params['D']
            ),
            tfirst=True,
            Dfun=jac
        )

        # Finalize the results
        ts_concat = t
        ys_concat = np.concatenate(sol, axis=0)

        # Calculate fluxes (fine-tuning the residuals by scaling)
        t_proportion = np.diff(ts_concat)
        ys_avg = np.convolve(ys_concat, np.ones(2), 'valid') / 2

        lateral_flux = np.zeros(ys_avg.shape)
        perc_lat_switch = ys_avg - reservoir['storage_threshold_primary_m'] > 0
        lateral_flux[perc_lat_switch] = reservoir['coeff_secondary'] * np.minimum(
            (ys_avg[perc_lat_switch] - reservoir['storage_threshold_primary_m']) / (
                        reservoir['storage_max_m'] - reservoir['storage_threshold_primary_m']), 1)
        lateral_flux_frac = lateral_flux * t_proportion

        perc_flux = np.zeros(ys_avg.shape)
        perc_flux[perc_lat_switch] = reservoir['coeff_primary'] * np.minimum(
            (ys_avg[perc_lat_switch] - reservoir['storage_threshold_primary_m']) / (
                        reservoir['storage_max_m'] - reservoir['storage_threshold_primary_m']), 1)
        perc_flux_frac = perc_flux * t_proportion

        et_from_soil = np.zeros(ys_avg.shape)
        ET_switch = ys_avg - cfe_state.soil_params['wltsmc'] * cfe_state.soil_params['D']> 0
        et_from_soil[ET_switch] = cfe_state.reduced_potential_et_m_per_timestep * np.minimum(
            (ys_avg[ET_switch] - cfe_state.soil_params['wltsmc']* cfe_state.soil_params['D']) / (reservoir['storage_threshold_primary_m'] - cfe_state.soil_params['wltsmc']* cfe_state.soil_params['D']), 1)
        et_from_soil_frac = et_from_soil * t_proportion

        infilt_to_soil = np.repeat(cfe_state.infiltration_depth_m, ys_avg.shape)
        infilt_to_soil_frac = infilt_to_soil * t_proportion

        # Scale fluxes
        sum_outflux = lateral_flux_frac + perc_flux_frac + et_from_soil_frac
        if sum_outflux.any() == 0:
            flux_scale = 0
        else:
            flux_scale = np.zeros(infilt_to_soil_frac.shape)
            flux_scale[sum_outflux != 0] = (np.diff(-ys_concat, axis=0)[sum_outflux != 0] + infilt_to_soil_frac[
                sum_outflux != 0]) / sum_outflux[sum_outflux != 0]
            flux_scale[sum_outflux == 0] = 0
        scaled_lateral_flux = lateral_flux_frac * flux_scale
        scaled_perc_flux = perc_flux_frac * flux_scale
        scaled_et_flux = et_from_soil_frac * flux_scale

        # Return the results
        cfe_state.primary_flux_m = math.fsum(scaled_perc_flux)
        cfe_state.secondary_flux_m = math.fsum(scaled_lateral_flux)
        cfe_state.actual_et_from_soil_m_per_timestep = math.fsum(scaled_et_flux)
        reservoir['storage_m'] = ys_concat[-1]

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
    def is_fabs_less_than_epsilon(self,a,epsilon):
        
        if np.abs(a) < epsilon:
            
            return True
        
        else:
            
            return False 
    
