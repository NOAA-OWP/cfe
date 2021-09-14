#ifndef PET_PENMAN_MONTEITH_METHOD_H
#define PET_PENMAN_MONTEITH_METHOD_H

// FUNCTION AND SUBROUTINE PROTOTYPES

double pevapotranspiration_penman_monteith_method(pet_model *model);

double penman_monteith_pet_calculation
(
  double delta,
  double gamma,
  double moist_air_density_kg_per_m3,
  double vapor_pressure_deficit_Pa,
  struct pet_model *model
);

//############################################################*
// subroutine to calculate evapotranspiration using Penman-   *
// Monteith FAO reference ET procedure.                       *
// Reference: http://www.fao.org/3/X0490E/x0490e06.htm        *
// F.L. Ogden, NOAA National Weather Service, 2020            *
//############################################################*
double pevapotranspiration_penman_monteith_method(pet_model *model)
{
  // local varibles
  double instantaneous_pet_rate_m_per_s;
  double psychrometric_constant_Pa_per_C;
  double slope_sat_vap_press_curve_Pa_s;
  double moist_air_density_kg_per_m3;
  double water_latent_heat_of_vaporization_J_per_kg;
  double moist_air_gas_constant_J_per_kg_K;
  double moist_air_specific_humidity_kg_per_m3;
  double vapor_pressure_deficit_Pa;
  double liquid_water_density_kg_per_m3;
  double lambda_pet;
  double delta;
  double gamma;

  calculate_intermediate_variables(model);

  liquid_water_density_kg_per_m3 = model->inter_vars.liquid_water_density_kg_per_m3;
  water_latent_heat_of_vaporization_J_per_kg=model->inter_vars.water_latent_heat_of_vaporization_J_per_kg;
  vapor_pressure_deficit_Pa=model->inter_vars.vapor_pressure_deficit_Pa;
  moist_air_gas_constant_J_per_kg_K=model->inter_vars.moist_air_gas_constant_J_per_kg_K;
  moist_air_density_kg_per_m3=model->inter_vars.moist_air_density_kg_per_m3;
  slope_sat_vap_press_curve_Pa_s=model->inter_vars.slope_sat_vap_press_curve_Pa_s;
  water_latent_heat_of_vaporization_J_per_kg=model->inter_vars.water_latent_heat_of_vaporization_J_per_kg;
  psychrometric_constant_Pa_per_C=model->inter_vars.psychrometric_constant_Pa_per_C;

  delta=slope_sat_vap_press_curve_Pa_s;
  gamma=psychrometric_constant_Pa_per_C;

  if(model->pet_options.use_penman_monteith_method==TRUE)
  {
    lambda_pet = penman_monteith_pet_calculation(delta,gamma,moist_air_density_kg_per_m3,vapor_pressure_deficit_Pa,model);
  }

  instantaneous_pet_rate_m_per_s= lambda_pet/(liquid_water_density_kg_per_m3*water_latent_heat_of_vaporization_J_per_kg);

  return(instantaneous_pet_rate_m_per_s);  // meters per second

}

//#####################################################################################*
// subroutine to calculate latent heat flux (lambda*pet)                                *
// using the Penman-Monteith equation as described by FAO                              *
//   see:   http://www.fao.org/3/X0490E/x0490e06.htm#aerodynamic%20resistance%20(ra)   *
// After Allen et al. ASCE Reference ET calculation method                             *
// F.L. Ogden, NOAA National Weather Service, 2020                                     *
//#####################################################################################*
double penman_monteith_pet_calculation
(
  double delta,
  double gamma,
  double moist_air_density_kg_per_m3,
  double vapor_pressure_deficit_Pa,
  struct pet_model *model
)

{
  // local varibles
  double pm_numerator;
  double pm_denominator;
  double aerodynamic_resistance_s_per_m;

  // this method requires more calculations 
  if(is_fabs_less_than_eps(model->pet_params.vegetation_height_m,1.0e-06)==TRUE)
  {
    // the vegetation height was not specified.  TODO should warn??
    fprintf(stderr,"WARNING: Vegetation height not specified in the Penman-Monteith routine.  Using 0.5m.\n");
    model->pet_params.vegetation_height_m=0.5;  // use a reasonable assumed value
  }
  
  // use approximations from UN FAO: http://www.fao.org/3/X0490E/x0490e06.htm#aerodynamic%20resistance%20(ra)
  model->pet_params.zero_plane_displacement_height_m=2.0/3.0*model->pet_params.vegetation_height_m;
  model->pet_params.momentum_transfer_roughness_length_m=0.123*model->pet_params.vegetation_height_m;
  model->pet_params.heat_transfer_roughness_length_m=0.1*model->pet_params.vegetation_height_m;

  aerodynamic_resistance_s_per_m = calculate_aerodynamic_resistance(model);

  // all the ingredients have been prepared.  Make Penman-Monteith soufle...
  // from: http://www.fao.org/3/X0490E/x0490e06.htm#aerodynamic%20resistance%20(ra)

  pm_numerator = delta* (model->pet_forcing.net_radiation_W_per_sq_m - model->pet_forcing.ground_heat_flux_W_per_sq_m) + 
              moist_air_density_kg_per_m3 * CP *
              vapor_pressure_deficit_Pa/aerodynamic_resistance_s_per_m;

  pm_denominator = delta + gamma * (1.0+model->pet_forcing.canopy_resistance_sec_per_m/aerodynamic_resistance_s_per_m);
          
  return(pm_numerator/pm_denominator);  // Latent heat flux in Watts per sq. m., or J per (s m2)
}

#endif // PET_PENMAN_MONTEITH_METHOD_H
