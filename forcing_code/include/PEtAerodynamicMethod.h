#ifndef PET_AERODYNAMIC_METHOD_H
#define PET_AERODYNAMIC_METHOD_H

// FUNCTION AND SUBROUTINE PROTOTYPES

double pevapotranspiration_aerodynamic_method(pet_model *model);


//############################################################*
// subroutine to calculate evapotranspiration using           *
// Chow, Maidment, and Mays textbook                          *
// F.L. Ogden, NOAA National Weather Service, 2020            *
//############################################################*
double pevapotranspiration_aerodynamic_method(pet_model *model)
{
  // local varibles
  double psychrometric_constant_Pa_per_C;
  double slope_sat_vap_press_curve_Pa_s;
  double moist_air_density_kg_per_m3;
  double water_latent_heat_of_vaporization_J_per_kg;
  double moist_air_gas_constant_J_per_kg_K;
  double vapor_pressure_deficit_Pa;
  double liquid_water_density_kg_per_m3;
  double aerodynamic_method_pevapotranspiration_rate_m_per_s;
  double mass_flux;
  double von_karman_constant_squared=(double)KV2;  // a constant equal to 0.41 squared

  calculate_intermediate_variables(model);

  liquid_water_density_kg_per_m3 = model->inter_vars.liquid_water_density_kg_per_m3;
  water_latent_heat_of_vaporization_J_per_kg=model->inter_vars.water_latent_heat_of_vaporization_J_per_kg;
  vapor_pressure_deficit_Pa=model->inter_vars.vapor_pressure_deficit_Pa;
  moist_air_gas_constant_J_per_kg_K=model->inter_vars.moist_air_gas_constant_J_per_kg_K;
  moist_air_density_kg_per_m3=model->inter_vars.moist_air_density_kg_per_m3;
  slope_sat_vap_press_curve_Pa_s=model->inter_vars.slope_sat_vap_press_curve_Pa_s;
  psychrometric_constant_Pa_per_C=model->inter_vars.psychrometric_constant_Pa_per_C;

  if( model->pet_options.use_penman_monteith_method == FALSE)  // we don't use this term in Penman-Monteith method
  {
    if (model->bmi.verbose >1)
      printf("Use Penman Monteith method is FALSE\n");

    // This is equation 3.5.16 from Chow, Maidment, and Mays textbook.
    mass_flux = 0.622*von_karman_constant_squared*moist_air_density_kg_per_m3*      // kg per sq. meter per sec.
                vapor_pressure_deficit_Pa*model->pet_forcing.wind_speed_m_per_s/
                (model->pet_forcing.air_pressure_Pa*
                pow(log(model->pet_params.wind_speed_measurement_height_m/model->pet_params.zero_plane_displacement_height_m),2.0));
    aerodynamic_method_pevapotranspiration_rate_m_per_s=mass_flux/liquid_water_density_kg_per_m3;  
  }

  return(aerodynamic_method_pevapotranspiration_rate_m_per_s);
}

#endif // PET_AERODYNAMIC_METHOD_H
