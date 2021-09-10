#ifndef PET_ENERGY_BALANCE_METHOD_H
#define PET_ENERGY_BALANCE_METHOD_H

// FUNCTION AND SUBROUTINE PROTOTYPES

double pevapotranspiration_energy_balance_method(pet_model *model);


//############################################################*
// subroutine to calculate evapotranspiration using           *
// Chow, Maidment, and Mays textbook                          *
// F.L. Ogden, NOAA National Weather Service, 2020            *
//############################################################*
double pevapotranspiration_energy_balance_method(pet_model *model)
{
  // local varibles
  double water_latent_heat_of_vaporization_J_per_kg;
  double liquid_water_density_kg_per_m3;
  double lambda_pet;
  double radiation_balance_pevapotranspiration_rate_m_per_s;

  // from FAO document: cp specific heat at constant pressure, 1.013 10-3 [MJ kg-1 �C-1],

  // IF SOIL WATER TEMPERATURE NOT PROVIDED, USE A SANE VALUE
  if(100.0 > model->pet_forcing.water_temperature_C) model->pet_forcing.water_temperature_C=22.0; // growing season

  // CALCULATE VARS NEEDED FOR THE ALL METHODS:

  liquid_water_density_kg_per_m3 = calc_liquid_water_density_kg_per_m3(model->pet_forcing.water_temperature_C); // rho_w

  water_latent_heat_of_vaporization_J_per_kg=2.501e+06-2370.0*model->pet_forcing.water_temperature_C;  // eqn 2.7.6 Chow etal. 
                                                                                              // aka 'lambda'

  // We need this in all options except for aerodynamic or Penman-Monteith methods.
  // Radiation balance is the simplest method.  Involves only radiation calculations, no aerodynamic calculations.

  lambda_pet=0.0;
  if( (model->pet_options.use_aerodynamic_method == FALSE ) && (model->pet_options.use_penman_monteith_method==FALSE) )
  {
    // This is equation 3.5.9 from Chow, Maidment, and Mays textbook.
    lambda_pet=model->pet_forcing.net_radiation_W_per_sq_m;
    radiation_balance_pevapotranspiration_rate_m_per_s=lambda_pet/
                                  (liquid_water_density_kg_per_m3*water_latent_heat_of_vaporization_J_per_kg);
  }
  return(radiation_balance_pevapotranspiration_rate_m_per_s);
}

#endif // PET_ENERGY_BALANCE_METHOD_H
