#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/cfe.h"
#include "../include/bmi_cfe.h"
#include "../forcing_code/include/pet.h"
#include "../forcing_code/include/bmi_pet.h"
#include "../forcing_code/include/aorc.h"
#include "../forcing_code/include/bmi_aorc.h"

/***************************************************************
    Function to pass PET to CFE using BMI.
***************************************************************/
void pass_pet_to_cfe(Bmi *cfe_bmi_model, Bmi *pet_bmi_model){
    double var_val;
    double *var_ptr = &var_val;

    pet_bmi_model->get_value(pet_bmi_model, "water_potential_evaporation_flux", &var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, "water_potential_evaporation_flux", &var_ptr);
}

/***************************************************************
    Function to pass the forcing data from AORC to CFE using BMI.
    This requires a lot of getters and setters, 
    so no need to clutter up main program
***************************************************************/
void pass_forcing_from_aorc_to_cfe(Bmi *cfe_bmi_model, Bmi *aorc_bmi_model){

    /********************************************************************
        TODO: Get variable names through BMI, then loop through those
              so we don't re-write the get/set functions over and over
    ********************************************************************/

    double var_val;
    double *var_ptr = &var_val;
    
    aorc_bmi_model->get_value(aorc_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", var_ptr);

}

/***************************************************************
    Function to pass the forcing data from AORC to PET using BMI.
    This requires a lot of getters and setters, 
    so no need to clutter up main program
***************************************************************/
void pass_forcing_from_aorc_to_pet(Bmi *pet_bmi_model, Bmi *aorc_bmi_model){

    /********************************************************************
        TODO: Get variable names through BMI, then loop through those
              so we don't re-write the get/set functions over and over
    ********************************************************************/

    double var1_val;
    double *var1_ptr = &var1_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_air__temperature", var1_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_air__temperature", var1_ptr);

    double var2_val;
    double *var2_ptr = &var2_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_air__pressure", var2_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_air__pressure", var2_ptr);

    double var3_val;
    double *var3_ptr = &var3_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "atmosphere_air_water~vapor__relative_saturation", var3_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "atmosphere_air_water~vapor__relative_saturation", var3_ptr);

    double var4_val;
    double *var4_ptr = &var4_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", var4_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", var4_ptr);

    double var5_val;
    double *var5_ptr = &var5_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_radiation~incoming~shortwave__energy_flux", var5_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_radiation~incoming~shortwave__energy_flux", var5_ptr);

    double var6_val;
    double *var6_ptr = &var6_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_radiation~incoming~longwave__energy_flux", var6_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_radiation~incoming~longwave__energy_flux", var6_ptr);

    double var7_val;
    double *var7_ptr = &var7_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_wind__x_component_of_velocity", var7_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_wind__x_component_of_velocity", var7_ptr);

    double var8_val;
    double *var8_ptr = &var8_val;
    aorc_bmi_model->get_value(aorc_bmi_model, "land_surface_wind__y_component_of_velocity", var8_ptr);
    pet_bmi_model->set_value(pet_bmi_model, "land_surface_wind__y_component_of_velocity", var8_ptr);
}

/********************************************************************
*********************************************************************
This main program is a mock framwork.
This is not part of BMI, but acts as the driver that calls the model.
*********************************************************************
*********************************************************************/
int
 main(int argc, const char *argv[])
{

  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=3){
    printf("make sure to include a path to the config files for CFE, AORC & PET\n");
    exit(1);
  }

  /************************************************************************
      allocating memory to store the entire BMI structure for CFE and AORC
  ************************************************************************/
  printf("allocating memory to store entire BMI structure for CFE\n");
  Bmi *cfe_bmi_model = (Bmi *) malloc(sizeof(Bmi));
  printf("allocating memory to store entire BMI structure for AORC\n");
  Bmi *aorc_bmi_model = (Bmi *) malloc(sizeof(Bmi));
  printf("allocating memory to store entire BMI structure for PET\n");
  Bmi *pet_bmi_model = (Bmi *) malloc(sizeof(Bmi));

  /************************************************************************
      Registering the BMI model for CFE and AORC
  ************************************************************************/
  printf("Registering BMI CFE model\n");
  register_bmi_cfe(cfe_bmi_model);
  printf("Registering BMI AORC model\n");
  register_bmi_aorc(aorc_bmi_model);
  printf("Registering BMI PET model\n");
  register_bmi_pet(pet_bmi_model);

  /************************************************************************
      Initializing the BMI model for CFE and AORC
  ************************************************************************/
  printf("Initializeing BMI CFE model\n");
  const char *cfg_file_cfe = argv[1];
  cfe_bmi_model->initialize(cfe_bmi_model, cfg_file_cfe);
  printf("Initializeing BMI AORC model\n");
  const char *cfg_file_aorc = argv[2];
  printf("AORC config file %s\n", cfg_file_aorc);
  aorc_bmi_model->initialize(aorc_bmi_model, cfg_file_aorc);
  printf("Initializeing BMI PET model\n");
  const char *cfg_file_pet = argv[3];
  pet_bmi_model->initialize(pet_bmi_model, cfg_file_pet);

  /************************************************************************
    Get the information from the configuration here in Main
  ************************************************************************/
  printf("Get the information from the configuration here in Main\n");
  cfe_state_struct *cfe1;
  cfe1 = (cfe_state_struct *) cfe_bmi_model->data;
  printf("forcing file for the CFE module %s\n", cfe1->forcing_file);
  pet_model *pet1;
  pet1 = (pet_model *) pet_bmi_model->data;
  printf("forcing file for the PET module %s\n", pet1->forcing_file);
  aorc_model *aorc1;
  aorc1 = (aorc_model *) aorc_bmi_model->data;
  printf("forcing file for the AORC module %s\n", aorc1->forcing_file);

  /************************************************************************
    This is the basic process for getting the three things to talk through BMI
    1. Update the AORC forcing data
    2. Getting forcing from AORC and setting forcing for PET
    3. Update the PET model
    3. Getting forcing from AORC and setting forcing for CFE
    4. Getting PET from PET and setting for CFE
    5. Update the CFE model.
  ************************************************************************/

  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/
  printf("looping through and calling update\n");

  if (cfe1->verbosity > 0)
    print_cfe_flux_header();

  for (int i = 0; i < 720; i++){


    aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
    
    pass_forcing_from_aorc_to_pet(pet_bmi_model, aorc_bmi_model);   // Get and Set values

    pet_bmi_model->update(pet_bmi_model);

    pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
    pass_pet_to_cfe(cfe_bmi_model, pet_bmi_model);   // Get and Set values
    
    if (pet1->aorc.air_temperature_2m_K != aorc1->aorc.air_temperature_2m_K){
      printf("ERROR: Temperature values do not match from AORC and PET\n");
      printf("Temperature value from AORC is %lf\n", aorc1->aorc.air_temperature_2m_K);
      printf("Tempterature value from PET is %lf\n", pet1->aorc.air_temperature_2m_K);
    }
    if (cfe1->aorc.precip_kg_per_m2 != aorc1->aorc.precip_kg_per_m2){
      printf("ERROR: Precip values do not match from AORC and CFE\n");
      printf("precip value from AORC is %lf\n", aorc1->aorc.precip_kg_per_m2);
      printf("precip value from CFE is %lf\n", cfe1->aorc.precip_kg_per_m2);
    }

    if (cfe1->et_struct.potential_et_m_per_s != pet1->pet_m_per_s){
      printf("ERROR: PET values do not match from PET and CFE\n");
      printf("PET value from PET is %8.9lf\n", pet1->pet_m_per_s);
      printf("PET value from CFE is %8.9lf\n", cfe1->et_struct.potential_et_m_per_s);
  }

    if (cfe1->verbosity > 2){
      printf("PET value from PET is %8.9lf\n", pet1->pet_m_per_s);
      printf("PET value from CFE is %8.9lf\n", cfe1->et_struct.potential_et_m_per_s);
  }

    cfe_bmi_model->update(cfe_bmi_model);                           //Update model 2

    if (cfe1->verbosity > 0)
      print_cfe_flux_at_timestep(cfe1);

  }

  // Run the Mass Balance check
  mass_balance_check(cfe1);

  /************************************************************************
    Finalize the CFE, PET and AORC bmi models
  ************************************************************************/
  printf("Finalizing models\n");
  cfe_bmi_model->finalize(cfe_bmi_model);
  aorc_bmi_model->finalize(aorc_bmi_model);
  pet_bmi_model->finalize(pet_bmi_model);

  return 0;
}

