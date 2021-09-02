#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/cfe.h"
#include "../include/bmi_cfe.h"
#include "../companion_model_code/include/pet.h"
#include "../companion_model_code/include/bmi_pet.h"
#include "../companion_model_code/include/aorc.h"
#include "../companion_model_code/include/bmi_aorc.h"

/***************************************************************
    Function to pass PET to CFE using BMI.
***************************************************************/
void pass_pet_to_cfe(Bmi *cfe_bmi_model, Bmi *pet_bmi_model){
    double var_val;
    double *var_ptr = &var_val;

    const char* var_name1 = "water_potential_evaporation_flux";
    pet_bmi_model->get_value(pet_bmi_model, var_name1, &var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name1, &var_ptr);
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

//    printf("getting AORC from BMI and setting in ET\n");
    const char* var_name1 = "land_surface_air__temperature";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name1, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name1, var_ptr);

    const char* var_name2 = "land_surface_air__pressure";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name2, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name2, var_ptr);

    const char* var_name3 = "atmosphere_air_water~vapor__relative_saturation";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name3, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name3, var_ptr);

    const char* var_name4 = "atmosphere_water__time_integral_of_precipitation_mass_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name4, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name4, var_ptr);

    const char* var_name5 = "land_surface_radiation~incoming~shortwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name5, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name5, var_ptr);

    const char* var_name6 = "land_surface_radiation~incoming~longwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name6, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name6, var_ptr);

    const char* var_name7 = "land_surface_wind__x_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name7, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name7, var_ptr);

    const char* var_name8 = "land_surface_wind__y_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name8, var_ptr);
    cfe_bmi_model->set_value(cfe_bmi_model, var_name8, var_ptr);
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

    double var_val;
    double *var_ptr = &var_val;

//    printf("getting AORC from BMI and setting in PET\n");
    const char* var_name1 = "land_surface_air__temperature";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name1, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name1, var_ptr);

    const char* var_name2 = "land_surface_air__pressure";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name2, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name2, var_ptr);

    const char* var_name3 = "atmosphere_air_water~vapor__relative_saturation";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name3, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name3, var_ptr);

    const char* var_name4 = "atmosphere_water__time_integral_of_precipitation_mass_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name4, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name4, var_ptr);

    const char* var_name5 = "land_surface_radiation~incoming~shortwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name5, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name5, var_ptr);

    const char* var_name6 = "land_surface_radiation~incoming~longwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name6, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name6, var_ptr);

    const char* var_name7 = "land_surface_wind__x_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name7, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name7, var_ptr);

    const char* var_name8 = "land_surface_wind__y_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name8, var_ptr);
    pet_bmi_model->set_value(pet_bmi_model, var_name8, var_ptr);
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
  cfe_state_struct *cfe;
  cfe = (cfe_state_struct *) cfe_bmi_model->data;
  printf("forcing file for the CFE module %s\n", cfe->forcing_file);
  pet_model *pet;
  pet = (pet_model *) pet_bmi_model->data;
  printf("forcing file for the PET module %s\n", pet->forcing_file);
  aorc_model *aorc;
  aorc = (aorc_model *) aorc_bmi_model->data;
  printf("forcing file for the AORC module %s\n", aorc->forcing_file);

  /************************************************************************
    This is the basic process for getting the three things to talk through BMI
    1. Update the AORC forcing data
    2. Getting forcing from AORC and setting forcing for PET
    3. Update the PET model
    3. Getting forcing from AORC and setting forcing for CFE
    4. Getting PET from PET and setting for CFE
    5. Update the CFE model.
  ************************************************************************/

  printf("Updating AORC: \n");
  aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
  pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
  pass_forcing_from_aorc_to_pet(pet_bmi_model, aorc_bmi_model);   // Get and Set values
  printf("Updating PET: \n");
  pet_bmi_model->update(pet_bmi_model);
  pass_pet_to_cfe(cfe_bmi_model, pet_bmi_model);   // Get and Set values
  printf("precip_kg_per_m2 CFE: %e, AORC: %e \n", cfe->aorc.precip_kg_per_m2, aorc->aorc.precip_kg_per_m2);
  printf("PET m per s CFE: %e, PET: %e \n", cfe->et_struct.potential_et_m_per_s, pet->pet_m_per_s);
  printf("Updating CFE: \n");
  cfe_bmi_model->update(cfe_bmi_model);                           //Update model 2
  printf("CFE streamflow out: %8.3lf\n", cfe->flux_Qout_m);
  printf("PET value from CFE %8.6e\n", cfe->et_struct.potential_et_m_per_timestep);
  printf("AET value from CFE %8.6e\n", cfe->et_struct.actual_et_m_per_timestep);

//  printf("looping through and calling update\n");
  int i=0;
  for (i = 0; i < 719; i++){


    aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
    pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
    pass_forcing_from_aorc_to_pet(pet_bmi_model, pet_bmi_model);   // Get and Set values
    pet_bmi_model->update(pet_bmi_model);
    pass_pet_to_cfe(cfe_bmi_model, pet_bmi_model);   // Get and Set values
    printf("precip_kg_per_m2 CFE: %e, AORC: %e \n", cfe->aorc.precip_kg_per_m2, aorc->aorc.precip_kg_per_m2);
    printf("PET m per s CFE: %e, PET: %e \n", cfe->et_struct.potential_et_m_per_s, pet->pet_m_per_s);
    cfe_bmi_model->update(cfe_bmi_model);                           //Update model 2

    printf("Qout, soil_res_storage_m, PET and AET values from CFE %lf %lf %lf %lf\n", 
            *cfe->flux_Qout_m, 
            cfe->soil_reservoir.storage_m,
            cfe->et_struct.potential_et_m_per_timestep, 
            cfe->et_struct.actual_et_m_per_timestep);
//    printf("--------------------- END OF TIMESTEP ---------------------------------\n");
  }

  printf("Finalizing models\n");
  cfe_bmi_model->finalize(cfe_bmi_model);
  aorc_bmi_model->finalize(aorc_bmi_model);
  pet_bmi_model->finalize(pet_bmi_model);

  return 0;
}

