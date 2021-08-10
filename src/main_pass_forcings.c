#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/cfe.h"
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"

#include "../../AORC/include/aorc.h"
#include "../../AORC/include/bmi_aorc.h"

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

    double *var = NULL;
    var = (double*) malloc (sizeof (double)*1);

//    printf("getting AORC from BMI and setting in ET\n");
    const char* var_name1 = "land_surface_air__temperature";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name1, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name1, &(var[0]));

    const char* var_name2 = "land_surface_air__pressure";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name2, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name2, &(var[0]));

    const char* var_name3 = "atmosphere_air_water~vapor__relative_saturation";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name3, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name3, &(var[0]));

    const char* var_name4 = "atmosphere_water__liquid_equivalent_precipitation_rate";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name4, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name4, &(var[0]));

    const char* var_name5 = "land_surface_radiation~incoming~shortwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name5, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name5, &(var[0]));

    const char* var_name6 = "land_surface_radiation~incoming~longwave__energy_flux";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name6, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name6, &(var[0]));

    const char* var_name7 = "land_surface_wind__x_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name7, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name7, &(var[0]));

    const char* var_name8 = "land_surface_wind__y_component_of_velocity";
    aorc_bmi_model->get_value(aorc_bmi_model, var_name8, &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, var_name8, &(var[0]));
}


/************************************************************************
    This main program is a mock framwork.
    This is not part of BMI, but acts as the driver that calls the model.
************************************************************************/
int
 main(int argc, const char *argv[])
{

  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=1){
    printf("make sure to include a path to the CFE config file\n");
    exit(1);
  }
  if(argc<=2){
    printf("make sure to include a path to the AORC config file\n");
    exit(1);
  }

  /************************************************************************
      allocating memory to store the entire BMI structure for CFE and AORC
  ************************************************************************/
  printf("allocating memory to store entire BMI structure for CFE\n");
  Bmi *cfe_bmi_model = (Bmi *) malloc(sizeof(Bmi));
  printf("allocating memory to store entire BMI AORC structure\n");
  Bmi *aorc_bmi_model = (Bmi *) malloc(sizeof(Bmi));

  /************************************************************************
      Registering the BMI model for CFE and AORC
  ************************************************************************/
  printf("Registering BMI CFE model\n");
  register_bmi_cfe(cfe_bmi_model);
  printf("Registering BMI AORC model\n");
  register_bmi_aorc(aorc_bmi_model);

  /************************************************************************
      Initializing the BMI model for CFE and AORC
  ************************************************************************/
  printf("Initializeing BMI CFE model\n");
  const char *cfg_file = argv[1];
  cfe_bmi_model->initialize(cfe_bmi_model, cfg_file);
  printf("Initializeing BMI AORC model\n");
  const char *cfg_file_aorc = argv[2];
  aorc_bmi_model->initialize(aorc_bmi_model, cfg_file_aorc);

  /************************************************************************
    Get the information from the configuration here in Main
  ************************************************************************/
  printf("Get the information from the configuration here in Main\n");
  cfe_model *cfe;
  cfe = (cfe_model *) cfe_bmi_model->data;
  aorc_model *aorc;
  aorc = (aorc_model *) aorc_bmi_model->data;

  /************************************************************************
    This is the basic process for getting two things to talk through BMI
    1. Update the AORC forcing data
    2. Getting forcing from AORC and setting forcing for CFE
    3. Update the CFE model.
  ************************************************************************/
  aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
  pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
  printf("precip_kg_per_m2 CFE: %e, AORC: %e \n", cfe->aorc.precip_kg_per_m2, aorc->aorc.precip_kg_per_m2);
  printf("Updating CFE: \n");
  cfe_bmi_model->update(cfe_bmi_model);                           //Update model 2
  printf("CFE streamflow out: %8.3lf\n", cfe->flux_Qout_m);

  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/
  printf("looping through and calling updata\n");
  int i;
  for (i = 1; i < 15; i++){
    aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
    pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
    printf("incoming_longwave_W_per_m2 CFE: %e, AORC: %e \n", cfe->aorc.incoming_longwave_W_per_m2, aorc->forcing_data_incoming_longwave_W_per_m2[i]);
    printf("precip_kg_per_m2 CFE: %e, AORC: %e \n", cfe->aorc.precip_kg_per_m2, aorc->aorc.precip_kg_per_m2);
    cfe_bmi_model->update(cfe_bmi_model);                           // Update model 2
    printf("flux_Qout_m %8.4e, flux_giuh_runoff_m %8.4e, flux_lat_m %8.4e, and nash_lat_runoff_m: %8.4e\n", 
            cfe->flux_Qout_m,  cfe->flux_giuh_runoff_m,  cfe->flux_lat_m,  cfe->flux_nash_lateral_runoff_m);
  }

  /************************************************************************
    Finalize both the CFE and AORC bmi models
  ************************************************************************/
  printf("Finalizing BFE and AORC models\n");
  cfe_bmi_model->finalize(cfe_bmi_model);
  aorc_bmi_model->finalize(aorc_bmi_model);

  return 0;
}

