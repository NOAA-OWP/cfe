#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/cfe.h"
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"

#include "../forcing_code/include/aorc.h"
#include "../forcing_code/include/bmi_aorc.h"

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
    
    aorc_bmi_model->get_value(aorc_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", &(var[0]));
    cfe_bmi_model->set_value(cfe_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", &(var[0]));

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
  cfe_state_struct *cfe_model_data;
  cfe_model_data = (cfe_state_struct *) cfe_bmi_model->data;
  aorc_model *aorc;
  aorc = (aorc_model *) aorc_bmi_model->data;

  /*************************************************************************
   This will be used to advance the model with update_until
  **************************************************************************/
  double model_time_step_size;
  cfe_bmi_model->get_time_step(cfe_bmi_model, &model_time_step_size);
  printf("The model time step size is: %lf\n", model_time_step_size);

  /************************************************************************
    This is the basic process for getting two things to talk through BMI
    1. Update the AORC forcing data
    2. Getting forcing from AORC and setting forcing for CFE
    3. Update the CFE model.
  ************************************************************************/
  
  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/
  printf("looping through and calling updata\n");
  if (cfe_model_data->verbosity > 0)
    print_cfe_flux_header();
  for (int i = 0; i < 45; i++){

    aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
  
    pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
  
    if (cfe_model_data->aorc.precip_kg_per_m2 != aorc->aorc.precip_kg_per_m2){
      printf("Precip values do not match\n");
      printf("precip value from AORC is %lf\n", aorc->aorc.precip_kg_per_m2);
      printf("precip value from CFE is %lf\n", cfe_model_data->aorc.precip_kg_per_m2);
    }

    cfe_bmi_model->update_until(cfe_bmi_model, (i+1)*model_time_step_size);                           // Update model 2
  
    if (cfe_model_data->verbosity > 0)
      print_cfe_flux_at_timestep(cfe_model_data);
  }

  // Run the Mass Balance check
  mass_balance_check(cfe_model_data);

  /************************************************************************
    Finalize both the CFE and AORC bmi models
  ************************************************************************/
  printf("Finalizing BFE and AORC models\n");
  cfe_bmi_model->finalize(cfe_bmi_model);
  aorc_bmi_model->finalize(aorc_bmi_model);

  return 0;
}

