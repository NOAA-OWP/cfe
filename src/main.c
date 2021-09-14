#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"
#include "../include/cfe.h"

/*
This main program is a mock framwork.
This is not part of BMI, but acts as the driver that calls the model.
*/
int
 main(int argc, const char *argv[])
{

  // A configuration file is required for running this model through BMI
  if(argc<=1){
    printf("make sure to include a path to the config file\n");
    exit(1);
  }

   printf("allocating memory to store entire BMI structure for CFE\n");
  // allocating memory to store the entire BMI structure for CFE
  Bmi *cfe_bmi_model = (Bmi *) malloc(sizeof(Bmi));

  printf("Registering BMI CFE model\n");
  register_bmi_cfe(cfe_bmi_model);

  printf("Initializeing BMI CFE model\n");
  const char *cfg_file = argv[1];
  cfe_bmi_model->initialize(cfe_bmi_model, cfg_file);

  printf("Get the information from the configuration here in Main\n");
  // Get the information from the configuration here in Main  
  cfe_state_struct *cfe_main_data;
  cfe_main_data = (cfe_state_struct *) cfe_bmi_model->data;

  printf("looping through and calling update\n");
  if (cfe_main_data->verbosity > 0)
      print_cfe_flux_header();

  for (int i = 0; i < 720; i++){

    cfe_bmi_model->update(cfe_bmi_model);

    if (cfe_main_data->verbosity > 0)
        print_cfe_flux_at_timestep(cfe_main_data);

  }

  // Run the Mass Balance check
  mass_balance_check(cfe_main_data);

  printf("Finalizing CFE model\n");
  cfe_bmi_model->finalize(cfe_bmi_model);

  return 0;
}
