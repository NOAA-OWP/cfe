#!/bin/bash
gcc -lm ./src/main_cfe_aorc_pet.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c -o run_cfe_aorc_et_bmi
./run_cfe_aorc_et_bmi ./configs/cat_87_bmi_config_cfe_pass.txt ./configs/cat_87_bmi_config_aorc.txt ./configs/cat_87_bmi_config_pet_pass.txt
