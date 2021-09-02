#!/bin/bash
gcc -lm ./src/main_cfe_aorc_pet.c ./companion_model_code/src/pet.c ./companion_model_code/src/bmi_pet.c ./src/cfe.c ./src/bmi_cfe.c ./companion_model_code/src/aorc.c ./companion_model_code/src/bmi_aorc.c -o run_cfe_aorc_et_bmi
./run_cfe_aorc_et_bmi ./configs/cat_89_bmi_config_cfe.txt ./configs/cat_89_bmi_config_aorc.txt ./configs/cat_89_bmi_config_pet.txt
