#!/bin/bash
gcc -lm ./src/main_pass_forcings.c ./original_author_code/cfe.c ./src/bmi_cfe.c ../alt-modular/Modules/AORC/src/aorc.c ../alt-modular/Modules/AORC/src/bmi_aorc.c  -o run_cfe_bmi_pass_forcings
./run_cfe_bmi_pass_forcings ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt
