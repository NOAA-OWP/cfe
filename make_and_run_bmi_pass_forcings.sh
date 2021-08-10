#!/bin/bash
touch run_cfe_pass_forcings
mv -f run_cfe_pass_forcings z_trash
gcc -lm ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ../AORC/src/aorc.c ../AORC/src/bmi_aorc.c  -o run_cfe_bmi_pass_forcings
./run_cfe_bmi_pass_forcings ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt
