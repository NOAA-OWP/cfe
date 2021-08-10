#!/bin/bash
touch run_cfe_bmi
mv -f run_cfe_bmi z_trash
gcc -lm ./src/main.c ./src/cfe.c ./src/bmi_cfe.c -o run_cfe_bmi
./run_cfe_bmi ./configs/cat_89_bmi_config_cfe.txt
