#!/bin/bash
gcc -lm ./src/main.c ./src/cfe.c ./src/bmi_cfe.c -o run_cfe_bmi
./run_cfe_bmi ./configs/cat_58_bmi_config_cfe.txt
