#!/bin/bash
gcc -lm -Wall ./src/main.c ./src/bmi_cfe.c ./src/cfe.c -o run_cfe_bmi
./run_cfe_bmi ./configs/cat_87_bmi_config_cfe.txt
