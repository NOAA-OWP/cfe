#!/bin/bash
gcc -lm main_unit_test_bmi.c ../src/bmi_cfe.c ../src/cfe.c -o run_cfe_bmi_test
./run_cfe_bmi_test ../configs/cat_89_bmi_config_cfe_unit_test.txt
#./run_cfe_bmi_test ../configs/cat_89_bmi_config_cfe.txt