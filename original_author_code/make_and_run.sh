#!/bin/bash
gcc -lm cfe_driver.c -o run_cfe
gcc -lm t-shirt_0.99d_2.c -o run_cfe_d2
gcc -lm t-shirt_0.99f.c -o run_cfe_f
./run_cfe
./run_cfe_d2
./run_cfe_f
