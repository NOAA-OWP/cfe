#!/bin/bash
gcc -lm cfe_driver.c cfe.c -o run_cfe
gcc -lm t-shirt_0.99f.c -o run_cfe_f
./run_cfe
./run_cfe_f
