#!/bin/bash
gcc -lm ./src/cfe_driver.c ./src/cfe.c -o run_cfe
./run_cfe
