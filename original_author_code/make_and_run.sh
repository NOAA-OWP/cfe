#!/bin/bash
gcc -lm cfe_driver.c cfe.c -o run_cfe
gcc -lm t-shirt_0.99f.c -o run_cfe_f
gcc -lm t-shirt_0.99g.c -o run_cfe_g
gcc -lm CFE_1.0.c -o run_cfe_1
./run_cfe
./run_cfe_f
mv test.out test_f.out
./run_cfe_g
mv test.out test_g.out
./run_cfe_1
mv test.out test_cfe_1.out
