#!/bin/bash
touch run_cfe
mv -f run_cfe z_trash
gcc -lm t-shirt_0.99d_2.c -o run_cfe
./run_cfe
