#!/bin/bash
FILE=./run_ser_test
if test -f "$FILE"; then
    # If compilation fails, don't want to use old one
    rm "$FILE"
fi

gcc ./cfe_serialize_test.c ./serialize_state.c ../src/bmi_cfe.c ../src/cfe.c -o run_ser_test -lm -lmsgpackc 
./run_ser_test
