#!/bin/bash
export OMP_NUM_THREADS=1
./test_FGG 100  6   16  2>&1 | tee FGG_100_6_16.dat
./test_FGG 100  12  16  2>&1 | tee FGG_100_12_16.dat
./test_FGG 1000 6   16  2>&1 | tee FGG_1000_6_16.dat
./test_FGG 1000 12  16  2>&1 | tee FGG_1000_12_16.dat