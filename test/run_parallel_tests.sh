#!/bin/bash
echo "Running speedup tests for N = 1000, 10000, 100000 and 1000000"
echo "This should take a while."
echo "Make sure to close all running processes to get a reliable result."
./test_parallel 1000  	100	12	2>&1 | tee parallel_speedup_1000_100_12.dat
./test_parallel 10000  	100	12	2>&1 | tee parallel_speedup_10000_100_12.dat
./test_parallel 100000 	100	12	2>&1 | tee parallel_speedup_100000_100_12.dat
./test_parallel 1000000 100	12	2>&1 | tee parallel_speedup_1000000_100_12.dat
