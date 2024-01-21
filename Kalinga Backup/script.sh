#!/bin/bash
#SBATCH --job-name=XLR8-ser
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=serial_test_%j.log
#
echo "Job started at date"
date
#
# Load Intel compiler module
module load compilers/intel/parallel_studio_xe_2017.3.191
g++ -fopenmp -std=c++11 M7.cpp -o output
wait
./output
#
echo ----
echo "Job ended at date"
date
