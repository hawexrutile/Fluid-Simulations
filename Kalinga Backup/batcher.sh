#!/bin/bash
# Array of angles
angles=(15.0 20.0 25.0 30.0 35.0 40.0 50.0 60.0 90.0)
repeat=3
# Path to your original batch script
original_script="script.sh"
# Iterate over the angles
for angle in "${angles[@]}"
do
    # Copy the original batch script and cpp file
    cp "$original_script" "script-$angle.sh"
    cp "M7.cpp" "M7-$angle.cpp"
    # edit the 47th line in the copied cpp file
    sed -i "47c\double theta = $angle*(M_PI/180);" "M7-$angle.cpp"
    # Edit lines 13 and 15 in the copied scriptutput-$angle" "script-$angle.sh"
    sed -i "14c\g++ -fopenmp -std=c++11 M7-$angle.cpp -o output-$angle" "script-$angle.sh"
    sed -i "16c\./output-$angle" "script-$angle.sh"
    # for each angle, run the program 3 times
    for (( i=1; i<=$repeat; i++ ))
    do
        # Submit the modified batch script to Slurm
        sbatch "script-$angle.sh"
        sleep 3
    done
    # wait for 3 seconds
done