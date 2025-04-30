#!/bin/bash

echo "---Experiment Overviews---"

python expt_overview.py

echo "---aas auto---"

# Define the test cases
test_cases=("53x" "536_pos" "536_power")


# Loop over the test cases
for tc in "${test_cases[@]}"
do
    # Run the Python script with the current test case
    python aas_1d.py $tc
    python mwt_1d_noise.py $tc
done

echo "---mwt auto---"

#Convert these scripts to use positional args like aas? 
python mwt_1d.py
python mwt_2d.py
