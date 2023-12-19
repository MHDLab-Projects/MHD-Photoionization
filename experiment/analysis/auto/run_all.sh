#!/bin/bash

echo "---Experiment Overviews---"

python expt_overview.py

echo "---Absem auto---"

# Define the test cases
test_cases=("53x" "536_pos" "536_power")


# Loop over the test cases
for tc in "${test_cases[@]}"
do
    # Run the Python script with the current test case
    python absem_1d.py $tc
done

echo "---mws auto---"

#Convert these scripts to use positional args like absem? 
python mws_1d.py
python mws_2d.py
