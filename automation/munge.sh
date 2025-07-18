# Munge raw experimental data. 
# source munge.sh > output/output.txt to capture all error logs

source ../.env

mkdir -p $REPO_DIR/experiment/data/munged

# Process 2023-05-24 first as we use emulsion recipe test plan from there. 
date_list=('2023-05-24' '2023-04-07' '2023-05-12' '2023-05-18')

echo "--- Standard Post Processing --- " 
cd $REPO_DIR/PostProcessor/scripts

for date in ${date_list[@]}
do
    python post_processing.py -d $date 
done

echo " --- Additional Munging and Processing --- "
cd $REPO_DIR/experiment/data

python munge_lecroy.py
python munge_spe.py
python munge_spe_beamtiming.py

# Previously in process_munged.sh. but output is included in munged data folder
date_list=('2023-04-07' '2023-05-12' '2023-05-18' '2023-05-24')

for date in "${date_list[@]}"
do
    python aas_setup.py -d $date
    python combine_lecroy_time.py -d $date
done