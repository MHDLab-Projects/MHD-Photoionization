# Processes munged data to form intermediate processed dataset for analysis, and then conversion to final dataset

source ../.env

expt_data_dir=$REPO_DIR/experiment/data

pushd $expt_data_dir

date_list=('2023-04-07' '2023-05-12' '2023-05-18' '2023-05-24')

for date in "${date_list[@]}"
do
    python absem_setup.py -d $date
    python combine_lecroy_time.py -d $date
done

python collect_data.py
python proc_add_coord.py

popd