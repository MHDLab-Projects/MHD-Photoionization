# Processes munged data to form intermediate processed dataset for analysis, and then conversion to final dataset

source ../.env

expt_data_dir=$REPO_DIR/experiment/data

pushd $expt_data_dir

python collect_data.py
python proc_add_coord.py

popd