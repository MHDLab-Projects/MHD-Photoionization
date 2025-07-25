# supplementary/provisional scripts
# Scripts located throughout the repo that need to be run in addition to final_dataset.sh for paper and SI
# Intended that over time all included in paepr will move to the final dataset folder as things are dialed in. 
# a large part of the supplementary figure generation happens during the rendering of analysis scripts into notebooks. This catches anything that is not in a rendered notebook.


source ../.env

## SI Auxilliary datasets ##

pushd $REPO_DIR/experiment/analysis/auto
# python expt_overview.py
# python mwt_1d.py
source run_all.sh
popd 

pushd $REPO_DIR/experiment/analysis/mwt/resampling
python process.py
python mwt_sample_nors_compare.py
popd

pushd $REPO_DIR/experiment/notebook/2023-05-18/mwt
python analysis_silicon_power.py
popd

pushd $REPO_DIR/final/sim_input
source make_sim_input.sh
popd