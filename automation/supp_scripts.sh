# supplementary/provisional scripts
# Scripts located throughout the repo that need to be run in addition to final_dataset.sh for paper and SI
# Intended that over time all included in paepr will move to the final dataset folder as things are dialed in. 
# a large part of the supplementary figure generation happens during the rendering of analysis scripts into notebooks. This catches anything that is not in a rendered notebook.


source .env
### Auxilliary datasets ###

#TODO: move below into final dataset folder?
pushd $REPO_DIR/modeling/viability/dataset
source calc_all.sh
popd

pushd $REPO_DIR/modeling/viability/figure_panels
source gen_figs.sh
popd

pushd $REPO_DIR/experiment/notebook/2018-11-20
python munge_spe.py
python analysis_beam_timing.py
popd

pushd $REPO_DIR/experiment/analysis/topcam
python spe_perspective.py
python spe_analysis_2023-05-24.py
popd


## SI Auxilliary datasets ##

pushd $REPO_DIR/experiment/analysis/auto
python expt_overview.py
popd 

pushd $REPO_DIR/experiment/analysis/mws/resampling
python process.py
python mws_sample_nors_compare.py
popd

pushd $REPO_DIR/experiment/analysis/topcam
python spe_analysis_2023-05-18.py
popd

pushd $REPO_DIR/experiment/notebook/2023-05-18/mws
python analysis_silicon_power.py
popd