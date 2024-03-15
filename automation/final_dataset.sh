# Runs the sequence of scripts in order to generate final datasets

source .env

pushd $REPO_DIR/final/dataset

# Generates files needed by others
echo "CFD Line profile extraction"
python extract_line_profiles_Yeq.py
echo "MWS Nothing Transmission"
python mws_nothing.py
echo "CFD 2D Images"
python extract_cfd_2D.py


# Final scripts
echo "Seed dependence"
python 53x_final.py
echo "Position dependence"
python position_dependence.py

popd

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
