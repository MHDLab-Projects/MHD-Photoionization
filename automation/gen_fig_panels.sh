
source ../.env

pushd $REPO_DIR/final/figure_panels

python kwt_dependence.py
python position_dependence.py
python kwt_kinetics.py
python topcam_analysis_2023-05-24.py

popd