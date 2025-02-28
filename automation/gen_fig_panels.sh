
source ../.env

pushd $REPO_DIR/final/figure_panels

python kwt_dependence.py
python position_dependence.py
python topcam_analysis_2023-05-24.py
python 536_power.py

popd

pushd $REPO_DIR/modeling/viability/figure_panels

python gamma_panels.py

popd