
source .env

pushd $REPO_DIR/final/figure_panels

python kwt_dependence.py
python position_dependence.py

popd