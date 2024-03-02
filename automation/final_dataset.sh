# Runs the sequence of scripts in order to generate final datasets

source .env

pushd $REPO_DIR/final/dataset

# Generates files needed by others
echo "CFD Line profile extraction"
python extract_line_profiles_Yeq.py
echo "MWS Nothing Transmission"
python mws_nothing.py

# Final scripts
echo "Seed dependence"
python 53x_final.py
echo "Position dependence"
python position_dependence.py

popd