
source .env

mkdir -p experiment/data/munged

date_list=('2023-04-07' '2023-05-12' '2023-05-18' '2023-05-24')

echo "--- Standard Post Processing --- " 
cd $REPO_DIR/PostProcessor

for date in ${date_list[@]}
do
    python PostProcessing.py -d $date #-i absem absem_calib
done

echo " --- Additional Munging and Processing --- "
cd $REPO_DIR/experiment/data
./run_all.sh

echo " --- Rendering final notebooks --- "
cd $REPO_DIR/experiment/analysis
./render_all.sh