# This script collects final files from various locations and places them in the output directory.

source ../.env

# collect datasets for processed input data (reproducible from raw data, but can be added to repo to reproduce results to avoid reprocessing time)
# Not adding experiment\data\munged as the files are very large. Manually copy that folder if needed. Maybe have a separate script to copy all datasets together (to external drive)

echo "Collecting datasets"

output_dir=$REPO_DIR/output/datasets
mkdir -p $output_dir

mkdir -p $output_dir/final/dataset
cp -r $REPO_DIR/final/dataset/output $output_dir/final/dataset/

# This data should be the same as the final dataset, but addding here in case wanting to just drag only this in and test the pipeline...
mkdir -p $output_dir/tests/test_data_final 
cp -r $REPO_DIR/tests/test_data_final $output_dir/tests -r


# Collect final documents 

echo "Collecting documents"

output_dir=$REPO_DIR/output/docs
mkdir -p $output_dir

cp $REPO_DIR/final/figures/output $output_dir/figures -r

cd $REPO_DIR/doc/SI
cp SI.pdf $output_dir/Supporting\ Information.pdf

# Remove the output_dir/nb_renders directory
rm -rf "${output_dir}/nb_renders"

# Walk through directories in REPO_DIR
find "$REPO_DIR" -type d -name 'nb_render' -print0 | while IFS= read -r -d '' dir; do
    # Construct the destination directory path
    dest="${output_dir}/nb_renders/${dir#$REPO_DIR}"

    # Create the destination directory
    mkdir -p "$dest"

    # Copy the nb_render directory to the destination directory
    cp -r "$dir/"* "$dest/"
done

# Convert notebooks to PDF
# find "${output_dir}/nb_renders" -name '*.ipynb' -exec jupyter nbconvert --to pdf {} \;