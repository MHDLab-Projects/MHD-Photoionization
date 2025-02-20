# This script collects final files from various locations and places them in the output directory.

source ../.env

rm -rf $REPO_DIR/output


# Collect final documents 

echo "Collecting documents"

output_dir=$REPO_DIR/output
mkdir -p $output_dir

cp $REPO_DIR/final/figures/output $output_dir/figures -r

cd $REPO_DIR/doc/SI
cp SI.pdf $output_dir/Supporting\ Information.pdf

#final\sim_input\output\sim_input_mean_table.csv
cp $REPO_DIR/final/sim_input/output/sim_input_mean_table.csv $output_dir/sim_input_mean_table.csv

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