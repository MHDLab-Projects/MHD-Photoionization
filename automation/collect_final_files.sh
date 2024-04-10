# This script collects final files from various locations and places them in the output directory.

source .env

output_dir=$REPO_DIR/final/output
mkdir -p $output_dir

# Collect final documents 

echo "Collecting final documents"


# mkdir -p output
# source pandoc_convert.sh

# cp main.pdf $output_dir/Photoionization.pdf
# cp output/Photoionization.docx $output_dir/Photoionization.docx


cd $REPO_DIR/doc/SI_man
cp SI_man.pdf $output_dir/Supporting\ Information.pdf

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
find "${output_dir}/nb_renders" -name '*.ipynb' -exec jupyter nbconvert --to pdf {} \;