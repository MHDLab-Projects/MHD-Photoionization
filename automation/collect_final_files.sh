# This script collects final files from various locations and places them in the output directory.

source .env

output_dir=$REPO_DIR/final/output
mkdir -p $output_dir

# Collect final documents 

echo "Collecting final documents"

cd $REPO_DIR/doc/paper
cp main.pdf $output_dir/Photoionization\ Draft.pdf
cp output/Photoionization.docx $output_dir/Photoionization.docx
cp ../SI_man/SI_man.pdf $output_dir/Supporting\ Information.pdf


