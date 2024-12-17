# source munge.sh

auto_dir=$(pwd)

source ../.env

# cd $auto_dir
# source render_all_expt.sh

# cd $auto_dir
# source render_all_expt.sh final

# cd $auto_dir
# source supp_scripts.sh

cd $auto_dir
source render_figures.sh final/figures/schematics

cd $auto_dir
source render_figures.sh final/figures/SI

# Prepare for Latex PDF generation (manually done through ctrl alt b in VSCode)
cd $REPO_DIR/doc
# source gen_docs.sh
python setup_repodir.py