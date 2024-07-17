# source munge.sh

auto_dir=$(pwd)

source ../.env

source munge.sh

cd $auto_dir
source process_munged.sh

cd $auto_dir
source final_dataset.sh

cd $auto_dir
source render_all_expt.sh

cd $auto_dir
source render_all_expt.sh final

cd $auto_dir
source supp_scripts.sh

cd $auto_dir
source gen_fig_panels.sh

cd $auto_dir
source render_figures.sh

cd $REPO_DIR/doc
# source gen_docs.sh
python setup_repodir.py

# cd $auto_dir
# source collect_final_files.sh