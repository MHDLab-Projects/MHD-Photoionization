# source munge.sh

auto_dir=$(pwd)

source ../.env

cd $auto_dir
source render_all_expt.sh

cd $auto_dir
source render_all_expt.sh final

cd $auto_dir
source supp_scripts.sh

cd $REPO_DIR/doc
# source gen_docs.sh
python setup_repodir.py

# cd $auto_dir
# source collect_final_files.sh