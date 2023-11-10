
source ../.env

echo "---Experiment Date (Notebook) Processing---"
cd $REPO_DIR/experiment/notebook
source process_all.sh


cd $REPO_DIR/experiment/data

echo "--- Consolidating Data ---"

python collect_data.py
python proc_add_coord.py

cd $REPO_DIR/experiment/analysis

echo "--- Rendering Notebooks ---"
source render_all.sh
