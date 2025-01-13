python setup.py develop # For pi_paper_utils package

git submodule update --init

cd mhdlab
pip install -e .

cd mhdlab/fileio
git submodule update --init