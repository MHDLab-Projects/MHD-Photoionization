python setup.py develop # For pi_paper_utils package

git submodule update --init

cd mhdpy
pip install -e .

cd mhdpy/fileio
git submodule update --init