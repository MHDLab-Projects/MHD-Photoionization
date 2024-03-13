pip install setuptools # TODO: add to requirements?

python setup.py develop # For pi_paper_utils package

git submodule init
git submodule update

cd mhdpy
pip install -e .


cd mhdpy/fileio
git submodule update