pip install -r requirements.txt

cd mhdpy
pip install -e .

cd ..
cd modeling

python setup.py develop

cd mhdcantera

python setup.py develop

