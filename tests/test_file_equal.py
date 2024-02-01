
from pprint import pprint
import pytest
import warnings

from mhdpy.analysis.standard_import import *
from mhdpy.fileio.tdms import tdms2ds


import logging
import os
import pytest
import glob

test_data_path = pjoin(REPO_DIR, 'tests', 'test_data') 
input_data_folder = pjoin(REPO_DIR, 'experiment','data','proc_data')

# Get all TDMS files in the test_data folder
tdms_files = glob.glob(os.path.join(test_data_path, '*.tdms'))
# strip test_data_path from the file names
tdms_files = [os.path.basename(f) for f in tdms_files]

# Get all cdf files in the test_data folder, and subfolders
cdf_files = glob.glob(os.path.join(test_data_path, '**/*.cdf'), recursive=True)
# strip test_data_path from the file names
cdf_files = [os.path.relpath(f, test_data_path) for f in cdf_files]

# Set up logging
logging.basicConfig(level=logging.INFO)


#TDMS tests

@pytest.fixture(scope='function')
def data_tuple_tdms(request, tdms_file) -> (xr.Dataset, xr.Dataset):

    fp_old = os.path.join(test_data_path, tdms_file)
    fp_new = os.path.join(input_data_folder, tdms_file)

    # Log the files being processed
    logging.info(f'Processing files: {fp_old}, {fp_new}')

    ds_old = tdms2ds(fp_old)
    ds_new = tdms2ds(fp_new)

    return (ds_new, ds_old) 

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_equals_tdms(data_tuple_tdms):
    ds_new, ds_old = data_tuple_tdms

    assert ds_new.equals(ds_old)

# CDF tests

@pytest.fixture(scope='function')
def data_tuple_cdf(request, cdf_file) -> (xr.Dataset, xr.Dataset):

    fp_old = os.path.join(test_data_path, cdf_file)
    fp_new = os.path.join(input_data_folder, cdf_file)

    # Log the files being processed
    logging.info(f'Processing files: {fp_old}, {fp_new}')

    ds_old = xr.load_dataset(fp_old)
    ds_new = xr.load_dataset(fp_new)

    return (ds_new, ds_old)

@pytest.mark.parametrize('cdf_file', cdf_files)
def test_equals_cdf(data_tuple_cdf):
    ds_new, ds_old = data_tuple_cdf

    assert ds_new.equals(ds_old)