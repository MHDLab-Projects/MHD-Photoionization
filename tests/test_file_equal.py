
from pprint import pprint
import pytest
import warnings

from mhdpy.analysis.standard_import import *
from mhdpy.fileio.tdms import TFxr


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

    dsst_old = TFxr(fp_old).as_dsst()
    dsst_new = TFxr(fp_new).as_dsst()

    return (dsst_new, dsst_old) 


# This is to stop these tests themselves from being run (expecting 'data_tuple')
from mhdpy.xr_utils.testing import test_same_groups as same_groups
from mhdpy.xr_utils.testing import test_same_channels as same_channels
from mhdpy.xr_utils.testing import test_rounded_values_equal as rounded_values_equal
from mhdpy.xr_utils.testing import test_exact_values_equal as exact_values_equal
from mhdpy.xr_utils.testing import test_output_equal as output_equal

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_same_groups_tdms(data_tuple_tdms):
    same_groups(data_tuple_tdms)

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_same_channels_tdms(data_tuple_tdms):
    same_channels(data_tuple_tdms)

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_rounded_values_equal_tdms(data_tuple_tdms):
    rounded_values_equal(data_tuple_tdms)

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_exact_values_equal_tdms(data_tuple_tdms):
    exact_values_equal(data_tuple_tdms)

@pytest.mark.parametrize('tdms_file', tdms_files)
def test_output_equal_tdms(data_tuple_tdms):
    output_equal(data_tuple_tdms)

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

    if not ds_new.equals(ds_old):
        assertion_message = f'New and old datasets are not equal \n\n--New--\n{ds_new}\n\n--Old--{ds_old}'
        raise AssertionError(assertion_message)