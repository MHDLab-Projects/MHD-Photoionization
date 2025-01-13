
from pprint import pprint
import pytest
import warnings

from mhdlab.analysis.standard_import import *
from mhdlab.fileio.tdms import TFxr


import logging
import os
import pytest
import glob

#TODO: had to do hacky things with filepaths to get multiple test data folders to work. Need to fix this

test_data_path_exp = pjoin(REPO_DIR, 'tests', 'test_data') 
input_data_folder_exp = pjoin(REPO_DIR, 'experiment','data','proc_data')

test_data_path_modeling = pjoin(REPO_DIR, 'tests', 'test_data_modeling')
input_data_folder_modeling = pjoin(REPO_DIR, 'modeling','viability','dataset', 'output')

test_data_path_final = pjoin(REPO_DIR, 'tests', 'test_data_final')
#TODO: final figure panels should be both experiment and modeling
input_data_folder_final = pjoin(REPO_DIR, 'final', 'dataset', 'output')

# Get all TDMS files in the test_data folder
tdms_files = glob.glob(os.path.join(test_data_path_exp, '*.tdms'))
# strip test_data_path from the file names
tdms_files = [os.path.basename(f) for f in tdms_files]

# Get all cdf files in the test_data folder, and subfolders
cdf_files = glob.glob(os.path.join(test_data_path_exp, '**/*.cdf'), recursive=True)
# strip test_data_path from the file names
cdf_files = [os.path.relpath(f, test_data_path_exp) for f in cdf_files]

cdf_files_modeling = glob.glob(pjoin(test_data_path_modeling, '**/*.cdf'), recursive=True)
cdf_files_modeling = [os.path.relpath(f, test_data_path_modeling) for f in cdf_files_modeling]


cdf_files_final = glob.glob(pjoin(test_data_path_final, '**/*.cdf'), recursive=True)
cdf_files_final = [os.path.relpath(f, test_data_path_final) for f in cdf_files_final]


# Set up logging
logging.basicConfig(level=logging.INFO)


#TDMS tests

@pytest.fixture(scope='function')
def data_tuple_tdms(request, tdms_file) -> (xr.Dataset, xr.Dataset):

    fp_old = os.path.join(test_data_path_exp, tdms_file)
    fp_new = os.path.join(input_data_folder_exp, tdms_file)

    # Log the files being processed
    logging.info(f'Processing files: {fp_old}, {fp_new}')

    dsst_old = TFxr(fp_old).as_dsst(convert_to_PT=False)
    dsst_new = TFxr(fp_new).as_dsst(convert_to_PT=False)

    return (dsst_new, dsst_old) 


# This is to stop these tests themselves from being run (expecting 'data_tuple')
from mhdlab.xr_utils.testing import test_same_groups as same_groups
from mhdlab.xr_utils.testing import test_same_channels as same_channels
from mhdlab.xr_utils.testing import test_rounded_values_equal as rounded_values_equal
from mhdlab.xr_utils.testing import test_exact_values_equal as exact_values_equal
from mhdlab.xr_utils.testing import test_output_equal as output_equal

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
def data_tuple_cdf(request, cdf_file, new_data_dir, test_dir) -> (xr.Dataset, xr.Dataset):

    fp_test = os.path.join(test_dir, cdf_file)
    fp_new = os.path.join(new_data_dir, cdf_file)

    # Log the files being processed
    logging.info(f'Processing files: {fp_test}, {fp_new}')

    ds_old = xr.load_dataset(fp_test)
    ds_new = xr.load_dataset(fp_new)

    return (ds_new, ds_old)


params_exp = [(cdf_file, input_data_folder_exp, test_data_path_exp) for cdf_file in cdf_files]

@pytest.mark.parametrize('cdf_file,new_data_dir,test_dir', params_exp)
def test_equals_cdf_exp(data_tuple_cdf):
    ds_new, ds_old = data_tuple_cdf

    if not ds_new.equals(ds_old):
        assertion_message = f'New and old datasets are not equal \n\n--New--\n{ds_new}\n\n--Old--{ds_old}'
        raise AssertionError(assertion_message)

params_modeling = [(cdf_file, input_data_folder_modeling, test_data_path_modeling) for cdf_file in cdf_files_modeling]

@pytest.mark.parametrize('cdf_file,new_data_dir,test_dir', params_modeling)
def test_equals_cdf_modeling(data_tuple_cdf):
    ds_new, ds_old = data_tuple_cdf

    if not ds_new.equals(ds_old):
        assertion_message = f'New and old datasets are not equal \n\n--New--\n{ds_new}\n\n--Old--{ds_old}'
        raise AssertionError(assertion_message)


params_final = [(cdf_file, input_data_folder_final, test_data_path_final) for cdf_file in cdf_files_final]

@pytest.mark.parametrize('cdf_file,new_data_dir,test_dir', params_final)
def test_equals_cdf_final(data_tuple_cdf):
    ds_new, ds_old = data_tuple_cdf

    if not ds_new.equals(ds_old):
        assertion_message = f'New and old datasets are not equal \n\n--New--\n{ds_new}\n\n--Old--{ds_old}'
        raise AssertionError(assertion_message)

sim_input_file = 'sim_input_mean.xlsx'
sim_input_file_new = pjoin(REPO_DIR, 'final', 'sim_input', 'output', sim_input_file)
sim_input_file_test = pjoin(test_data_path_final, sim_input_file)

def test_sim_input_mean():
    ds_new = pd.read_excel(sim_input_file_new)
    ds_old = pd.read_excel(sim_input_file_test)

    if not ds_new.equals(ds_old):
        assertion_message = f'New and old datasets are not equal \n\n--New--\n{ds_new}\n\n--Old--{ds_old}'
        raise AssertionError(assertion_message)
