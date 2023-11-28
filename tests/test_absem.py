
from pprint import pprint
import pytest
import warnings

from mhdpy.analysis.standard_import import *

test_data_path = 'tests/test_data' 
input_data_folder = pjoin(REPO_DIR, 'experiment','data','proc_data')

@pytest.fixture(scope='module')
def data_tuple(request) -> (xr.Dataset, xr.Dataset):

    fp_old = pjoin(test_data_path, 'ds_absem.cdf')
    fp_new = pjoin(input_data_folder, 'ds_absem.cdf')

    ds_old = xr.load_dataset(fp_old)
    ds_new = xr.load_dataset(fp_new)

    return (ds_new, ds_old) 

def test_equals(data_tuple):
    ds_new, ds_old = data_tuple

    assert ds_new.equals(ds_old)
