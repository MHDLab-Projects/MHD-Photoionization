#%%
from mhdpy.analysis.standard_import import *

from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred, pipe_fit_alpha_1
from mhdpy.analysis.absem.fit_prep import interp_alpha
from mhdpy.analysis import mws
from mhdpy.xr_utils import fit_da_lmfit

import pytest    

@pytest.fixture
def ds_mws():
    DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

    tc = '53x'

    ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

    seldict = dict(date='2023-05-18', run_num=1)
    ds_lecroy = ds_lecroy.sel(seldict)

    ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
    ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

    ds_lecroy = ds_lecroy.mean('mnum')

    return ds_lecroy

def test_pipe_fit_mws_1(ds_mws):
    ds_mws_fit, ds_p, ds_p_stderr = mws.fitting.pipe_fit_mws_1(ds_mws['AS'])

    assert 'AS_fit' in ds_mws_fit

    ne0_fit_values = ds_p['ne0'].values

    #replace nan with 0 
    ne0_fit_values = np.nan_to_num(ne0_fit_values)

    expected = [           0, 2.56266691e+12, 1.07160749e+13, 6.17915301e+12,
       1.13125377e+13, 1.19845333e+13, 1.00853221e+13, 9.72104615e+12]

    assert np.allclose(ne0_fit_values, expected, rtol=1e-4)