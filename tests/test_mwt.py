#%%
from mhdlab.analysis.standard_import import *

from mhdlab.analysis.aas.fitting import gen_model_alpha_blurred, pipe_fit_alpha_1
from mhdlab.analysis.aas.fit_prep import interp_alpha
from mhdlab.analysis import mwt
from mhdlab.xr_utils import fit_da_lmfit
from pint import Quantity

import pytest    

@pytest.fixture
def ds_mwt():
    DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

    tc = '53x'
    ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

    seldict = dict(date='2023-05-24', run_num=1)
    ds_lecroy = ds_lecroy.sel(seldict)

    ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected


    return ds_lecroy

from mhdlab import xr_utils

@pytest.fixture
def ds_mwt_all_mnum():

    DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

    tc = '53x'
    ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

    seldict = dict(date='2023-05-24', run_num=1)
    ds_lecroy = ds_lecroy.sel(seldict)
    # ds_lecroy = ds_lecroy.groupby('kwt').apply(lambda x: x.xr_utils.assign_mnum('mnum'))

    ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
    ds_lecroy = ds_lecroy.mwt.calc_AS_rel()

    return ds_lecroy


import pi_paper_utils as ppu
from pint import Unit
def test_calc_AS_rel(ds_mwt):
    ds_mwt = ds_mwt.mwt.calc_AS_rel()

    assert 'dAS_rel' in ds_mwt
    
    assert ds_mwt['dAS_rel'].pint.quantify().pint.units == Unit('dimensionless')

def test_calc_mag_phase_AS_abs(ds_mwt):
    ds_nothing = ppu.fileio.load_mwt_T0()
    ds_mwt = ds_mwt.mwt.calc_AS_abs(mag_0=ds_nothing)

    assert 'AS_abs' in ds_mwt

    assert ds_mwt['AS_abs'].pint.quantify().pint.units == Unit('dimensionless')


def test_pipe_fit_mwt_1_avgbef(ds_mwt):

    ds_mwt = ds_mwt.mean('mnum') 
    ds_mwt = ds_mwt.mwt.calc_AS_rel()

    ds_mwt_fit, ds_p, ds_p_stderr = mwt.fitting.pipe_fit_mwt_1(ds_mwt['dAS_rel'])

    assert 'AS_fit' in ds_mwt_fit

    ne0_fit_values = ds_p['ne0'].values

    #replace nan with 0 
    ne0_fit_values = np.nan_to_num(ne0_fit_values)

    # expected = [           0, 2.56266691e+12, 1.07160749e+13, 6.17915301e+12,
    #    1.13125377e+13, 1.19845333e+13, 1.00853221e+13, 9.72104615e+12]

    # assert np.allclose(ne0_fit_values, expected, rtol=1e-2)

#TODO: this fails when taking average after calculating mnum. Need to investigate why. 
def test_pipe_fit_mwt_1_avgafter(ds_mwt):

    ds_mwt = ds_mwt.mwt.calc_AS_rel()
    ds_mwt = ds_mwt.mean('mnum') 

    ds_mwt_fit, ds_p, ds_p_stderr = mwt.fitting.pipe_fit_mwt_1(ds_mwt['dAS_rel'])

    assert 'AS_fit' in ds_mwt_fit


def test_pipe_fit_mwt_2(ds_mwt_all_mnum):
    ds_mwt_fit, ds_p, ds_p_stderr = mwt.fitting.pipe_fit_mwt_2(ds_mwt_all_mnum['dAS_rel'].mean('mnum'))


def test_pipe_fit_mwt_3(ds_mwt):
    ds_mwt = ds_mwt.mean('mnum') 
    ds_mwt = ds_mwt.mwt.calc_AS_rel()

    ds_mwt_fit, ds_p, ds_p_stderr = mwt.fitting.pipe_fit_mwt_3(ds_mwt['dAS_rel'])

# from mhdlab.analysis.mwt.fitting import gen_model_dnedt
# from mhdlab.xr_utils import fit_da_lmfit_global
# def test_mwt_fit_global_ne0(ds_mwt_all_mnum):

#     da_fit = ds_mwt_all_mnum['dAS_rel']
#     da_fit = da_fit.sel(kwt=slice(0.1, 10)) #Low seed takes long time, but passes

#     mod, params = gen_model_dnedt(take_log=True)

#     ds_mwt_fit, ds_p, ds_p_stderr = da_fit.mwt.perform_fit(mod, params, method='global')

# def test_mwt_fit_global_kr(ds_mwt_all_mnum):

#     da_fit = ds_mwt_all_mnum['dAS_rel']
#     da_fit = da_fit.sel(kwt=slice(0.1, 10)) #Low seed takes long time, but passes

#     mod, params = gen_model_dnedt(take_log=True)

#     params['ne0'].value = Quantity(2e13, 'cm**-3').to('um**-3').magnitude
#     params['ne0'].vary = False
#     params['kr'].vary = True

#     ds_mwt_fit, ds_p, ds_p_stderr = da_fit.mwt.perform_fit(mod, params, method='global')