#%%
from mhdpy.analysis.standard_import import *

from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred, pipe_fit_alpha_1
from mhdpy.analysis.absem.fit_prep import interp_alpha
from mhdpy.analysis import absem
from mhdpy.xr_utils import fit_da_lmfit

import pytest    

@pytest.fixture
def ds_absem():
    tc = '53x'
    ds_absem = xr.load_dataset(pjoin(REPO_DIR, 'tests', 'test_data', 'absem', '{}.cdf'.format(tc)))

    ds_absem = ds_absem.sel(mp='barrel').sel(date='2023-05-24').sel(run_num= 1)
    ds_absem = ds_absem.dropna('mnum','all').mean('mnum')
    ds_absem = ds_absem.absem.calc_alpha()
    return ds_absem

@pytest.fixture
def ds_absem_all_mnum():
    tc = '53x'
    ds_absem = xr.load_dataset(pjoin(REPO_DIR, 'tests', 'test_data', 'absem', '{}.cdf'.format(tc)))

    ds_absem = ds_absem.sel(mp='barrel').sel(date='2023-05-24').sel(run_num= 1)
    ds_absem = ds_absem.absem.calc_alpha()
    return ds_absem

def test_absem_fit_interpolated(ds_absem):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_absem['alpha'].sel(kwt=1,method='nearest')

    da_alpha = da_alpha.sel(wavelength=slice(730,790))
    da_alpha = interp_alpha(da_alpha)

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=True)

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.32787e+22, rel=1e-4)


def test_absem_fit_uninterpolated(ds_absem):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_absem['alpha'].sel(kwt=1,method='nearest')

    da_alpha = da_alpha.sel(wavelength=slice(730,790))

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.2287e+22, rel=1e-4)

def test_pipe_fit_alpha_1(ds_absem):
    spectral_reduction_params_fp = os.path.join(REPO_DIR, 'experiment', 'metadata', 'spectral_reduction_params.csv')
    spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

    ds_absem = ds_absem.sel(kwt=1,method='nearest')
    ds_alpha_fit, ds_p, ds_p_stderr = absem.fitting.pipe_fit_alpha_1(ds_absem, spect_red_dict)

def test_wing_cut(ds_absem):
    ds_absem = ds_absem.sel(kwt=1,method='nearest') # Wing cut needs just wavelength dim
    ds_absem = ds_absem.absem.reduce_keep_wings()

from mhdpy.xr_utils import fit_da_lmfit_global
def test_absem_fit_global(ds_absem_all_mnum):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_absem_all_mnum['alpha']

    da_alpha = da_alpha.sel(wavelength=slice(730,790))

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit_global(da_alpha, final_model, pars, 'wavelength', wls)

    # assert ds_p['nK_m3'].item() == pytest.approx(1.2287e+22, rel=1e-4