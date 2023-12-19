#%%
from mhdpy.analysis.standard_import import *

from mhdpy.analysis.spectral import interp_alpha
from mhdpy.analysis.spectral import gen_model_alpha_blurred
from mhdpy.analysis.xr import fit_da_lmfit

import pytest    

@pytest.fixture
def da_alpha():
    
    tc = '53x'

    ds_absem = xr.load_dataset(pjoin(REPO_DIR, 'tests', 'test_data', 'absem', '{}.cdf'.format(tc)))
    # Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

    ds_absem['diff'] = ds_absem['led_on'] - ds_absem['led_off']
    ds_absem['alpha'] = 1 - ds_absem['diff']/ds_absem['calib']

    ds_sel = ds_absem.sel(mp='barrel').sel(date='2023-05-24').sel(run_num= 1).sel(kwt=1,method='nearest')
    da_alpha = ds_sel['alpha'].dropna('mnum','all').mean('mnum')

    return da_alpha

def test_absem_fit_interpolated(da_alpha):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = da_alpha.sel(wavelength=slice(730,790))
    da_alpha = interp_alpha(da_alpha)

    final_model, pars = gen_model_alpha_blurred()

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.32787e+22, rel=1e-4)


def test_absem_fit_uninterpolated(da_alpha):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = da_alpha.sel(wavelength=slice(730,790))

    final_model, pars = gen_model_alpha_blurred()

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.2287e+22, rel=1e-4)