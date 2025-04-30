#%%
from mhdlab.analysis.standard_import import *

from mhdlab.analysis.aas.fitting import gen_model_alpha_blurred, pipe_fit_alpha_1, pipe_fit_alpha_2
from mhdlab.analysis.aas.fit_prep import interp_alpha
from mhdlab.analysis import aas
from mhdlab.xr_utils import fit_da_lmfit

from pi_paper_utils.constants import *

import pytest    

@pytest.fixture
def ds_aas():
    tc = '53x'
    ds_aas = xr.load_dataset(pjoin(REPO_DIR, 'tests', 'test_data', 'aas', '{}.cdf'.format(tc)))

    ds_aas = ds_aas.sel(mp='barrel').sel(date='2023-05-24').sel(run_num= 1)
    ds_aas = ds_aas.dropna('mnum','all').mean('mnum')
    ds_aas = ds_aas.aas.calc_alpha()
    return ds_aas

@pytest.fixture
def ds_aas_all_mnum():
    tc = '53x'
    ds_aas = xr.load_dataset(pjoin(REPO_DIR, 'tests', 'test_data', 'aas', '{}.cdf'.format(tc)))

    ds_aas = ds_aas.sel(mp='barrel').sel(date='2023-05-24').sel(run_num= 1)
    ds_aas = ds_aas.aas.calc_alpha()
    return ds_aas

def test_calc_alpha(ds_aas):
    ds_aas = ds_aas.aas.calc_alpha()

    pass

def test_aas_fit_interpolated(ds_aas):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_aas['alpha'].sel(kwt=1,method='nearest')

    da_alpha = da_alpha.sel(wavelength=slice(730,790))
    da_alpha = interp_alpha(da_alpha)

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=True)

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.32787e+22, rel=1e-4)


def test_aas_fit_uninterpolated(ds_aas):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_aas['alpha'].sel(kwt=1,method='nearest')

    da_alpha = da_alpha.sel(wavelength=slice(730,790))

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)

    assert ds_p['nK_m3'].item() == pytest.approx(1.2287e+22, rel=1e-4)

def test_pipe_fit_alpha_1(ds_aas):
    spectral_reduction_params_fp = os.path.join(REPO_DIR, 'experiment', 'metadata', 'spectral_reduction_params.csv')
    spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

    ds_aas = ds_aas.sel(kwt=1,method='nearest')
    ds_alpha_fit, ds_p, ds_p_stderr = aas.fitting.pipe_fit_alpha_1(ds_aas, fit_prep_kwargs={'spect_red_dict': spect_red_dict})

    assert ds_p['nK_m3'].item() == pytest.approx(1.17677e+22, rel=1e-4)

from mhdlab.analysis.aas.fitting import pipe_fit_alpha_num_1
import pi_paper_utils as ppu
def test_pipe_fit_alpha_num_1(ds_aas):
    ds_aas = ds_aas.drop(0, 'kwt')
    ds_cfd_beam = ppu.fileio.load_cfd_beam(kwt_interp=ds_aas.coords['kwt'].values)

    ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0).sel(motor=180,method='nearest').sel(mp='barrel')

    da_cfd_beam = ds_cfd_beam[CFD_K_SPECIES_NAME]
    da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

    ds_fit_aas = ds_aas.sel(kwt= da_cfd_beam.kwt.values) #TODO: downselecting as we don't have cfd for all kwt. Remove once we do
    # da_cfd_beam = da_cfd_beam.dropna('dist', how='all') # Having to do this with coarse CFD data?
    ds_fit_aas, ds_p_aas, ds_p_stderr_aas = pipe_fit_alpha_num_1(ds_fit_aas, perform_fit_kwargs={'nK_profile':da_cfd_beam})



def test_pipe_fit_alpha_2(ds_aas_all_mnum):
    ds_alpha_fit, ds_p, ds_p_stderr = aas.fitting.pipe_fit_alpha_2(ds_aas_all_mnum)

def test_wing_cut(ds_aas):
    ds_aas = ds_aas.sel(kwt=1,method='nearest') # Wing cut needs just wavelength dim
    ds_aas = ds_aas.aas.reduce_keep_wings()

from mhdlab.xr_utils import fit_da_lmfit_global
def test_aas_fit_global(ds_aas_all_mnum):
    """Tests a fit of a single absorption emssision spectrum and compares the fitted nK value to a previously obtained value"""

    da_alpha = ds_aas_all_mnum['alpha']

    da_alpha = da_alpha.sel(wavelength=slice(730,790))

    final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

    da_fit = da_alpha.pint.dequantify()
    wls = da_alpha.coords['wavelength'].values
    fits, ds_p, ds_p_stderr = fit_da_lmfit_global(da_fit, final_model, pars, 'wavelength', wls)

    # assert ds_p['nK_m3'].item() == pytest.approx(1.2287e+22, rel=1e-4