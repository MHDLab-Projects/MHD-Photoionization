#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.mws_utils import calc_mag_phase_AS
from mhdpy.analysis.xr import WeightedMeanAccessor

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem = ds_absem.stack(run = ['date','run_num']).dropna('run', how='all')

ds_absem['diff'] = ds_absem['led_on'] - ds_absem['led_off']
ds_absem['alpha'] = 1 - ds_absem['diff']/ds_absem['calib']

ds_absem
# %%

ds_sel = ds_absem.sel(mp='barrel').sel(run=('2023-05-24', 1)).sel(kwt=1,method='nearest')

da_sel = ds_sel['alpha'].dropna('mnum','all')

da_sel

#%%

da_sel.plot(hue='mnum')

plt.ylim(-0.1,1.1)

# add horizontal line at 0 
plt.axhline(0, color='k', linestyle='--', linewidth=0.5)

#%%

from mhdpy.analysis.spectral import interp_alpha


da_alpha = interp_alpha(da_sel)


#%%

from mhdpy.analysis.spectral import model_blurredalpha_2peak, alpha_cut

final_model, pars = model_blurredalpha_2peak()


spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()


alpha_tc = da_alpha

# beta = -np.log(1-da_alpha)/pars['L'].value
# beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
# alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

# peak1_low, peak1_high = spect_red_dict['peak1_low'], spect_red_dict['peak1_high']
# peak2_low, peak2_high = spect_red_dict['peak2_low'], spect_red_dict['peak2_high']
# wls = alpha_tc.sel(wavelength=slice(peak1_low, peak1_high)).coords['wavelength']
# alpha_tc = alpha_tc.drop_sel(wavelength=wls)

# wls = alpha_tc.sel(wavelength=slice(peak2_low, peak2_high)).coords['wavelength']
# alpha_tc = alpha_tc.drop_sel(wavelength=wls)    

alpha_tc_red = alpha_tc.dropna('wavelength', how='all')

alpha_tc_red = alpha_tc_red.sel(wavelength=slice(750,790))

# %%

# spect_red_dict

# alpha_tc_red.sel(wavelength=slice(peak1_low,peak1_high))

# alpha_tc_red.dropna('wavelength', how='any')

#%%

alpha_tc_red.plot(hue='mnum')

plt.ylim(-0.1,1.1)
# plt.xlim(760,775)

plt.gca().get_legend().remove()


#%%

from lmfit import minimize, Parameters, report_fit, Model
from mhdpy.analysis.spectral import blurred_alpha

def objective(params, x, data):
    """Calculate total residual for fits of blurredalpha_2peak to several data sets."""
    mnums = data.shape[0]
    resid = 0.0*data[:]


    param_dict=  params.valuesdict()
    param_dict = {k: v for k, v in param_dict.items() if k != 'nK_m3'}

    model_vals = blurred_alpha(x, **param_dict)


    # make residual per data set
    for i in range(mnums):
        resid[i, :] = data[i,:] - model_vals


    # now flatten this to a 1D array, as minimize() needs

    return resid.flatten()

# ... (your data generation code here) ...

x = alpha_tc_red.coords['wavelength'].values
data_vals = alpha_tc_red.values.T


final_model, pars = model_blurredalpha_2peak()
# add parameters for your model here

out = minimize(objective, pars, args=(x, data_vals))
report_fit(out)

#%%

plt.figure()

param_dict=  out.params.valuesdict()
param_dict = {k: v for k, v in param_dict.items() if k != 'nK_m3'}


for i in range(data_vals.shape[0]):
    y_fit = blurred_alpha(x, **param_dict)
    plt.plot(x, data_vals[i, :], 'o', x, y_fit, '-')



# %%[markdown]

# compare to fits of each individual spectrum

#%%

from lmfit import minimize, Parameters, report_fit, Model

da_sel = ds_sel['alpha'].dropna('mnum','all')

da_sel = interp_alpha(da_sel)

da_sel = da_sel.sel(wavelength=slice(750,790))

da_sel.plot(hue='mnum')

#%%

from mhdpy.analysis.xr import fit_da_lmfit

final_model, pars = model_blurredalpha_2peak()

alpha_tc = da_sel

wls = alpha_tc.coords['wavelength'].values
fits, ds_p, ds_p_stderr = fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

#%%

da = xr.merge([alpha_tc, fits]).to_array('var')

g = da.plot(hue='var', row='mnum')

for i , ax in enumerate(g.axes.flatten()):
    # ax.set_ylim(-0.1,1.1)
    # ax.set_xlim(750,790)
    mnum = da.coords['mnum'].values[i]
    nK_m3 = ds_p.sel(mnum=mnum)['nK_m3'].item()
    nK_m3_stderr = ds_p_stderr.sel(mnum=mnum)['nK_m3'].item()

    ax.set_title('mnum = {}, nK_m3 = {:.2e} +/- {:.2e}'.format(mnum, nK_m3, nK_m3_stderr))

#%%

from mhdpy.plot.common import xr_errorbar

xr_errorbar(ds_p['nK_m3'], ds_p_stderr['nK_m3'])

#%%
nK_mean = ds_p['nK_m3'].mean('mnum').item()
nK_std = ds_p['nK_m3'].std('mnum').item()

print('nK_mean = {:.2e}, nK_std = {:.2e}'.format(nK_mean, nK_std))

#%%[markdown]

# Compare to average before alpha calculation

#%%

led_on = ds_sel['led_on'].mean('mnum')
led_off = ds_sel['led_off'].mean('mnum')

alpha = 1 - (led_on - led_off)/ds_sel['calib'].mean('mnum')

alpha = interp_alpha(alpha)

alpha = alpha.sel(wavelength=slice(750,790))

alpha.plot()

plt.ylim(-0.1,1.1)


#%%

alpha_tc_red = alpha

wls = alpha_tc_red.coords['wavelength'].values
fits, ds_p, ds_p_stderr = fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'
alpha.name = 'data'

#%%

da = xr.merge([alpha, fits]).to_array('var')

da.plot(hue='var')

#%%

print("nK_m3 = {:.2e} +/- {:.2e}".format(ds_p['nK_m3'].item(), ds_p_stderr['nK_m3'].item()))
