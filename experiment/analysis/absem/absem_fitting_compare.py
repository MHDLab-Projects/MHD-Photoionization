#%%[markdown]

# # Absorption emission fitting algorithm comparison


#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.absem.fit_prep import interp_alpha
from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred 
from mhdpy.analysis import absem

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem = ds_absem.stack(run = ['date','run_num']).dropna('run', how='all')

ds_absem = ds_absem.absem.calc_alpha()

ds_sel = ds_absem.sel(mp='barrel').sel(run=('2023-05-24', 1)).sel(kwt=1,method='nearest')

ds_sel = ds_sel.dropna('wavelength', how='all')
ds_sel = ds_sel.sel(wavelength=slice(750,790))

da_sel = ds_sel['alpha'].dropna('mnum','all')

da_sel

#%%

da_sel.plot(hue='mnum')

plt.ylim(-0.1,1.1)

# add horizontal line at 0 
plt.axhline(0, color='k', linestyle='--', linewidth=0.5)

#%%[markdown]

# ## Global Least squares reduction. 

# All measurement numbers are concatenated into one array, and fitted in one optimization process
# 


#%%

from mhdpy.analysis.absem.fitting import alpha_2peak
from mhdpy.xr_utils import fit_da_lmfit_global

da_alpha = interp_alpha(da_sel)

# Just using this to get defaul parameters
final_model, pars = gen_model_alpha_blurred()

x = da_alpha.coords['wavelength'].values
out = fit_da_lmfit_global(da_alpha, final_model, pars, 'wavelength', x)
out

#%%

plt.figure()

param_dict=  out.params.valuesdict()

da_alpha.plot(hue='mnum')

plt.plot(x, final_model.eval(out.params, x=x))


# %%[markdown]

# ## Fitting each individual spectrum


#%%

from mhdpy.xr_utils import fit_da_lmfit

final_model, pars = gen_model_alpha_blurred()

wls = da_alpha.coords['wavelength'].values
fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

#%%

da = xr.merge([da_alpha, fits]).to_array('var')

g = da.plot(hue='var', row='mnum')

for i , ax in enumerate(g.axes.flatten()):
    # ax.set_ylim(-0.1,1.1)
    # ax.set_xlim(750,790)
    mnum = da.coords['mnum'].values[i]
    nK_m3 = ds_p.sel(mnum=mnum)['nK_m3'].item()
    nK_m3_stderr = ds_p_stderr.sel(mnum=mnum)['nK_m3'].item()

    ax.set_title('mnum = {}, nK_m3 = {:.2e} +/- {:.2e}'.format(mnum, nK_m3, nK_m3_stderr))

#%%

ds_p_stderr['nK_m3']

#%%

plt.errorbar(ds_p.coords['mnum'] , ds_p['nK_m3'], ds_p_stderr['nK_m3'], fmt='o')

#%%
nK_mean = ds_p['nK_m3'].mean('mnum').item()
nK_std = ds_p['nK_m3'].std('mnum').item()

print('nK_mean = {:.2e}, nK_std = {:.2e}'.format(nK_mean, nK_std))

#%%[markdown]

# ## Compare to performing average before alpha calculation

#%%

led_on = ds_sel['led_on'].mean('mnum')
led_off = ds_sel['led_off'].mean('mnum')

da_alpha = 1 - (led_on - led_off)/ds_sel['calib'].mean('mnum')

da_alpha = interp_alpha(da_alpha)

da_alpha.plot()

plt.ylim(-0.1,1.1)


#%%


wls = da_alpha.coords['wavelength'].values
fits, ds_p, ds_p_stderr = fit_da_lmfit(da_alpha, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'
da_alpha.name = 'data'

#%%

da = xr.merge([da_alpha, fits]).to_array('var')

da.plot(hue='var')

#%%

print("nK_m3 = {:.2e} +/- {:.2e}".format(ds_p['nK_m3'].item(), ds_p_stderr['nK_m3'].item()))

# %%
