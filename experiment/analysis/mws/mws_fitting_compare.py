
#%%[markdown]

# # Microwave Scattering fit comparison

#  TODO: one date is offset in time
#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.mws import calc_mag_phase_AS

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num']).dropna('run',how='all')
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = calc_mag_phase_AS(ds_lecroy)

# %%

da_sel = ds_lecroy['AS'].sel(kwt=1, method='nearest')

#%%[markdown]

# ## fitting after averaging mnum

#%%

da_fit = da_sel.mean('mnum')

tc_dim = 'run'
pulse_max = da_fit.sel(time=slice(-1,1)).max('time')
da_fit = da_fit.where(pulse_max > 5e-4) # Targeting low power...
da_fit = da_fit.dropna(tc_dim, how='all')

da_fit = da_fit/pulse_max


# %%


da_fit.plot(hue='run_plot', x='time')
plt.yscale('log')

plt.xlim(-1,40)


#%%

from mhdpy.analysis.mws.fitting import fit_fn 
from mhdpy.xr_utils import fit_da_lmfit
from lmfit import Model


mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0, take_log=True))

params = mod.make_params()
params['dne'].value = 1e14
params['ne0'].value = 3e12
params['dne'].min = 0
params['ne0'].min = 0

#%%

da_fit = np.log(da_fit).dropna('time', 'all')

da_fit_region = da_fit.sel(time=slice(0,30))
x_eval = da_fit.sel(time=slice(0,30)).coords['time'].values

# da_fit_coarse = da_fit_region.coarsen(time=10).mean()
fits, ds_p, ds_p_stderr = fit_da_lmfit(da_fit_region, mod, params, 'time', x_eval)

fits = np.exp(fits)
da_fit = np.exp(da_fit)
#%%


ds = xr.merge([fits, da_fit])

da = ds.to_array('var')

# for run in ds.coords['run_plot']:
#     da_sel = da.sel(run=run)
g= da.plot( x='time',hue='var', row='run')

# plt.xscale('log')
plt.yscale('log')


#%%

da_mean = ds_p['ne0']
da_mean.name = 'mean'
da_std = ds_p_stderr['ne0']
da_std.name='std'

df = xr.merge([da_mean, da_std]).to_dataframe()[['mean','std']]

df

#%%[markdown]

# ## Global optimization

# TODO: This cant handle nans when taking log. Just taking one rest case for now. 


#%%

da_fit = da_sel.copy()

da_fit = da_fit.isel(run=-2).dropna('mnum', how='all')

pulse_max = da_fit.sel(time=slice(-1,1)).max('time')
da_fit = da_fit/pulse_max

#%%

# da_fit.isel(mnum=[0,10,20,30,40,50]).plot(hue='mnum')

da_fit#.dropna('mnum', how='all')


#%%


from lmfit import minimize, Parameters, report_fit, Model


mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0, take_log=False))

def objective(params, x, data):
    """Calculate total residual for fits of blurredalpha_2peak to several data sets."""
    mnums = data.shape[0]
    resid = 0.0*data[:]


    param_dict=  params.valuesdict()

    model_vals = fit_fn(x, **param_dict)


    # make residual per data set
    for i in range(mnums):
        resid[i, :] = data[i,:] - model_vals


    # now flatten this to a 1D array, as minimize() needs

    return resid.flatten()

# da_fit2 = np.log(da_fit).dropna('time', 'all')
da_fit2 = da_fit

da_fit_region = da_fit2.sel(time=slice(0,30))
x_eval = da_fit2.sel(time=slice(0,30)).coords['time'].values

data_vals = da_fit_region.values.T

# Just using this to get defaul parameters

out = minimize(objective, params, args=(x_eval, data_vals))
report_fit(out)

#%%

param_dict=  out.params.valuesdict()
y_fit = fit_fn(x_eval, **param_dict)

#%%

da_fit_region.mean('mnum').plot()

plt.plot(x_eval, y_fit)

plt.yscale('log')

