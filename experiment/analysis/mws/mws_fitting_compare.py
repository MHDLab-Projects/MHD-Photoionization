
#%%[markdown]

# # Microwave Scattering fit comparison

#  TODO: one date is offset in time
#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num']).dropna('run',how='all')
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

# %%

da_sel = ds_lecroy['AS'].sel(kwt=1, method='nearest')

#%%[markdown]

# ## fitting after averaging mnum

# %%

da_fit = da_sel.mean('mnum')

da_fit.plot(hue='run_plot', x='time')
plt.yscale('log')

plt.xlim(-1,40)


#%%

from mhdpy.analysis.mws.fitting import gen_model_dnedt, pipe_fit_mws_1 

pre_norm_cutoff = 5e-4
fit_timewindow = slice(0,30)
da_fit = da_fit.where(da_fit.mws._pulse_max() > pre_norm_cutoff) # Targeting low power...

mod, params = gen_model_dnedt()

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(mod, params, fit_timewindow=fit_timewindow)

#%%

da = ds_mws_fit.to_array('var')

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

df.plot(marker='o', yerr='std', capsize=5)

#%%


pre_norm_cutoff = 5e-4
fit_timewindow = slice(0,30)
da_fit = da_fit.where(da_fit.mws._pulse_max() > pre_norm_cutoff) # Targeting low power...

mod, params = gen_model_dnedt()

params['ne0'].value = 1e13 #TODO: value from cfd
params['ne0'].vary = False
params['kr'].vary = True

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(mod, params, fit_timewindow=fit_timewindow)
#%%

da_mean = ds_p['kr']
da_mean.name = 'mean'
da_std = ds_p_stderr['kr']
da_std.name='std'

df = xr.merge([da_mean, da_std]).to_dataframe()[['mean','std']]
df.plot(marker='o', yerr='std', capsize=5)
plt.ylabel('kr')


#%%[markdown]

# ## Global optimization

# TODO: This cant handle nans when taking log. Just taking one test case for now. 

#%%

da_fit = da_sel.copy()

# da_fit = da_fit.isel(run=-2).dropna('mnum', how='all')

da_fit = da_fit/da_fit.mws._pulse_max()

da_fit = da_fit.sel(time=slice(0,30))

#TODO: make fit_da_lmfit_global able to handle multindex
da_fit = da_fit.unstack('run')

from mhdpy.xr_utils import fit_da_lmfit_global

#TODO: having to take log 
mod, params = gen_model_dnedt(take_log=False)

params['ne0'].value = 1e13 #TODO: value from cfd
params['ne0'].vary = False
params['kr'].vary = True

xs_eval = da_fit.coords['time'].values
da_mws_fit, ds_p, ds_p_stderr = fit_da_lmfit_global(da_fit, mod, params, 'time', xs_eval)

da_mws_fit = da_mws_fit.stack(run=('date','run_num')).dropna('run',how='all')
ds_p = ds_p.stack(run=('date','run_num')).dropna('run',how='all')
ds_p_stderr = ds_p_stderr.stack(run=('date','run_num')).dropna('run',how='all')
da_fit = da_fit.stack(run=('date','run_num')).dropna('run',how='all')

#%%

fig, axes = plt.subplots(len(da_fit.indexes['run']), figsize=(5,10))

for i, (group, da) in enumerate(da_fit.groupby('run')):
    ax=axes[i]

    da.isel(mnum=[0,10,20,30,40,50]).plot(hue='mnum', ax=ax)

    # plt.plot(xs_eval, y_fit)

    da_mws_fit.sel(run=group).plot(ax=ax, color='black')

    ax.set_yscale('log')
    ax.get_legend().remove()

#%%

da_mean = ds_p['kr']
da_mean.name = 'mean'
da_std = ds_p_stderr['kr']
da_std.name='std'

df = xr.merge([da_mean, da_std]).to_dataframe()[['mean','std']]
df.plot(marker='o', yerr='std', capsize=5)
plt.ylabel('kr')