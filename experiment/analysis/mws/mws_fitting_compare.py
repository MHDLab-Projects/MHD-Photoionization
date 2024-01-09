
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
ds_lecroy = ds_lecroy.xr_utils.stack_run()

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

from mhdpy.analysis.mws.fitting import gen_model_dnedt, pipe_fit_mws_2

da_fit = da_sel.copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)

# ds_mws_fit = ds_mws_fit.mean('mnum')


#%%

ds_mws_fit


#%%

fig, axes = plt.subplots(len(ds_mws_fit.indexes['run']), figsize=(5,10))

for i, (group, ds) in enumerate(ds_mws_fit.groupby('run')):
    ax=axes[i]

    ds['AS'].isel(mnum=[0,10,20,30,40,50]).plot(hue='mnum', ax=ax)

    # plt.plot(xs_eval, y_fit)

    ds_mws_fit['AS_fit'].sel(run=group).plot(ax=ax, color='black')

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
# %%
