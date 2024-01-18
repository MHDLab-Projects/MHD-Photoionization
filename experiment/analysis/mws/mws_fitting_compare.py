
#%%[markdown]

# # Microwave Scattering fit comparison

#  TODO: one date is offset in time
#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis.mws.fitting import gen_model_dnedt

dss_p = []

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

ds_lecroy = ds_lecroy.drop(0,'kwt')

#Make a general mws method? 
counts = ds_lecroy['AS'].isel(time=0).count('mnum')
ds_lecroy = ds_lecroy.where(counts > 100)

# %%

da_sel = ds_lecroy['AS']#.sel(kwt=1, method='nearest')

#%%[markdown]

# ## fitting after averaging mnum

# %%

da_fit = da_sel.mean('mnum')

da_fit.plot(hue='run_plot', row='kwt', x='time')
plt.yscale('log')

plt.xlim(-1,40)


#%%

# pre_norm_cutoff=5e-4
# da_fit = da_fit.where(da_fit.mws._pulse_max() > pre_norm_cutoff) # Targeting low power...
da_fit = da_fit.mws.fit_prep(pre_norm_cutoff=5e-4, remove_negative=False, min_mnum=None)

mod, params = gen_model_dnedt(take_log=False)

params['ne0'].value = Quantity(2e12, 'cm**-3').to('um**-3').magnitude #TODO: value from cfd
params['ne0'].vary = False
params['kr'].vary = True

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(mod, params)

#%%


ds_mws_fit.to_array('var').plot(row='kwt',hue='var', col='run')

plt.yscale('log')



#%%

ds_p2 = xr.merge([
    ds_p['kr'].rename('mean'),
    ds_p_stderr['kr'].rename('std')
    ])


dss_p.append(ds_p2.assign_coords(method='avgmnum_1'))

#%%
ds_p


#%%[markdown]

# ## Global optimization

# TODO: This cant handle nans when taking log. Just taking one test case for now. 

#%%

da_fit = da_sel.copy()

da_fit = da_fit.unstack('run')

da_fit = da_fit.mws.fit_prep(pre_norm_cutoff=5e-4, remove_negative=False)

mod, params = gen_model_dnedt(take_log=False)

params['ne0'].value = Quantity(2e12, 'cm**-3').to('um**-3').magnitude #TODO: value from cfd
params['ne0'].vary = False
params['kr'].vary = True

fits, ds_p, ds_p_stderr = da_fit.mws.perform_fit(mod, params, method='global')

ds_mws_fit = xr.merge([fits, da_fit])

ds_mws_fit = ds_mws_fit.xr_utils.stack_run()
ds_p = ds_p.xr_utils.stack_run()
ds_p_stderr = ds_p_stderr.xr_utils.stack_run()


#%%

ds_mws_fit.mean('mnum').to_array('var').plot(row='kwt',hue='var', col='run')

plt.yscale('log')

#%%

# da_fit.xr_utils.stack_run().mean('mnum').plot(hue='run_plot', row='kwt', x='time')

#%%

# fig, axes = plt.subplots(len(ds_mws_fit.indexes['run']), figsize=(5,15))

# for i, (group, ds) in enumerate(ds_mws_fit.groupby('run')):
#     ax=axes[i]

#     ds['AS'].isel(mnum=[0,10,20,30,40,50]).plot(hue='mnum', ax=ax)

#     # plt.plot(xs_eval, y_fit)

#     ds_mws_fit['AS_fit'].sel(run=group).plot(ax=ax, color='black')

#     ax.set_yscale('log')
#     ax.get_legend().remove()
#     ax.set_ylim(1e-3,)

#%%

ds_p2 = xr.merge([
    ds_p['kr'].rename('mean'),
    ds_p_stderr['kr'].rename('std')
    ])

dss_p.append(ds_p2.assign_coords(method='pipe_2'))

#%%[markdown]

# ## exponential fit

#%%

from mhdpy.analysis.mws.fitting import gen_model_exp

mod, params = gen_model_exp()

da_fit = da_sel.mean('mnum')

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(
    mod, params,
    fit_timewindow=slice(Quantity(5, 'us'),Quantity(20, 'us')) 
    )

ds_mws_fit = xr.merge([ds_mws_fit, da_fit.rename('AS_all')])

#%%
ds_mws_fit.to_array('var').plot(col='run', row='kwt', hue='var')

plt.yscale('log')

#%%

ne0 = Quantity(2e12, 'cm**-3').to('um**-3').magnitude

# tau = 1/(2*kr*ne0)
# kr = 1/(2*tau*ne0)

ds_p['kr'] = 1/(2*ne0*ds_p['decay'])

ds_p['kr']

#%%

ds_p2 = xr.merge([
    ds_p['kr'].rename('mean'),
    ds_p['kr'].where(False).fillna(0).rename('std') # kr std not defined yet for exponential fit...
    ])

dss_p.append(ds_p2.assign_coords(method='exp'))

#%%

ds_p = xr.concat(dss_p, dim='method')

ds_p

#%%

# Define your custom plotting function
def custom_plot(x, mean, std, **kwargs):
    x = [str(t) for t in x]
    plt.errorbar(x, mean, yerr=std, capsize=4, marker='o')


fg = xr.plot.FacetGrid(ds_p, row='run', col='method')

# Apply the custom function to each facet

fg.map(custom_plot, 'kwt', 'mean', 'std')


#%%

# ds_p_plot = ds_p.unstack('run').stack

fig, axes = plt.subplots(len(ds_p.indexes['run']), figsize=(5,15))

for i, (method, ds) in enumerate(ds_p.groupby('run')):

    ax=axes[i]

    ds['mean'].plot(hue='method', ax=ax, marker='o')

# %%
ds