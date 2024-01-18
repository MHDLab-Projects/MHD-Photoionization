
#%%[markdown]

# # Microwave Scattering fit comparison

#  TODO: one date is offset in time
#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

dss_p = []

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
# pre_norm_cutoff=5e-4
# da_fit = da_fit.where(da_fit.mws._pulse_max() > pre_norm_cutoff) # Targeting low power...
da_fit = da_fit.mws.fit_prep(pre_norm_cutoff=5e-4, remove_negative=False)

mod, params = gen_model_dnedt(take_log=False)

params['ne0'].value = Quantity(2e12, 'cm**-3').to('um**-3').magnitude #TODO: value from cfd
params['ne0'].vary = False
params['kr'].vary = True

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(mod, params)

#%%

da_fit
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

from mhdpy.analysis.mws.fitting import gen_model_dnedt, pipe_fit_mws_2

da_fit = da_sel.copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)



#%%

ds_mws_fit


#%%

fig, axes = plt.subplots(len(ds_mws_fit.indexes['run']), figsize=(5,15))

for i, (group, ds) in enumerate(ds_mws_fit.groupby('run')):
    ax=axes[i]

    ds['AS'].isel(mnum=[0,10,20,30,40,50]).plot(hue='mnum', ax=ax)

    # plt.plot(xs_eval, y_fit)

    ds_mws_fit['AS_fit'].sel(run=group).plot(ax=ax, color='black')

    ax.set_yscale('log')
    ax.get_legend().remove()
    ax.set_ylim(1e-3,)

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
ds_mws_fit.to_array('var').plot(row='run',hue='var')

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


fg = xr.plot.FacetGrid(ds_p, row='method')

# Apply the custom function to each facet

fg.map(custom_plot, 'run', 'mean', 'std')


#%%

for method, ds in ds_p.groupby('method'):
    ds = ds.squeeze() # having to do this 
    plt.errorbar(ds['run_plot'], ds['mean'], yerr=ds['std'], capsize=4, marker='o', label=method)

plt.legend()
