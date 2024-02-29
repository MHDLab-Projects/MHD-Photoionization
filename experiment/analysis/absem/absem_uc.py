#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.xr_utils import WeightedMeanAccessor

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray


#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem

#%%

from mhdpy.analysis.absem import calc_alpha_scipp

ds = calc_alpha_scipp(ds_absem)

ds = ds.xr_utils.stack_run()

ds

#%%

from mhdpy.plot.common import xr_errorbar_axes


fig, axes = plt.subplots(len(ds.coords['run']),2 , figsize=(5,10), sharex=True, sharey=True)

for i, run in enumerate(ds.coords['run']):
    ds_sel = ds.sel(run=run).sel(mp='barrel').dropna('kwt',how='all')

    #TODO: move in ds_sel
    ds_sel = ds_sel.drop([c for c in ds_sel.coords if c not in ds_sel.indexes])

    ax = axes[i, 0]

    if len(ds_sel.coords['kwt']):
        xr_errorbar_axes(ds_sel['mean'], ds_sel['std'], ax, huedim='kwt')

        ax.get_legend().remove()

    ds_sel = ds.sel(run=run).sel(mp='mw_horns').dropna('kwt',how= 'all')

    #TODO: move in ds_sel
    ds_sel = ds_sel.drop([c for c in ds_sel.coords if c not in ds_sel.indexes])

    ax = axes[i, 1]

    if len(ds_sel.coords['kwt']):
        xr_errorbar_axes(ds_sel['mean'], ds_sel['std'], ax, huedim='kwt')
        ax.get_legend().remove()

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

# %%


ds = ds.assign_coords(run_plot = ('run', ds.indexes['run'].values))

ds_sel = ds.sel(mp='barrel').dropna('kwt',how='all')

ds_sel['mean'].plot(col='run', row='kwt')

plt.ylim(-0.1,1.1)


#%%

ds_sel = ds.sel(mp='barrel').dropna('kwt',how='all')

# Define your custom plotting function
def custom_plot(wavelength, mean, std, **kwargs):
    lower_bound = mean - std
    upper_bound = mean + std
    plt.plot(wavelength, mean, **kwargs)
    plt.fill_between(wavelength, lower_bound, upper_bound, alpha=0.5, color='gray')  # alpha sets the transparency

# Create a FacetGrid
g = xr.plot.FacetGrid(ds_sel, col='run', row='kwt')

# Apply the custom function to each facet
g.map(custom_plot, 'wavelength', 'mean', 'std')

plt.ylim(-0.1,1.1)

#%%

# Create a FacetGrid
g = xr.plot.FacetGrid(ds_sel, col='run', row='kwt')

# Apply the custom function to each facet
g.map(custom_plot, 'wavelength', 'mean', 'std')

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

#%%

#%%

# Define your custom plotting function
def custom_plot(wavelength, mean, std, count, **kwargs):
    lower_bound = mean - std
    upper_bound = mean + std
    plt.plot(wavelength, mean, **kwargs)
    plt.fill_between(wavelength, lower_bound, upper_bound, alpha=0.5, color='gray')  # alpha sets the transparency
    # count should be a single value, add it as text to the plot in the upper left
    if not np.isnan(count):
        plt.text(0.05, 0.9, '{}'.format(int(count)), transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')



# Create a FacetGrid
g = xr.plot.FacetGrid(ds_sel, col='run', row='kwt')

# Apply the custom function to each facet
g.map(custom_plot, 'wavelength', 'mean', 'std', 'count')

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1 

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds.rename(mean='alpha')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_1(ds_fit, spect_red_dict)

#%%

ds_nK = xr.merge([
    ds_p['nK_m3'].to_dataset(name='mean'),
    ds_p_stderr['nK_m3'].to_dataset(name='stderr'),
    ds['count']  #TODO: better naming
])

ds_nK['std'] = ds_nK['stderr']*np.sqrt(ds_nK['count'])

ds_nK

#%%
# da.values

g = xr.plot.FacetGrid(ds_nK, col='mp', row='run')

g.map(plt.errorbar, 'kwt', 'mean', 'std', capsize=5)

plt.yscale('log')
#%%


from mhdpy.plot.common import xr_errorbar_axes

ds_sel = ds_nK.sel(mp='barrel').dropna('kwt',how='all').drop('mp')

fig, axes = plt.subplots()

xr_errorbar_axes(ds_sel['mean'], ds_sel['std'], axes, huedim='run')

plt.yscale('log')
plt.xscale('log')
#%%

ds_nK2 = ds_nK.wma.calc_weighted_mean('run')

fig, axes = plt.subplots()

xr_errorbar_axes(ds_nK2['mean'], ds_nK2['std'], axes, huedim='mp')

plt.yscale('log')
plt.xscale('log')


#%%

# Compare to simply calculating new mean/std over run 
ds_nK2 = ds_nK.wma.initialize_stat_dataset('mean', 'run')

fig, axes = plt.subplots()

xr_errorbar_axes(ds_nK2['mean'], ds_nK2['std'], axes, huedim='mp')

plt.yscale('log')
plt.xscale('log')



