#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.mws_utils import calc_mag_phase_AS

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

# Scipp does not write non indexed coordinates back to xarray well...
ds_absem = ds_absem.drop([c for c in ds_absem.coords if c not in ds_absem.indexes])

#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem

#%%


ds_mean = ds_absem.mean('mnum')
ds_std = ds_absem.std('mnum')

# %%

import scipp as sc


sds = sc.compat.from_xarray(ds_mean)
# sds['led_on']
sds['led_on'].variances = sc.compat.from_xarray(ds_std['led_on']).values
sds['led_off'].variances = sc.compat.from_xarray(ds_std['led_off']).values

sds

#%%

sds['diff'] = sds['led_on'] - sds['led_off']
sds['alpha'] = 1 - sds['diff']/sds['calib']
sds

#%%

da_mean = sc.compat.to_xarray(sds['alpha'])
da_mean.name = 'mean'
da_std = sc.compat.to_xarray(sc.stddevs(sds['alpha']))
da_std.name = 'std'

ds = xr.merge([da_mean, da_std])

ds = ds.stack(run = ['date','run_num']).dropna('run', how='all')

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

#%%

ds_sel.coords['kwt']
# %%


ds = ds.assign_coords(run_plot = ('run', ds.indexes['run'].values))

ds_sel = ds.sel(mp='barrel').dropna('kwt',how='all')

ds_sel['mean'].plot(col='run', row='kwt')

plt.ylim(-0.1,1.1)

#%%

da_sel = ds.to_array('stat')
da_sel

#%%
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

da_counts = ds_absem.count('mnum').mean('wavelength')['led_on']

da_counts = da_counts.stack(run = ['date','run_num']).dropna('run', how='all')
da_counts = da_counts.assign_coords(run_plot = ('run', da_counts.indexes['run'].values))
da_counts = da_counts.where(da_counts> 0, drop=True)

da_counts_sel = da_counts.sel(mp='barrel').dropna('kwt',how='all')

da_counts_sel

# da_counts_sel.plot(hue='run_plot', x='kwt', marker='o')


#%%

ds_sel['count'] = da_counts_sel

ds_sel

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

from mhdpy.analysis.spectral import perform_fit_alpha 

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds.rename(mean='alpha')

ds_p, ds_p_stderr, ds_alpha = perform_fit_alpha(ds_fit, spect_red_dict)

#%%
ds_nK = xr.merge([
    ds_p['nK_m3'].to_dataset(name='mean'),
    ds_p_stderr['nK_m3'].to_dataset(name='std'),
])

ds_nK

#%%

# testing uncertianties package
#TODO: move somewhere else, but refactor fitting first?

from uncertainties import unumpy as unp

ds_nK_unstacked = ds_nK.unstack('run')

arr = unp.uarray(ds_nK_unstacked['mean'].values, ds_nK_unstacked['std'].values)

# For some reason this appears to have the same order as arr
# dims = list(ds_nK_unstacked.indexes.keys())
dims = ds_nK_unstacked['mean'].dims #TODO: appears conversion to datarray is needed for correct dim order..

coords = {dim: ds_nK_unstacked.coords[dim].values for dim in dims}
da_nK_uc = xr.DataArray(arr, dims=dims, coords=coords)

da_nK_uc

#%%

da = da_nK_uc.sel(kwt=1, method='nearest').sel(mp='barrel').stack(run=['date','run_num'])

arr = da.values

np.nanmean(arr, axis=-1)

#%%

# TODO: can't make sense of the std dev

# narr = [v.nominal_value for v in arr]
# np.nanmean(narr)
# np.nanstd(narr)

narr = [v.std_dev**2 for v in arr]
np.sqrt(np.nanmean(narr))

#%%

da = da_nK_uc.stack(run=['date','run_num'])

da = xr.apply_ufunc(np.nanmean, da, input_core_dims=[['run']], output_core_dims=[[]], kwargs={'axis':-1})

da_nominal = xr.apply_ufunc(unp.nominal_values, da)
da_std = xr.apply_ufunc(unp.std_devs, da)

ds_nK_uc = xr.merge([da_nominal.rename('mean'), da_std.rename('std')])
ds_nK_uc    

#%%
# da.values

g = xr.plot.FacetGrid(ds_nK_uc, col='mp')

g.map(plt.errorbar, 'kwt', 'mean', 'std', capsize=5)

plt.yscale('log')
# %%
