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
