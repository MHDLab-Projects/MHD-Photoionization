#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem

#%%[markdown]

# # Absem

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

# %%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='kwt', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit = ds_absem

ds_absem_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)


# %%
da_plot = ds_absem_fit[['alpha_red','alpha_fit']].to_array('var')

g = da_plot.sel(mp='barrel').plot(col='kwt',hue='var', row='run', ylim = (1e-3,2), yscale='log', figsize=(10,10))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')

#%%

da = ds_p['nK_m3']

# da = da.dropna('run', how='all')

g  = da.plot(hue='run_plot', x='kwt',col='mp', marker='o')

dropna(g)

plt.yscale('log')
plt.xscale('log')

#%%

ds_nK = xr.merge([
    ds_p['nK_m3'].to_dataset(name='mean'),
    ds_p_stderr['nK_m3'].to_dataset(name='stderr'),
    ds_absem.mean('wavelength').count('mnum')['alpha'].to_dataset(name='count')
])

ds_nK['std'] = ds_nK['stderr']*np.sqrt(ds_nK['count'])

# ds_nK2 = ds_nK.wma.calc_weighted_mean('run')
# ds_nK2 = ds_nK.wma.initialize_stat_dataset('mean', 'run')
# ds_nK2

#%%[markdown]

# # Lecroy

#%%

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

da_fit = ds_lecroy['AS']

# %%

g = da_fit.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

dropna(g)

#%%[markdown]

# ## Fixed kr, vary ne0, old pipeline

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_1 

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_1(da_fit.mean('mnum'))

#%%

da = ds_mws_fit.to_array('var')

# for run in ds.coords['run_plot']:
#     da_sel = da.sel(run=run)
g= da.plot( x='time',hue='var', row='kwt', col='run')

# plt.xscale('log')
plt.yscale('log')

dropna(g)

#%%

ds_ne0 = ds_p['ne0'].to_dataset(name='mean')
ds_ne0['std'] = ds_p_stderr['ne0']

plot.common.xr_errorbar(ds_ne0['mean'], ds_ne0['std'], huedim='run')

#%%[markdown]

# ## Fixed ne0, vary kr, new pipeline

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)

#TODO: sterr is nan where ds_p is not?
ds_p['kr'] = ds_p['kr'].where(~ds_p_stderr['kr'].isnull())
#%%

ds_kr
#%%
ds_kr = ds_p['kr'].to_dataset(name='mean')
ds_kr['std'] = ds_p_stderr['kr']

plot.common.xr_errorbar(ds_kr['mean'], ds_kr['std'], huedim='run')

plt.yscale('log')
# plt.ylim(1e-15,2e-14)


#%%[markdown]

# # Compare Lecroy and MWS


# %%

ds_lecroy['AS'].dropna('run', how='all').dropna('kwt', how='all')

#%%

da_stats = ds_lecroy['AS'].sel(time=slice(-1,1))

mws_max = da_stats.mean('mnum').max('time')

mws_std = da_stats.std('mnum').max('time')

#%%

ds_params = xr.merge([
    ds_nK['mean'].sel(mp='barrel').drop('mp').rename('nK_m3'),
    ds_ne0['mean'].rename('ne0'),
    ds_kr['mean'].rename('kr'),
    mws_max.rename('AS_max'),
    mws_std.rename('AS_std')
    ])

ds_params


#%%

# Perform averages over run
# TODO: combine with wma accessor

count = ds_params['nK_m3'].count('run')

ds_p_stats = xr.Dataset(coords=count.coords)

for var in ds_params.data_vars:
    mean = ds_params[var].mean('run')
    std = ds_params[var].std('run')
    stderr = std/np.sqrt(count)

    ds_p_stats["{}_mean".format(var)] = mean
    ds_p_stats["{}_std".format(var)] = std
    ds_p_stats["{}_stderr".format(var)] = stderr
#%%

ds_p_stats


#%%

plt.errorbar(
    ds_p_stats['nK_m3_mean'], 
    ds_p_stats['AS_max_mean'], 
    xerr=ds_p_stats['nK_m3_stderr'], 
    yerr=ds_p_stats['AS_max_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

plt.ylabel("AS Maximum")
plt.xlabel("nK_m3")

#%%

plt.errorbar(
    ds_p_stats['nK_m3_mean'], 
    ds_p_stats['ne0_mean'], 
    xerr=ds_p_stats['nK_m3_stderr'], 
    yerr=ds_p_stats['ne0_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

plt.ylim(3e12,2e13)
plt.xlim(4e20,2e22)

plt.ylabel("Ne0")
plt.xlabel("nK_m3")

#%%


plt.errorbar(
    ds_p_stats['nK_m3_mean'], 
    ds_p_stats['kr_mean'], 
    xerr=ds_p_stats['nK_m3_stderr'], 
    yerr=ds_p_stats['kr_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

# plt.ylim(3e12,2e13)
plt.xlim(4e20,2e22)

plt.ylabel("Kr")
plt.xlabel("nK_m3")

# %%
