#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
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

#%%

from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','cuttimes_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

def downselect_acq_time(ds, df_cuttimes, tc):
    """
    Function to downselect ds to times in df_cuttimes
    this is used for ds that was generated with one broad timewindow, then downselect to more specific plateau times
    """
    acq_times = ds['acq_time'].unstack('run').stack(temp=[...]).dropna('temp').values

    acq_times_set = set(acq_times)
    acq_times_in_ct = []

    for _, ct in df_cuttimes.iterrows():
        acq_times_in_ct.extend(time for time in acq_times_set if ct['Start Time'] <= time <= ct['Stop Time'])

    acq_times_in_ct = list(acq_times_in_ct)

    ds = ds.where(ds['acq_time'].isin(acq_times_in_ct), drop=True)

    print("Downselecting to times in cuttimes")
    print("Length before: ", len(acq_times))
    print("Length after: ", len(acq_times_in_ct))

    return ds


ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs, tc)

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

ds_lecroy = ds_lecroy.drop(0,'kwt')

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs, tc)

ds_fit = ds_lecroy.mean('mnum')
da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

#%%


g = da_fit.plot(hue='run_plot', row='kwt', x='time')

dropna(g)

#%%[markdown]

# ## Fixed kr, vary ne0, old pipeline

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_1 


ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_1(da_fit)

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

da_fit = ds_lecroy['AS']

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)

#TODO: sterr is nan where ds_p is not?
ds_p['kr'] = ds_p['kr'].where(~ds_p_stderr['kr'].isnull())

#%%
ds_kr = ds_p['kr'].to_dataset(name='mean')
ds_kr['std'] = ds_p_stderr['kr']

plot.common.xr_errorbar(ds_kr['mean'], ds_kr['std'], huedim='run')

plt.yscale('log')
# plt.ylim(1e-15,2e-14)


#%%[markdown]

# # Compare Lecroy and MWS
#%%

cantera_data_dir = os.path.join(REPO_DIR, 'modeling','dataset','output')
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf')).sel({'phi': 0.7})
kr = ds_TP_params['kr']
kr = kr.pint.quantify('cm^3/s').pint.to('um^3/us')
kr_sel = kr.sel(P=1e5, method='nearest').sel(T=[1525, 1750, 1975])

kr_sel = kr_sel.assign_coords(Kwt = kr_sel.coords['Kwt']*1e2)

kr_sel



# %%

ds_lecroy['AS'].dropna('run', how='all').dropna('kwt', how='all')

#%%

da_stats = ds_lecroy['AS'].sel(time=slice(-1,1))

mws_max = da_stats.mean('mnum').max('time')

mws_std = da_stats.std('mnum').max('time')

#%%

ds_params = xr.merge([
    ds_nK['mean'].sel(mp='barrel').drop('mp').rename('nK_m3_barrel'),
    ds_nK['mean'].sel(mp='mw_horns').drop('mp').rename('nK_m3_mw_horns'),
    ds_ne0['mean'].rename('ne0'),
    ds_kr['mean'].rename('kr'),
    mws_max.rename('AS_max'),
    mws_std.rename('AS_std')
    ])

ds_params


#%%

# Perform averages over run
# TODO: combine with wma accessor

count = ds_params['nK_m3_barrel'].count('run')

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

vars = ['nK_m3_barrel', 'nK_m3_mw_horns', 'kr', 'AS_max']

fig, axes = plt.subplots(len(vars), 1, figsize=(5,10), sharex=True)

for i, var in enumerate(vars):

    ax = axes[i]

    ax.errorbar(
        ds_p_stats.coords['kwt'], 
        ds_p_stats['{}_mean'.format(var)], 
        yerr=ds_p_stats['{}_stderr'.format(var)], 
        marker='o', capsize=5,
        label=var
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(var)


axes[0].set_ylabel("nK_m3 Barrel [#/m^3]")
axes[1].set_ylabel("nK_m3 MW Horns [#/m^3]")
axes[2].set_ylabel("kr [um^3/us]")
axes[3].set_ylabel("AS Maximum")

kr_sel.sel(Kwt=slice(0.1,1)).plot(hue='T', marker='o', ax=axes[2])
axes[2].set_title('')

axes[-1].set_xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params.png'), dpi=300, bbox_inches='tight')

#%%

plt.errorbar(
    ds_p_stats['nK_m3_barrel_mean'], 
    ds_p_stats['AS_max_mean'], 
    xerr=ds_p_stats['nK_m3_barrel_stderr'], 
    yerr=ds_p_stats['AS_max_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

plt.ylabel("AS Maximum")
plt.xlabel("nK_m3")

#%%

plt.errorbar(
    ds_p_stats['nK_m3_barrel_mean'], 
    ds_p_stats['ne0_mean'], 
    xerr=ds_p_stats['nK_m3_barrel_stderr'], 
    yerr=ds_p_stats['ne0_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

# plt.ylim(3e12,2e13)
plt.xlim(4e20,2e22)

plt.ylabel("Ne0 [#/um**3]")
plt.xlabel("nK_m3_barrel Barrel [#/m**3]")

#%%


plt.errorbar(
    ds_p_stats['nK_m3_barrel_mean'], 
    ds_p_stats['kr_mean'], 
    xerr=ds_p_stats['nK_m3_barrel_stderr'], 
    yerr=ds_p_stats['kr_stderr'], 
    marker='o', capsize=5
    )

plt.xscale('log')
plt.yscale('log')

# plt.ylim(3e12,2e13)
# plt.xlim(4e20,2e22)

plt.ylabel("kr [cm**3/s]")
plt.xlabel("nK Barrel [#/m^3]")

# plt.axhline(kr_sel.item().magnitude, linestyle = '--')
plt.twiny()
kr_sel.plot(hue='T', marker='o')
plt.gca().get_legend().set_bbox_to_anchor((0.1, 0.6))
plt.xscale('log')


# %%


