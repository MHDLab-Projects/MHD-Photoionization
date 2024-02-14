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
from mhdpy.coords.ct import downselect_acq_time

#%%[markdown]

# # Absem

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem = ds_absem.drop(0, 'kwt')

#%%

from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs)

# %%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='kwt', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)

#%%

# # Examine calibration offset. This does not appear to be off as much as in the 5xx notebook. TODO: investigate and quantify off-peak calibraiton offset. 
# ds_absem['alpha'].mean('mnum').sel(mp='mw_horns').isel(run=-1).plot(hue='kwt')


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
da.coords['kwt'].attrs = dict(long_name="Kwt", units='%')

# da = da.dropna('run', how='all')

g  = da.plot(hue='run_plot', x='kwt',col='mp', marker='o')

dropna(g)

plt.yscale('log')
plt.xscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, '53x_AES_nK.png'), dpi=300, bbox_inches='tight')

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

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

ds_fit = ds_lecroy.mean('mnum')
da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

#%%

g = da_fit.plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')

dropna(g)

#%%[markdown]

# # Exponential Fit

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))

#%%

ds_mws_fit[['AS_all','AS_fit']].to_array('var').plot(hue='var', row='kwt', col='run')

plt.yscale('log')

#%%


#%%

# Load cfd K+

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)

ds_cfd = ds_cfd.interp(kwt=ds_lecroy.coords['kwt']).dropna('kwt', how='all')

ds_cfd = ds_cfd.cfd.convert_species_rho()


goldi_pos = ds_cfd['x'].min().item() + 0.18
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')


ds_cfd['Yeq_K+'].plot()


#%%

# Cantera data

cantera_data_dir = os.path.join(REPO_DIR, 'modeling','dataset','output')
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf')).sel({'phi': 0.7})
kr = ds_TP_params['kr']
kr = kr.pint.quantify('cm^3/s').pint.to('um^3/us')
kr_sel = kr.sel(P=1e5, method='nearest').sel(T=[1525, 1750, 1975])

kr_sel = kr_sel.assign_coords(Kwt = kr_sel.coords['Kwt']*1e2)

kr_sel = kr_sel.interp(Kwt=ds_lecroy.coords['kwt'], kwargs={'fill_value': 'extrapolate'})
kr_sel = kr_sel.pint.quantify('um^3/us')

kr_sel

#%%

from mhdpy.plot.common import xr_errorbar_axes
from mhdpy.plot.common import xr_errorbar

fig, axes = plt.subplots(3, 1, figsize=(5,10), sharex=True)


da_mean = ds_p.mean('run')['decay']
da_std = ds_p.std('run')['decay']

xr_errorbar_axes(da_mean, da_std, axes=axes[0])

plt.xscale('log')

axes[0].set_ylabel('Decay Constant [us]')

kr_sel.pint.to('cm^3/us').plot(label='Cantera', hue='T', marker='o', ax=axes[1])
axes[1].set_title('')


Kp_decay = 1/(2*ds_p['decay'].pint.quantify('us')*kr_sel)
Kp_decay = Kp_decay.pint.to('1/cm^3')

da_mean = Kp_decay.mean('run')
da_std = Kp_decay.std('run')

xr_errorbar_axes(da_mean, da_std, huedim='T', axes=axes[2])
axes[2].get_legend().remove()

ds_cfd['Yeq_K+'].plot(label='CFD Yeq_K+', marker='o', ax=axes[2])
# ds_cfd['Yeq_OH'].plot(label='CFD Yeq_K+', marker='o', ax=axes[2])

axes[2].legend()

plt.yscale('log')
plt.title('')


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

ds_cfd['Yeq_K+'].pint.quantify('1/cm^3').pint.to('1/um^3').plot(label='CFD', marker='o')

plt.yscale('log')

#%%[markdown]

# ## Fixed ne0, vary kr, new pipeline

#%%

da_fit = ds_lecroy['AS']

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

dss_mws_fit = []
dss_p = []
dss_p_stderr = []

for kwt in ds_cfd.coords['kwt']:
    da_fit_sel = da_fit.sel(kwt=kwt)

    Kp_val = ds_cfd['Yeq_K+'].sel(kwt=kwt).item()
    ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit_sel, take_log=False, ne0=Kp_val)

    #TODO: sterr is nan where ds_p is not?
    ds_p['kr'] = ds_p['kr'].where(~ds_p_stderr['kr'].isnull())

    dss_mws_fit.append(ds_mws_fit)
    dss_p.append(ds_p)
    dss_p_stderr.append(ds_p_stderr)

ds_mws_fit = xr.concat(dss_mws_fit, dim='kwt')
ds_p = xr.concat(dss_p, dim='kwt')
ds_p_stderr = xr.concat(dss_p_stderr, dim='kwt')

#%%
ds_kr = ds_p['kr'].to_dataset(name='mean')
ds_kr['std'] = ds_p_stderr['kr']

plot.common.xr_errorbar(ds_kr['mean'], ds_kr['std'], huedim='run')

plt.yscale('log')
# plt.ylim(1e-15,2e-14)


#%%[markdown]

# # Compare Lecroy and MWS
#%%





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

fig, axes = plt.subplots(3, 1, figsize=(5,10), sharex=True)

var = 'nK_m3_barrel'
axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='Barrel'
    )

var = 'nK_m3_mw_horns'
axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='MW Horns'
    )

var = 'kr'
axes[1].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5
    )

kr_sel.sel(kwt=slice(0.1,1)).plot(hue='T', marker='o', ax=axes[1])
axes[1].set_title('')

var = 'AS_max'
axes[2].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5
    )

for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')


axes[0].set_ylabel("nK_m3 [#/m^3]")
axes[1].set_ylabel("kr [um^3/us]")
axes[2].set_ylabel("AS Maximum")


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


