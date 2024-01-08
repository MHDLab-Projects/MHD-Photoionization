#%%


from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

#%%

# ds_absem['calib'].isel(run=[0,1]).mean('mnum').sel(mp='barrel').plot(hue='run_plot', x='wavelength', row='kwt')

# %%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='kwt', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds_absem.mean('mnum')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_1(ds_fit, spect_red_dict)

#%%

ds = ds_alpha_fit[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds

ds_p.coords['kwt'].attrs = dict(long_name='K Mass Frac.', units='%')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

# %%
da_plot = ds_alpha_fit[['alpha','alpha_fit']].to_array('var')

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
ds_nK2 = ds_nK.wma.initialize_stat_dataset('mean', 'run')
ds_nK2

#%%[markdown]

# # Lecroy

#%%

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

# %%

g = ds_lecroy.mean('mnum')['AS'].plot(hue='run_plot', row='kwt', x='time')

dropna(g)



#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('AS', 'mnum')
ds_l_unc_2 = ds_l_unc.wma.calc_weighted_mean('run')

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_1 

da_fit = ds_l_unc_2['mean'].copy().rename('AS')

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_1(da_fit)

da = ds_mws_fit.to_array('var')

# for run in ds.coords['run_plot']:
#     da_sel = da.sel(run=run)
g= da.plot( x='time',hue='var', row='kwt')

# plt.xscale('log')
plt.yscale('log')

dropna(g)

#%%

ds_ne0 = ds_p['ne0'].to_dataset(name='mean')
ds_ne0['std'] = ds_p_stderr['ne0']

plt.errorbar(ds_ne0['kwt'], ds_ne0['mean'], ds_ne0['std'], marker='o', capsize=5)

#%%[markdown]

# # Compare Lecroy and MWS


# %%

mws_max = ds_l_unc_2.max('time')

max_time = ds_l_unc_2['mean'].argmax('time')

#TODO: pretty sure isel is right here
mws_max['std'] = ds_l_unc_2['std'].isel(time=max_time)

ds_nK2_barrel = ds_nK2.sel(mp='barrel')

#%%


ds = xr.merge([
    mws_max.rename({'mean':'AS'})['AS'],
    ds_nK2_barrel['mean'].rename('nK_m3')
    ]).sortby('kwt')

#%%

plt.errorbar(ds_nK2_barrel['mean'], mws_max['mean'], xerr=ds_nK2_barrel['std'], yerr=mws_max['std'], marker='o', capsize=5)

plt.xscale('log')
plt.yscale('log')

plt.ylabel("AS Maximum")
plt.xlabel("nK_m3")

#%%

plt.errorbar(ds_nK2_barrel['mean'], ds_ne0['mean'], xerr=ds_nK2_barrel['std'], yerr=ds_ne0['std'], marker='o', capsize=5)

plt.xscale('log')


plt.ylabel("ne_0 fit")
plt.xlabel("nK_m3")