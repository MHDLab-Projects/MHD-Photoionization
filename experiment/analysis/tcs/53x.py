#%%


from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.mws import calc_mag_phase_AS
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.stack(run=['date','run_num']).dropna('run',how='all')
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_absem = ds_absem.absem.calc_alpha()

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

ds_p, ds_p_stderr, ds_alpha = pipe_fit_alpha_1(ds_fit, spect_red_dict)

#%%

ds = ds_alpha[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds

ds_p.coords['kwt'].attrs = dict(long_name='K Mass Frac.', units='%')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

# %%
da_plot = ds_alpha[['alpha','alpha_fit']].to_array('var')

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

ds_nK


#%%

g = xr.plot.FacetGrid(ds_nK, col='mp', row='run')

g.map(plt.errorbar, 'kwt', 'mean', 'std', capsize=5)

plt.yscale('log')
plt.xscale('log')

dropna(g)

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

#%%[markdown]

# # Lecroy

#%%

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num']).dropna('run',how='all')
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = calc_mag_phase_AS(ds_lecroy)

# %%

g = ds_lecroy.mean('mnum')['AS'].plot(hue='run_plot', row='kwt', x='time')

dropna(g)

#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('mag', 'mnum')
ds_l_unc['count'] = ds_l_unc['count'].isel(time=0)
ds_l_unc

#%%

# make facetgrid plot with confidence intervals

# Define your custom plotting function
def custom_plot(time, mean, std, count, **kwargs):
    lower_bound = mean - std
    upper_bound = mean + std
    plt.plot(time, mean, **kwargs)
    plt.fill_between(time, lower_bound, upper_bound, alpha=0.5, color='gray')  # alpha sets the transparency
    # count should be a single value, add it as text to the plot in the upper left
    if not np.isnan(count):
        plt.text(0.05, 0.9, '{}'.format(int(count)), transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')



g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

plt.xlim(-10,10)

#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('AS', 'mnum')
ds_l_unc['count'] = ds_l_unc['count'].isel(time=0)
ds_l_unc

#%%

g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

plt.xlim(-10,10)

#%%

g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'stderr', 'count')


#%%

ds_l_unc_2 = ds_l_unc.wma.calc_weighted_mean('run')

#%%

g = xr.plot.FacetGrid(ds_l_unc_2, row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

# plt.ylim(0,1)


#%%


from mhdpy.analysis.mws.fitting import fit_fn 
from mhdpy.xr_utils import fit_da_lmfit
from lmfit import Model

da_fit = ds_l_unc_2['mean'].copy()


pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

tc_dim = [dim for dim in da_fit.dims if dim != 'time'][0]

da_fit = da_fit.where(pulse_max > 5e-4) # Targeting low power...
da_fit = da_fit.dropna(tc_dim, how='all')
pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

da_fit = da_fit/pulse_max

da_fit = np.log(da_fit).dropna('time', 'all')

mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0, take_log=True))

params = mod.make_params()
params['dne'].value = 1e14
params['ne0'].value = 3e12
params['dne'].min = 0
params['ne0'].min = 0


da_fit_region = da_fit.sel(time=slice(0,30))
x_eval = da_fit.sel(time=slice(0,30)).coords['time'].values

# da_fit_coarse = da_fit_region.coarsen(time=10).mean()
fits, ds_p, ds_p_stderr = fit_da_lmfit(da_fit_region, mod, params, 'time', x_eval)

fits = np.exp(fits)
da_fit = np.exp(da_fit)
# %%


ds = xr.merge([fits, da_fit])

da = ds.to_array('var')

# for run in ds.coords['run_plot']:
#     da_sel = da.sel(run=run)
g= da.plot( x='time',hue='var', row='kwt')

# plt.xscale('log')
plt.yscale('log')

dropna(g)

#%%

ds_p['ne0'].plot( x='kwt', marker='o')

#%%

ds_p['ne0'].drop(0,'kwt').plot(marker='o')

plt.yscale('log')
plt.xscale('log')

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

ds.plot.scatter(x='nK_m3', y='AS')


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