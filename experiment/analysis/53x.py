#%%


from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.mws_utils import calc_mag_phase_AS

# %%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.stack(run=['date','run_num'])
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_absem['diff'] = ds_absem['led_on'] - ds_absem['led_off']
ds_absem['alpha'] = 1 - ds_absem['diff']/ds_absem['calib']

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num'])
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = calc_mag_phase_AS(ds_lecroy)['AS']


#%%

ds_absem

# %%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='kwt', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)
# %%
from mhdpy.plot import dropna

g = ds_lecroy.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

dropna(g)

#%%



ds_fit = ds_absem.mean('mnum')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

#%%

#Downselect to remove data where there are no peaks... TODO: revisit
# alpha_tc = alpha_tc.where(alpha_tc.sel(wavelength=slice(760,775)).max('wavelength') > 0.5).dropna('time', how='all')
# alpha_tc = alpha_tc.sel(time=alpha_tc.coords['time'].values[::100])


#%%

spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength', how='all')
alpha_tc_red.name = 'alpha_red'

wls = alpha_tc.coords['wavelength'].values
fits, ds_p, ds_p_stderr = analysis.xr.fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, alpha_tc_red, fits]).sel(wavelength=slice(760,775))


#%%

ds = ds_alpha[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds

ds_p.coords['kwt'].attrs = dict(long_name='K Mass Frac.', units='%')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

# %%
g = ds_alpha[['alpha','alpha_fit']].to_array('var').sel(mp='barrel').plot(col='kwt',hue='var', row='run', ylim = (1e-3,2), yscale='log', figsize=(10,6))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')

#%%

da = ds_p['nK_m3']

da = da.dropna('run', how='all')

g  = da.plot(hue='run_plot', x='kwt',col='mp', marker='o')

dropna(g)

plt.yscale('log')
plt.xscale('log')

# %%

mws_max = ds_lecroy.mean('mnum').max('time')
nK = ds_p['nK_m3'].sel(mp='barrel')

#%%

ds = xr.merge([mws_max,nK]).sortby('kwt')

ds.unstack('run').plot.scatter(x='nK_m3', y='AS', hue='date')

#%%

for run in ds.coords['run_plot']:

    ds_sel = ds.sel(run=run).dropna('kwt')
    if len(ds_sel.coords['kwt']):
        ds_sel.set_coords('nK_m3')['AS'].sortby('nK_m3').plot(x='nK_m3', label=run.item(), marker='o')
        # plt.plot(ds_sel['nK_m3'],ds_sel['AS'], label=run.item() ,marker='o')

plt.legend()

plt.xscale('log')

#%%

ds.sel(run=('2023-05-18', 1)).dropna('kwt')

#%%




from mhdpy.mws_utils.fitting import fit_fn 
from mhdpy.analysis.xr import fit_da_lmfit
from lmfit import Model

da_fit = ds_lecroy.mean('mnum').copy()


pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

tc_dim = [dim for dim in da_fit.dims if dim != 'time'][0]

da_fit = da_fit.where(pulse_max > 5e-4) # Targeting low power...
da_fit = da_fit.dropna(tc_dim,'all')
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
g= da.plot(col='run', x='time',hue='var', row='kwt')

# plt.xscale('log')
plt.yscale('log')

dropna(g)

#%%

ds_p['ne0'].plot(hue='run_plot', x='kwt', marker='o')
plt.yscale('log')
plt.xscale('log')