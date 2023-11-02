#%%

from mhdpy.analysis.standard_import import *

from mhdpy.mws_utils import calc_mag_phase_AS

from mhdpy.plot import dropna
# %%

tc = '5x6_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.stack(run=['date','run_num'])
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num'])
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = calc_mag_phase_AS(ds_lecroy)['AS']


da_lecroy = da_lecroy.drop('run') # Only one run for this tc

#%%

ds_absem
# %%


da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')


#%%

da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']

# da = da.max('wavelength').count('mnum')
da = da.max('wavelength').mean('mnum')



da
g = da.plot(hue='phi', marker='o')


dropna(g)

# plt.ylabel("measurement count")
#%%
da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']

da.sel(motor=100, method='nearest', phi=1).dropna('mnum','all')


#%%

## Examine time data...

ds_absem_time = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_absem.cdf'))['alpha']
dsst = mhdpy.io.TFxr(pjoin(DIR_PROC_DATA, 'dsst.tdms')).as_dsst()

# tw = slice(Timestamp('2023-05-24 19:45:01.091800832'), Timestamp('2023-05-24 20:39:19.309871616'), None)
tw = slice(Timestamp('2023-05-24 20:12:07.301042944'), Timestamp('2023-05-24 20:12:46.490260736'), None)


da = ds_absem_time.sel(time=tw).dropna('time','all')

da.mean('wavelength').plot(hue='mp', marker='o')

# ds_absem['alpha'].mean('wavelength')

plt.twinx()

# dsst['motor']['Motor C Relative'].sel(time=tw).plot()

dsst['hvof']['CC_equivalenceRatio'].sel(time=tw).plot()


#%%
#%%

ds_fit = ds_absem.mean('mnum')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)


spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength','all')
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


# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)
# %%

#%%

# %%

tc = '5x3_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.stack(run=['date','run_num'])
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num'])
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = calc_mag_phase_AS(ds_lecroy)['AS']


da_lecroy = da_lecroy.drop('run') # Only one run for this tc


# %%

da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)



#%%

da_max = da_lecroy_mean.sel(time=slice(-1,1)).max('time')

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')


#%%

ds_fit = ds_absem.mean('mnum')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)


spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength','all')
alpha_tc_red.name = 'alpha_red'

wls = alpha_tc.coords['wavelength'].values
fits, ds_p, ds_p_stderr = analysis.xr.fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, alpha_tc_red, fits]).sel(wavelength=slice(760,775))



#%%


# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)