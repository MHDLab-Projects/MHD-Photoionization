#%%


from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '516_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]

# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds_absem.mean('mnum')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_1(ds_fit, spect_red_dict)

#%%

mws_max = ds_lecroy['AS'].mean('mnum').max('time').rename('AS_max')
nK_mw_horns = ds_p['nK_m3'].sel(mp='mw_horns').rename('nK_mw_horns').drop('mp')
nK_barrel = ds_p['nK_m3'].sel(mp='barrel').rename('nK_barrel').drop('mp')

mws_pp = ds_lecroy['mag_pp'].mean('mnum').rename('mag_pp')
mws_pp_std = ds_lecroy['mag'].sel(time=slice(-50,-1)).std('time').mean('mnum').rename('mag_pp_std')

mws_max_std_ratio = mws_max/mws_pp_std
mws_max_std_ratio = mws_max_std_ratio.rename('AS_max_std_ratio')


delta_pd1 = ds_lecroy['pd1'] - ds_lecroy['pd1'].sel(time=slice(-1,0)).mean('time')
delta_pd1 = delta_pd1.dropna('run', how='all')
delta_pd1 = delta_pd1.mean('mnum').max('time')
# delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')

#%%

ds_lecroy['pd1'].dropna('run', how='all')

#%%

ds = xr.merge([
    mws_max, 
    mws_pp,
    mws_pp_std,
    mws_max_std_ratio,
    nK_mw_horns,
    nK_barrel,
    delta_pd1.rename('delta_pd1') 
    ]).sortby('motor').dropna('run', how='all')

ds['AS_max'].attrs = dict(long_name='AS Max')
ds['mag_pp'].attrs = dict(long_name='Mag. Pre Pulse (PP)', units='V')
ds['mag_pp_std'].attrs = dict(long_name='Pre Pulse (PP) Std Dev.', units='V')
ds['AS_max_std_ratio'].attrs = dict(long_name='AS Max / PP Std Dev.', units='1/V')

ds

ds.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, '516_pos_ds_p_stats.cdf'))
# %%
