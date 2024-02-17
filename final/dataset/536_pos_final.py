#%%


from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]


# add 536 photodiode data from 5x6_pos run. Do not have photodiode data for 536_pos run

tc = '5x6_pos'

ds_lecroy_5x6 = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy_5x6 = ds_lecroy_5x6.xr_utils.stack_run()

ds_lecroy_5x6 = ds_lecroy_5x6.sel(phi=0.8, method='nearest')

ds_lecroy_5x6 = ds_lecroy_5x6[['pd1','pd2']]

ds_lecroy = xr.merge([ds_lecroy_5x6, ds_lecroy])

ds_lecroy.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, '536_pos_lecroy.cdf'))

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1, pipe_fit_alpha_2

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds_absem

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

ds = ds_alpha_fit[['alpha_red', 'alpha_fit']]
ds = ds.rename({'alpha_red': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%


mws_max = ds_lecroy['AS'].mean('mnum').max('time').rename('AS_max')
nK = ds_p['nK_m3'].sel(mp='mw_horns')

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
    nK,
    delta_pd1.rename('delta_pd1') 
    ]).sortby('motor').dropna('run', how='all')

ds['AS_max'].attrs = dict(long_name='AS Max')
ds['mag_pp'].attrs = dict(long_name='Mag. Pre Pulse (PP)', units='V')
ds['mag_pp_std'].attrs = dict(long_name='Pre Pulse (PP) Std Dev.', units='V')
ds['AS_max_std_ratio'].attrs = dict(long_name='AS Max / PP Std Dev.', units='1/V')

ds

ds.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, '536_pos_ds_p_stats.cdf'))