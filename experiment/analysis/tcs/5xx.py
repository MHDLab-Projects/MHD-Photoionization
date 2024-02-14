#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem

from mhdpy.plot import dropna

#%%
fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

fp_dsst = pjoin(DIR_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()


# %%

tc = '5x6_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = ds_lecroy.mws.calc_mag_phase_AS()['AS']

da_lecroy = da_lecroy.drop('run') # Only one run for this tc

#%%

da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')

#%%


da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)

plt.xlim(-1,)

#%%

da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi')

plt.yscale('log')


plt.xlim(-1,)


#%%

da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']
# da = da.max('wavelength').count('mnum')
da = da.max('wavelength').mean('mnum')

g = da.plot(hue='phi', marker='o')

dropna(g)

#%%
da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']

da.sel(motor=100, method='nearest', phi=1).dropna('mnum', how='all')


#%%

## Examine time data...
# I believe this is redundant with other time analysis, probably can remove 

ds_absem_time = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_absem.cdf'))

ds_absem_time['diff'] = ds_absem_time['led_on'] - ds_absem_time['led_off']
ds_absem_time['alpha'] = 1 - ds_absem_time['diff']/ds_absem_time['calib']

ds_absem_time = ds_absem_time.set_index(acq=['time','mp']).unstack('acq')

ds_absem_time = ds_absem_time['alpha']

ds_absem_time

#%%
dsst = mhdpy.fileio.TFxr(pjoin(DIR_PROC_DATA, 'dsst.tdms')).as_dsst()

# tw = slice(Timestamp('2023-05-24 19:45:01.091800832'), Timestamp('2023-05-24 20:39:19.309871616'), None)
tw = slice(Timestamp('2023-05-24 20:12:07.301042944'), Timestamp('2023-05-24 20:12:46.490260736'), None)

da = ds_absem_time.sel(time=tw).dropna('time', how='all')
da.mean('wavelength').plot(hue='mp', marker='o')

plt.twinx()

dsst['hvof']['CC_equivalenceRatio'].sel(time=tw).plot()

#%%

# Check assigned coord values

da = ds_absem_time.sel(time=tw).dropna('time', how='all')
da.mean('wavelength').plot(hue='mp', marker='o')

plt.twinx()

dst_coords['phi'].sel(time=tw).plot(marker='o')

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


# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)
# %%
# %%

tc = '5x3_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

# ds_absem = ds_absem.drop(34.81, 'motor')

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = ds_lecroy.mws.calc_mag_phase_AS()['AS']


da_lecroy = da_lecroy.drop('run') # Only one run for this tc

#%%


da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)

plt.xlim(-1,)

#%%

da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi')

plt.yscale('log')


plt.xlim(-1,)

#%%[markdown]

# # quick partial fix calibraiton data for 5x3_pos. 


# the calibration data appears to have a significant offset on the high wavlength side.
# I looked at the alpha in 53x notebook and it doesn't appear to be as bad as here
# TODO: some method of correcting calibration for this case. beta_offfset_line? or multiply calibration by a linear function of wavelength?

#%%

fp_calib = pjoin(REPO_DIR, 'experiment', 'data','proc_data', 'ds_calib.cdf')

ds_calib = xr.load_dataset(fp_calib)

ds_calib

#%%

da_calib = ds_calib.sel(mp='mw_horns').sel(time=Timestamp('2023-05-24T18:12:49.047832386'))
# da_calib = ds_calib.sel(mp='mw_horns').sel(time=Timestamp('2023-05-24T22:42:25.960157725'))

da_calib

ds_absem['calib'] = da_calib['diff']

ds_absem = ds_absem.absem.calc_alpha()

#%%

ds_sel = ds_absem.mean('run').sel(mp='mw_horns').mean('mnum')

# ds_sel[['diff','calib']].to_array('var').plot(col='motor', row='phi', hue='var')

#%%

ds_sel[['diff', 'calib']].to_array('var').sel(motor=180, method='nearest').plot(hue='var', row='phi')


#%%

ratio = (ds_sel['diff']/ds_sel['calib'])

ratio.sel(motor=50, method='nearest').plot(row='phi')

plt.ylim(0.8,1.2)


#%%
ds_sel['alpha'].sel(motor=50, method='nearest').plot(hue='phi')

plt.ylim(0,1)


# %%

da_lecroy_mean = da_lecroy.mean('mnum')

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)



#%%

da_max = da_lecroy_mean.sel(time=slice(-1,1)).max('time')

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')


#%%


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit = ds_absem.isel(run=0) #only one run
ds_fit = ds_fit.sel(wavelength=slice(750,790))

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

#%%

ds_alpha_fit['alpha_red'].sel(mp='mw_horns').plot(hue='phi', row='motor')

# plt.ylim

# plt.ylim()




# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)

plt.yscale('log')
# %%

da_count = ds_absem['alpha'].isel(wavelength=0).count('mnum')

da_count.sel(mp='mw_horns').plot(hue='phi', marker='o')
# %%
