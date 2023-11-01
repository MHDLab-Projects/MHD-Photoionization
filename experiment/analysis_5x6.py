#%%

from mhdpy.analysis.standard_import import *

from mhdpy.mws_utils import calc_mag_phase_AS

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
# %%


da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')

# %%
