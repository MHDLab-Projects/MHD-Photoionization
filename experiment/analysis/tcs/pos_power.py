#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem

# %%

tc = '536_pos_power'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = ds_lecroy.mws.calc_mag_phase_AS()['AS']

# %%

da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')
#%%
from mhdpy.plot import dropna

g = da_max.plot(col='motor',hue ='run_plot', x='power', col_wrap=3, marker='o')


# %%

# counts = da_lecroy.isel(time=0).count('mnum')


# counts.plot(hue='run_plot', col='power', x='motor')

# #%%

# da_sel = da_lecroy.where(counts > 10)
# da_sel = da_sel.dropna('motor', how='all')


# da_max = da_sel.mean('mnum').sel(time=slice(-1,1)).max('time')


# g = da_max.plot(col='motor',hue ='run_plot', x='power', col_wrap=3, marker='o')
#%%

da_sel = da_max.isel(motor=[0,2,4])
g = da_sel.plot(col='motor',hue ='run_plot', x='power', marker='o')

#%%

da_sel = da_max.isel(motor=[0,2,4])
g = da_sel.plot(hue='motor',row ='run', x='power', marker='o')

plt.xscale('log')
plt.yscale('log')
