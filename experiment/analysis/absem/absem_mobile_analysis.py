#%%



from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]

# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)


# %%

ds = ds_absem.sel(mp='mw_horns').mean('mnum').dropna('run', how='all')

ds = ds.drop(34.81, 'motor')

ds
# %%

ds[['led_on', 'led_off']].to_array('var').plot(col='run', row='motor', hue='var', x='wavelength')

#%%

ds[['diff', 'calib']].to_array('var').plot(col='run', row='motor', hue='var', x='wavelength')

#%%

ds['alpha'].plot(col='run', row='motor', x='wavelength')

plt.ylim(-0.2,1.1)

#%%

wl_range = slice(740,750)

rat = ds['diff']/ds['calib']

rat = rat.sel(wavelength=wl_range).mean('wavelength')

rat
# %%

ds_fix = ds.copy()

ds_fix['calib'] = ds_fix['calib'] * rat

ds_fix = ds_fix.absem.calc_alpha()

#%%

ds_fix

#%%

ds_fix[['diff', 'calib']].to_array('var').plot(col='run', row='motor', hue='var', x='wavelength')

#%%

ds_fix['alpha'].plot(col='run', row='motor', x='wavelength')

plt.ylim(-0.2,1.1)

# %%
