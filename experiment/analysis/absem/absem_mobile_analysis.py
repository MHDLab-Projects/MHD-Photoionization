#%%



from mhdlab.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

from mhdlab.analysis import mwt
from mhdlab.analysis import absem
from mhdlab.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_absem = ppu.fileio.load_absem(tc)


# %%

ds = ds_absem.sel(mp='mw_horns').mean('mnum').dropna('run', how='all')


ds
# %%

ds[['led_on', 'led_off']].to_array('var').plot(col='run', row='motor', hue='var', x='wavelength')

#%%

ds[['diff', 'calib']].to_array('var').plot(col='run', row='motor', hue='var', x='wavelength')

#%%

ds['alpha'].plot(col='run', row='motor', x='wavelength')

plt.ylim(-0.2,1.1)

#%%

wl_range = slice(750,755)

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
