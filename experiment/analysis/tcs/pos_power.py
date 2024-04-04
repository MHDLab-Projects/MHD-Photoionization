#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.analysis import absem

# %%

tc = '536_pos_power'

ds_absem = ppu.fileio.load_absem(tc)

ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')
da_lecroy = ds_lecroy['AS_abs'] 

# %%

da_max = da_lecroy.sel(time=slice(-1,1)).max('time')
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

# %%
