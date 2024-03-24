

#%%[markdown]

# # Investigating the shape of the mws data for individual mnum

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis import mws

# %%

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='relative')
da_sel = ds_lecroy['dAS_rel'].sel(kwt=1, method='nearest')
# %%

da_sel.sel(mnum=slice(0,10)).plot(hue='mnum', row='run')

plt.xlim(-1,10)
plt.ylim(1e-3,)
plt.yscale('log')


#%%

da_2 = da_sel/da_sel.mws._pulse_max()

da_2.sel(mnum=slice(0,10)).plot(hue='mnum', row='run')


plt.xlim(-1,10)
plt.ylim(1e-1,)
plt.yscale('log')

#%%

da_2.mean('mnum').plot(hue='run_plot', x='time')
plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')

#%%

da_3 = da_sel.where(da_sel.mws._pulse_max() > 5e-4)

da_3 = da_3/da_3.mws._pulse_max()

da_3.mean('mnum').plot(hue='run_plot', x='time')
plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')
# %%

da_sel

#%%

run_max = da_sel.sel(time=slice(-2,2)).groupby('run_plot').max('time').max('mnum')

run_max

#%%

da_hp = da_sel/run_max

da_hp.isel(mnum=slice(0,10)).plot(row='run', hue='mnum', x='time')

plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')
#%%

da_hp = da_hp.where(da_hp.mws._pulse_max() > 0.3)
# %%

da_hp.mean('mnum').plot(hue='run_plot', x='time')
# plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')