#%%
from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu


# %%

ds_notorch = ppu.fileio.load_lecroy('NoTorch_pos', avg_mnum=True, AS_calc='absolute')

ds_5x0 = ppu.fileio.load_lecroy('5x0_pos', avg_mnum=True, AS_calc='absolute')

ds_5x3 = ppu.fileio.load_lecroy('5x3_pos', avg_mnum=True, AS_calc='absolute')
ds_5x6 = ppu.fileio.load_lecroy('5x6_pos', avg_mnum=True, AS_calc='absolute')

ds_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=True, AS_calc='absolute')





#%%

fig,axes = plt.subplots(1,5, figsize=(15,5), sharey=True)

ds_notorch['dAS_abs'].mws._pulse_max().mean('run').plot(ax=axes[0])
axes[0].set_title('No Torch')

ds_5x0['dAS_abs'].mws._pulse_max().mean('run').plot(hue='phi', ax=axes[1])
axes[1].set_title('5x0')
ds_5x3['dAS_abs'].mws._pulse_max().mean('run').plot(hue='phi', ax=axes[2])

axes[2].set_title('5x3')

ds_5x6['dAS_abs'].mws._pulse_max().mean('run').plot(hue='phi', ax=axes[3])

axes[3].set_title('5x6')


ds_536['dAS_abs'].mws._pulse_max().mean('run').plot(hue='phi', ax=axes[4])
axes[4].set_title('536')

for ax in axes:
    ax.set_ylabel('$\Delta AS_{max}$')


#%%


    
ds_notorch['AS_abs'].mws._pulse_max().mean('run').plot(label='No Torch')
ds_5x0['AS_abs'].mws._pulse_max().mean('run').sel(phi=0.8, method='nearest').plot(label='5x0')
ds_5x3['AS_abs'].mws._pulse_max().mean('run').sel(phi=0.8, method='nearest').plot(label='5x3')

ds_536['AS_abs'].mws._pulse_max().mean('run').plot(label='536')


plt.legend()

plt.title('phi=0.8')

plt.ylabel('$\Delta AS_{max}$')

# plt.yscale('log')
# %%
ds_notorch = ppu.fileio.load_lecroy('NoTorch_pos', avg_mnum=False, AS_calc='absolute')

ds_5x0 = ppu.fileio.load_lecroy('5x0_pos', avg_mnum=False, AS_calc='absolute')

ds_5x3 = ppu.fileio.load_lecroy('5x3_pos', avg_mnum=False, AS_calc='absolute')
ds_5x6 = ppu.fileio.load_lecroy('5x6_pos', avg_mnum=False, AS_calc='absolute')

ds_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')

#%%

from mhdpy.plot.common import xr_errorbar_axes

fig,axes = plt.subplots(1, figsize=(15,5), sharey=True)

da_mean = ds_notorch['AS_abs'].mws._pulse_max().isel(run=0).mean('mnum')
da_std = ds_notorch['AS_abs'].mws._pulse_max().isel(run=0).std('mnum')

xr_errorbar_axes(da_mean, da_std, label='No Torch', axes=axes)

da_mean = ds_5x0['AS_abs'].mws._pulse_max().isel(run=0).sel(phi=0.8,method='nearest').mean('mnum')
da_std = ds_5x0['AS_abs'].mws._pulse_max().isel(run=0).sel(phi=0.8,method='nearest').std('mnum')

xr_errorbar_axes(da_mean, da_std, label='5x0', axes=axes)

da_mean = ds_5x3['AS_abs'].mws._pulse_max().isel(run=0).sel(phi=0.8,method='nearest').mean('mnum')
da_std = ds_5x3['AS_abs'].mws._pulse_max().isel(run=0).sel(phi=0.8,method='nearest').std('mnum')

xr_errorbar_axes(da_mean, da_std, label='5x3', axes=axes)

da_mean = ds_536['AS_abs'].mws._pulse_max().mean('mnum').mean('run')
da_std = ds_536['AS_abs'].mws._pulse_max().std('mnum').mean('run')

xr_errorbar_axes(da_mean, da_std, label='536', axes=axes)


plt.legend()

#%%

ds_536['mag_peak'] = ds_536['mag'].sel(time=slice(0,.1)).mean('time')

ds_536['delta_mag'] = ds_536['mag_pp'] - ds_536['mag_peak'] 

da_mean = ds_536['delta_mag'].mean('mnum').mean('run')
da_std = ds_536['delta_mag'].std('mnum').mean('run')

fig,axes = plt.subplots(1, figsize=(15,5), sharey=True)

xr_errorbar_axes(da_mean, da_std, label='536', axes=axes)

#%%

ds_536['delta_mag_time'] = ds_536['mag_pp'] - ds_536['mag']

ds_536['delta_mag_time'].mean('mnum').sel(motor=180, method='nearest').mean('run').plot()

plt.yscale('log')

#%%

sfr = ds_536['delta_mag']/ds_536['mag_fluct']

sfr.mean('mnum').mean('run').plot()

plt.ylabel('SFR (delta mag / mag_pp_std)')


#%%

da_plot = ds_536['T_abs'].mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,200,250], method='nearest')

da_plot.plot(hue='motor')

#%%
AS = 1 - ds_536['T_abs']

da_plot = AS.mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,180,250], method='nearest')

da_plot.plot(hue='motor')

plt.ylabel('AS')

#%%

g = da_plot.sel(motor=[100,180], method='nearest').plot(col='motor', sharey=False)

for ax in g.axes.flatten():
    # ax.set_yscale('log')
    ax.set_ylabel('AS')

#%%

T_pp = ds_536['T_abs'].sel(time=slice(-50,-1)).mean('time')
AS_pp = 1 - T_pp

AS_pp.mean('mnum').mean('run').plot(label='AS Pre Pulse')

T_laser = ds_536['T_abs'].sel(time=slice(0,.1)).mean('time')
AS_laser = 1 - T_laser

AS_laser.mean('mnum').mean('run').plot(label='AS Laser')

plt.ylabel('AS')

plt.legend()


#%%

import numpy as np

# Pre Pulse
AS_pp_mean = AS_pp.mean('mnum').mean('run')
AS_pp_std = AS_pp.mean('mnum').std('run')
plt.errorbar(np.arange(len(AS_pp_mean)), AS_pp_mean, yerr=AS_pp_std, label='AS Pre Pulse', capsize=5, marker='o')

# Laser
AS_laser_mean = AS_laser.mean('mnum').mean('run')
AS_laser_std = AS_laser.mean('mnum').std('run')
plt.errorbar(np.arange(len(AS_laser_mean)), AS_laser_mean, yerr=AS_laser_std, label='AS Laser', capsize=5, marker='o')

plt.ylabel('AS')
plt.legend()

#%%

# Pre Pulse
AS_pp_mean = AS_pp.mean('mnum').mean('run')
AS_pp_std = AS_pp.mean('mnum').std('run')

# Laser
AS_laser_mean = AS_laser.mean('mnum').mean('run')
AS_laser_std = AS_laser.mean('mnum').std('run')

fig, ax = plt.subplots()


xr_errorbar_axes(AS_pp_mean, AS_pp_std, label='AS Pre Pulse', axes=ax)
xr_errorbar_axes(AS_laser_mean, AS_laser_std, label='AS Laser', axes=ax)

plt.ylabel('AS')

#%%


# ds_536_200 = ds_536.where(ds_536['mnum'] < 200, drop=True)

ds_plot = ds_536

AS_pp_fluct = (1- ds_plot['T_abs']).sel(time=slice(-50,-1)).std('time')
AS_pp = (1- ds_plot['T_abs']).sel(time=slice(-50,-1)).mean('time')
AS_laser = (1- ds_plot['T_abs']).sel(time=slice(0,.1)).mean('time')

AS_pp_mean = AS_pp.mean('mnum').sel(run=('2023-05-18', 2)).dropna('motor', how='all')
AS_pp_std = AS_pp_fluct.mean('mnum').sel(run=('2023-05-18', 2)).dropna('motor', how='all')

AS_pp_std = AS_pp_fluct.mean('mnum').sel(run=('2023-05-18', 2)).dropna('motor', how='all')
AS_AS_pp_std = AS_pp_fluct.mean('mnum').sel(run=('2023-05-18', 2)).dropna('motor', how='all')

fig, ax = plt.subplots()

xr_errorbar_axes(AS_pp_mean, AS_pp_std, label='AS Pre Pulse', axes=ax, capsize=5)
xr_errorbar_axes(AS_laser_mean, AS_laser_std, label='AS Laser', axes=ax, capsize=5)


plt.ylabel('AS')


#%%

count = ds_536['T_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all').count('mnum').isel(time=0)

count.plot()

plt.ylabel("Count")

#%%

# T_laser = ds_536['T_abs'].sel(time=slice(0,.1)).mean('time')
# AS_laser = 1 - T_laser

AS = 1 - ds_536['T_abs']

g = AS.mean('mnum').mean('run').plot(row='motor', sharey=False)


for ax in g.axes.flatten():
    ax.set_yscale('log')
    ax.set_ylabel('AS')
# plt.yscale('log')

plt.xlim(-10,)

#%%

deltaAS = AS - AS.sel(time=slice(-50,-1)).mean('time')

g = deltaAS.mean('mnum').mean('run').plot(row='motor', sharey=False)


for ax in g.axes.flatten():
    ax.set_yscale('log')
    ax.set_ylabel('AS')
# plt.yscale('log')

plt.xlim(-10,)
# %%


dAS_max = ds_536['AS_abs'].mws._pulse_max().sel(run=[('2023-05-18', 1), ('2023-05-18', 2)]).mean('mnum')

mag_pp_std = ds_536['mag_fluct'].sel(run=[('2023-05-18', 1), ('2023-05-18', 2)]).mean('mnum')

rat = dAS_max/mag_pp_std

g= rat.plot(hue='run_plot', x='motor')

dropna(g)

plt.axvline(180, color='k', linestyle='--')

#%%

window_mins = [5,10,20,30,40,50]

ds_plot = ds_536['mag'].sel(run=('2023-05-18', 1)).dropna('motor', how='all')

for window_min in window_mins:
    mag_pp_std = ds_plot.sel(time=slice(-1,window_min)).std('time')
    mag_pp_std = mag_pp_std.mean('mnum')

    mag_pp_std.plot(label='{} us'.format(window_min))

plt.legend()
