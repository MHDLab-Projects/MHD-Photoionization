# %%
from mhdpy.analysis.standard_import import *

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

munged_dir = pjoin(REPO_DIR, 'experiment','data','munged')

#%%

"""
I am adding the coordinates to the data here because pd1 and pd2 are removed in normal post processing.
Keeping them in caused the lecry dataset to double in size and caused some sort of error when running proc_add_coord.py (terminal just said 'killed', probably ran out of memory)
This was pulled from add_coord_overview.py
"""

from mhdpy.fileio.tdms import tdms2ds
from mhdpy.coords import assign_coords_multi
from mhdpy.fileio.ct import load_df_cuttimes
from datetime import datetime

DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data')
coords_to_assign = tdms2ds(pjoin(DIR_PROC_DATA, 'dst_coords.tdms'))

coord_signal_dict = {k: coords_to_assign[k].dropna('time',how='all') for k in coords_to_assign.data_vars}

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'ct_sequence.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes)
df_cuttimes['date'] = df_cuttimes['Start Time'].dt.date
df_cuttimes = df_cuttimes.where(df_cuttimes['date'].astype(str) == '2023-05-24').dropna()
df_cuttimes = df_cuttimes.set_index('Event')

df_cuttimes

#%%
from mhdpy.analysis import mws

fp = pjoin(munged_dir, datestr, 'ds_lecroy_time.cdf')
ds_lecroy = xr.load_dataset(fp)

ds_lecroy = ds_lecroy.mws.calc_AS_rel()

ds_lecroy = ds_lecroy[['AS_rel','pd1','pd2']]

ds_lecroy = ds_lecroy.rename({'AS_rel':'mw'})

ds_lecroy = ds_lecroy.pint.dequantify()



#%%[markdown]

# # 53x

#%%


s = df_cuttimes.loc['53x_1']
tw = slice(s['Start Time'], s['Stop Time'])

ds_sel = ds_lecroy.sel(acq_time=tw)
ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)

ds = ds.sel(phi=0.8, method='nearest')

ds = ds.dropna('mnum', how='all')

#%%

ds.mean('mnum').to_array('var').plot(row='kwt', hue='var')

#%%


ds[['pd1','pd2']].mean('mnum').to_array('var').plot(row='kwt', hue='var')


#%%

ds_2 = ds.mws._pulse_max().mean('mnum')

ds_2['mw'].plot(marker='o')
plt.twinx()
ds_2['pd1'].plot(marker='o', label ='photodiode', color='r')

plt.legend()

#%%[markdown]

# # 536 Power

#%%


s = df_cuttimes.loc['536_power_1']
tw = slice(s['Start Time'], s['Stop Time'])

ds_sel = ds_lecroy.sel(acq_time=tw)
ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)

#%%

ds_2 = ds.mws._pulse_max().mean('mnum')

ds_2['mw'].plot(marker='o')
plt.twinx()
ds_2['pd1'].plot(marker='o', label ='photodiode', color='r')

plt.legend()

plt.xscale('log')
plt.yscale('log')

# %%[markdown]

# # 5x6 pos

#%%

s = df_cuttimes.loc['5x6_pos_1']
tw = slice(s['Start Time'], s['Stop Time'])

ds_sel = ds_lecroy.sel(acq_time=tw)
ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)

#%%

ds[['pd1','pd2']].mean('mnum').to_array('var').plot(hue='var', col='motor',row='phi', marker='o')

# plt.ylim(0,)
#%%

ds['pd1'].mean('mnum').plot(hue='var', col='motor',row='phi', marker='o', sharey=False)

plt.xlim(-2,10)

#%%

ds_2 = ds.mws._pulse_max()

ds_2 = ds_2[['pd1','pd2','mw']].to_array('var').mean('mnum')

ds_2.plot(col='var', hue='phi', marker='o', sharey=False)

#%%


ds['dpd1'] = ds['pd1'] - ds['pd1'].sel(time=slice(-1,0)).mean('time')

# For some reason can't use _pulse_max here
ds_2 = ds.sel(time=slice(-1,1)).max('time')

ds_2 = ds_2[['dpd1','mw']].to_array('var').mean('mnum')

ds_2.plot(col='var', hue='phi', marker='o', sharey=False)

ds_2
#%%



# %%[markdown]

# # 5x3 pos

#%%

s = df_cuttimes.loc['5x3_pos_1']
tw = slice(s['Start Time'], s['Stop Time'])

ds_sel = ds_lecroy.sel(acq_time=tw)
ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)

#%%

ds_2 = ds.mws._pulse_max()

ds_2 = ds_2[['pd','mw']].to_array('var').mean('mnum')

ds_2.plot(col='var', hue='phi', marker='o', sharey=False)

#%%[markdown]

# # 536 pos

#%%
s = df_cuttimes.loc['516_pos_1']
tw = slice(s['Start Time'], s['Stop Time'])

ds_sel = ds_lecroy.sel(acq_time=tw)
ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)
ds = ds.dropna('mnum', how='all')

ds



# %% [markdown]
# # Compare time series of 770 nm photodiode (pd) and microwave scattering (mw) 

# %%

ds['pd'].isel(mnum=[0,10,20,30]).plot(row='motor', hue='mnum')
# plt.yscale('log')
# plt.xlim(-2,2)

# %%


ds['mw'].isel(mnum=[0,10,20,30]).plot(row='motor', hue='mnum')
# plt.yscale('log')
# plt.xlim(-2,2)

#%%

#TODO: Happening twice for AS...move up to apply to pd?
pre_pulse_avg = ds.sel(time=slice(-2,-1)).mean('time')

ds_max = ds - pre_pulse_avg

ds_max = ds_max.sel(time=slice(-1,1)).max('time')

ds_max = ds_max/ds_max.max()

ds_max.mean('mnum').to_array('var').plot(hue='var', marker='o')

# %%

ds['pd'].sel(motor=25.12, method='nearest').plot()

# %%

# g = ds.isel(motor=0).sel(mnum=slice(0,10)).to_array('var').plot(row='mnum', hue='var')

ds_plot = ds.isel(motor=1)

n_plots = 10

fig, axes = plt.subplots(n_plots, figsize = (5, 3*n_plots))

for i in range(n_plots):
    ax = axes[i]
    ds_plot['mw'].sel(mnum=i).plot(ax=ax)

    tax = ax.twinx()

    ds_plot['pd'].sel(mnum=i).plot(ax=tax, color='r')

# %% [markdown]
# # Frequency Spectrum
# 
# Look before laser pulse

# %%

import xrft

ds_pre = ds.sel(time=slice(-50,-1))

ds_pre['pd'].isel(mnum=[0,10,20,30]).plot(row='motor', hue='mnum')

plt.yscale('log')


# %%

da_ft_pd = xrft.power_spectrum(ds_pre['pd'], dim='time')
da_ft_pd.name = 'pd'
da_ft_mw = xrft.power_spectrum(ds_pre['mw'], dim='time')
da_ft_mw.name = 'mw'

da_ft = xr.merge([da_ft_mw,da_ft_pd])

da_ft.sel(motor=25.12, method='nearest').mean('mnum').to_array('var').plot(hue='var')

plt.yscale('log')
plt.xlim(0,)
plt.xlabel("Frequency (MHz)")


# %%

da_ft['pd'].mean('mnum').plot(hue='motor')

plt.yscale('log')
plt.xlim(0,)

plt.xlabel("Frequency (MHz)")

# %% [markdown]
# # Convolution

# %%

ds_sel = ds.sel(time=slice(-50,-1)).sel(motor=25.12, method='nearest').sel(mnum=0)

a1 = ds_sel['pd']
a2 = ds_sel['mw']

out = np.convolve(a1, a2)

plt.plot(out)
# %%

s1 = (a1 - np.mean(a1))
s1 = s1/s1.max()

s2 = (a2 - np.mean(a2))
s2 = s2/s2.max()

plt.plot(s1)
plt.plot(s2)


# %%

out = np.convolve(s1, s2)

plt.plot(out)
# %%
