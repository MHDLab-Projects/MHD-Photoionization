# %%
from mhdpy.analysis.standard_import import *

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

munged_dir = pjoin(REPO_DIR, 'experiment','data','munged')
fp = pjoin(munged_dir, datestr, 'ds_lecroy_time.cdf')

#%%

"""
I am adding the coordinates to the data here because pd1 and pd2 are removed in normal post processing.
Keeping them in caused the lecry dataset to double in size and caused some sort of error when running proc_add_coord.py (terminal just said 'killed', probably ran out of memory)
This was pulled from add_coord_overview.py
"""

from mhdpy.fileio.tdms import tdms2ds
from mhdpy.coords import assign_coords_multi
from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list

DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data')
coords_to_assign = tdms2ds(pjoin(DIR_PROC_DATA, 'dst_coords.tdms'))

coord_signal_dict = {k: coords_to_assign[k].dropna('time',how='all') for k in coords_to_assign.data_vars}

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'cuttimes.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes)
cuttimes = extract_cuttime_list(df_cuttimes)

tw = cuttimes[16] # 516_pos

#%%

ds_in = xr.load_dataset(fp)
ds_sel = ds_in.sel(acq_time=tw)

ds = assign_coords_multi(ds_sel, coord_signal_dict, min_mnum=2)

ds = ds.dropna('mnum', how='all')

ds

# %%

from mhdpy.analysis.mws import calc_mag_phase_AS
ds_mw = calc_mag_phase_AS(ds)

dapd = ds['pd1']
dapd.name='pd'

damw = ds_mw['AS']
damw.name = 'mw'

dapd.coords['time'] = damw.coords['time']

ds = xr.merge([dapd, damw])

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
