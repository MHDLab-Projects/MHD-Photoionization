# %%

from mhdpy.analysis.standard_import import *

data_folder = mhdpy.fileio.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-05-24')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%

tc = '516_pos_1'

from mhdpy.mws_utils import calc_mag_phase_AS


fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)


# %%
ds_in

# %%

ds_mw = calc_mag_phase_AS(ds_in)

dapd = ds_in['pd1']
dapd.name='pd'

damw = ds_mw['AS']
damw.name = 'mw'

dapd.coords['time'] = damw.coords['time']

# dapd = dapd/dapd.max('time')
# damw = damw/damw.max('time')

ds = xr.merge([dapd, damw])


# %%

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


ds['pd'].sel(motor=25.12).plot()


# %%
from mhdpy.plot import dropna

# g = ds.isel(motor=0).sel(mnum=slice(0,10)).to_array('var').plot(row='mnum', hue='var')

ds_plot = ds.isel(motor=1)

n_plots = 10

fig, axes = plt.subplots(n_plots, figsize = (5, 3*n_plots))

for i in range(n_plots):
    ax = axes[i]
    ds_plot['mw'].sel(mnum=i).plot(ax=ax)

    tax = ax.twinx()

    ds_plot['pd'].sel(mnum=i).plot(ax=tax, color='r')
#

# dropna(g)

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

da_ft.sel(motor=25.12).mean('mnum').to_array('var').plot(hue='var')

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

ds_sel = ds.sel(time=slice(-50,-1)).sel(motor=25.12).sel(mnum=0)

a1 = ds_sel['pd']
a2 = ds_sel['mw']

import numpy as np

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
