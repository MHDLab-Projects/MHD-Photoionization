
#%%
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import pickle
import os
import scipp as sc
from pandas import Timestamp

from mhdpy.io import TFxr, gen_path_date
data_folder = gen_path_date('2023-04-07')
dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

ds_calib_mp = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_calib.cdf'))
ds_mp_mean = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_mean.cdf'))
ds_mp_std = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_std.cdf'))

#%%

# tw = slice(Timestamp('2023-04-07 20:28:32.165745664'), Timestamp('2023-04-07 20:33:42.461606400'), None)
tw = slice(Timestamp('2023-04-07 20:09:17.302072576'), Timestamp('2023-04-07 20:11:26.381779968'), None)


wavelength_sel = slice(80, 310)

ds_mean_region = ds_mp_mean.sel(time=tw)
ds_std_region = ds_mp_std.sel(time=tw)

#%%

ds_mean_region

#%%

ds_sc = sc.compat.from_xarray(ds_mean_region[['led_on','led_off']])

ds_sc['led_on'].variances = ds_std_region['led_on'].values**2
ds_sc['led_off'].variances = ds_std_region['led_off'].values**2

timesel = 4
ds_sc_sel = ds_sc['time',timesel]['wavelength', wavelength_sel]

ds_sc_sel.plot()
#%%

ds_sc['diff'] = ds_sc['led_on'] - ds_sc['led_off']

ds_sc_sel = ds_sc['time',timesel]['wavelength', wavelength_sel]
ds_sc_sel['diff'].plot()

#%%

ds_calib = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_calib.cdf'))

ds_sc_calib = sc.compat.from_xarray(ds_calib['diff'].to_dataset())

ds_sc_calib['diff'].unit = 'counts/ms'
ds_sc['alpha'] =1*sc.units.dimensionless - ds_sc['diff']/ds_sc_calib['diff']

#%%

plt.figure()

ds_sc_sel = ds_sc['time',timesel]['wavelength', wavelength_sel]
ds_sc_sel['alpha'].plot()

plt.savefig('proc_data/absem_uncertainty_scipp.png')

# %%
