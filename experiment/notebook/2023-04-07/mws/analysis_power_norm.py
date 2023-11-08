"""Examine power normalization for a specific time window
The time offset is varied and the acquisition time standard deviation is examined
"""

#%%
from mhdpy.analysis.standard_import import *

data_folder = mhdpy.io.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-04-07')

dsst = mhdpy.io.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%


from mhdpy.mws_utils import calc_mag_phase_AS

tc = 'ds_SeedRamp_Run2'

fp_in = pjoin(data_folder, 'Lecroy',  '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)

ds = calc_mag_phase_AS(ds_in)

ds.coords['acq_time'].attrs = dict(long_name='Acquisition Time')
ds.coords['time'].attrs = dict(long_name='Osc. Time', units = '$\mu s$')

#TODO: fix timeshift in setup script
ds = ds.assign_coords(acq_time= ds.coords['acq_time'].values - np.timedelta64(1, 'h'))

# %%

from pandas import Timestamp
timesl = slice(Timestamp('2023-04-07 20:17:38.667688704'), Timestamp('2023-04-07 20:20:13.299116288'), None)

ds_sel = ds.sortby('acq_time').sel(acq_time=timesl)
# %%

da = ds_sel['AS']
da.plot()

#%%

da_mean = da.mean('acq_time')
da_std = da.std('acq_time')

#%%

from mhdpy.plot.common import xr_errorbar

xr_errorbar(da_mean, da_std)
plt.xlim(-10,10)

#%%

da_power = dsst['lasen_meter1']['Power']
da_power = da_power.pint.dequantify()

da_power = da_power.sel(time=timesl)

da_max = da.max('time').rename(acq_time='time')
da_max = da_max.pint.dequantify()

# da_power = da_power.pint.magnitude

#%%

da_power.plot(label='Power Data')
plt.twinx()
da_max.plot(color='r', label='Osc. Max')

plt.legend()

#%%


timedelta = np.timedelta64(10, 'ms')

stds = []
timeshifts = []

for i in range(-300, 300):
    da_power_new = da_power.copy()
    time_shift = timedelta*i

    old_time_coords = da_power.coords['time']
    new_time_coords = old_time_coords + time_shift

    da_power_new = da_power_new.assign_coords(time=new_time_coords)

    da_power_new = da_power_new.interp(time=da_max.coords['time'])

    da_max_norm = da_max/da_power_new

    stds.append(da_max_norm.std('time').item())     
    timeshifts.append(time_shift)


# %%

plt.plot(timeshifts, stds)
plt.ylabel('Power-Normalized standard deviation')
plt.xlabel('Power and Pulse data time shift [ms]')

# %%
