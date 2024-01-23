#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

munged_dir = 'munged'

import re

regex = re.compile('^\d{4}-\d\d-\d\d$')


dates = []
for dir in os.listdir(munged_dir):
    m = regex.match(dir)
    if m:
        dates.append(dir)

# sort date strings by date. 
# The output datasets are not monotonic in time otherwise, is it better to just sort later?
from datetime import datetime
dates = sorted(dates, key=lambda x: datetime.strptime(x, '%Y-%m-%d'))
    
dates
# %%


## Absem

dss = []
for date in dates:
    fp = pjoin(munged_dir, date, 'Munged','Spectral', 'ds_absem_mp.cdf')
    ds = xr.load_dataset(fp)
    # ds = ds[['alpha']]
    dss.append(ds)


#%%

ds_absem = xr.concat(dss, 'acq')

ds_absem['mp'] = ds_absem['mp'].astype(str) #TODO: being truncated t to 'mw_hor' when writing to cdf...this fixes. 

ds_absem.to_netcdf(pjoin('proc_data','ds_absem.cdf'))

#%%

# Absem calib

dss = []
for date in dates:
    fp = pjoin(munged_dir, date, 'Munged','Spectral', 'ds_calib.cdf')
    ds = xr.load_dataset(fp)
    # ds = ds[['diff']]
    dss.append(ds)


ds_calib = xr.concat(dss, 'time')
ds_calib['mp'] = ds_calib['mp'].astype(str) #TODO: being truncated t to 'mw_hor' when writing to cdf...this fixes. 
ds_calib.to_netcdf(pjoin('proc_data','ds_calib.cdf'))

#%%

## MWS


dss = []

for date in dates:
    fp = pjoin(munged_dir, date, 'ds_lecroy_time.cdf')
    ds = xr.load_dataset(fp)
    ds = ds[['i', 'q']]

    #TODO: the starting time of each pulse is not the same for each run.
    # Could shift the starting time for each date, but dont sure how that would work with time_offset
    # for now just trimming a few us off the start and end of each pulse
    ds = ds.sel(time=slice(-48e-6, 48e-6))

    dss.append(ds)


ds_lecroy = xr.concat(dss, 'acq_time', join='override')

time_offset = 0.93
ds_lecroy = ds_lecroy.assign_coords(time=ds_lecroy.coords['time']*1e6 - time_offset)
ds_lecroy.coords['time'].attrs['units'] = 'us'
ds_lecroy.coords['time'].attrs['long_name'] = 'Time'

ds_lecroy.to_netcdf(pjoin('proc_data','ds_lecroy.cdf'))
#%%
from mhdpy.coords import assign_coords_multi
dss_hvof = []
dss_calor = []
dss_motor = []
dss_filterwheel =[]
for date in dates:
    data_folder = pjoin('munged',date) 

    dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

    dss_hvof.append(dsst['hvof'])
    dss_calor.append(dsst['calorimetry'])
    dss_motor.append(dsst['motor'])
    dss_filterwheel.append(dsst['filterwheel'])


# %%

dsst = {
    'hvof': xr.concat(dss_hvof, 'time'),
    'calorimetry': xr.concat(dss_calor, 'time'),
    'motor': xr.concat(dss_motor, 'time'),
    'filterwheel': xr.concat(dss_filterwheel, 'time')
}
# %%
from mhdpy.fileio.tdms import dsst_to_tdms

dsst_to_tdms(dsst, pjoin(DIR_PROC_DATA, 'dsst.tdms'), 'w')