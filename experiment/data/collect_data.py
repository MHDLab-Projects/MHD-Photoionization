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

    dss.append(ds)

ds_lecroy = xr.concat(dss, 'acq_time', join='override')

ds_lecroy.to_netcdf(pjoin('proc_data','ds_lecroy.cdf'))
#%%

# Narrow down dsst groups and concatenate across dates

keep_keys = ['hvof', 'calorimetry', 'motor', 'filterwheel', 'lasen_meter1', 'lasen_meter2']

dsst = {}

for key in keep_keys:
    dss = []
    for date in dates:
        data_folder = pjoin('munged',date) 
        dsst_date = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
        dss.append(dsst_date[key])
    dsst[key] = xr.concat(dss, 'time')

# %%
from mhdpy.fileio.tdms import dsst_to_tdms

dsst_to_tdms(dsst, pjoin(DIR_PROC_DATA, 'dsst.tdms'), 'w')