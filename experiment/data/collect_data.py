#%%

from mhdpy.analysis.standard_import import *

munged_dir = 'munged'

import re

regex = re.compile('^\d{4}-\d\d-\d\d$')


dates = []
for dir in os.listdir(munged_dir):
    m = regex.match(dir)
    if m:
        dates.append(dir)
    
dates
# %%


## Absem

dss = []
for date in dates:
    fp = pjoin(munged_dir, date, 'Munged','Spectral', 'ds_absem_mp.cdf')
    ds = xr.load_dataset(fp)
    ds = ds[['alpha']]
    dss.append(ds)


#%%

ds_absem = xr.concat(dss, 'acq')

ds_absem.to_netcdf(pjoin('proc_data','ds_absem.cdf'))
#%%

## MWS


dss = []

timecoords = None
for date in dates:
    fp = pjoin(munged_dir, date, 'ds_lecroy_time.cdf')
    ds = xr.load_dataset(fp)
    ds = ds[['i', 'q']]

    if timecoords is not None:
        ds = ds.assign_coords(time=timecoords)
    else:
        timecoords = ds.coords['time']
    dss.append(ds)


ds_lecroy = xr.concat(dss, 'acq_time')
ds_lecroy.to_netcdf(pjoin('proc_data','ds_lecroy.cdf'))
#%%
from mhdpy.mws_utils.coords import gen_coords_to_assign_1, assign_coords_multi
dss_hvof = []
dss_motor = []
dss_filterwheel =[]
for date in dates:
    data_folder = pjoin('munged',date) 

    dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

    dss_hvof.append(dsst['hvof'])
    dss_motor.append(dsst['motor'])
    dss_filterwheel.append(dsst['filterwheel'])


# %%

dsst = {
    'hvof': xr.concat(dss_hvof, 'time'),
    'motor': xr.concat(dss_motor, 'time'),
    'filterwheel': xr.concat(dss_filterwheel, 'time')
}
# %%
from mhdpy.fileio.tdms import dsst_to_tdms

dsst_to_tdms(dsst, pjoin(DIR_PROC_DATA, 'dsst.tdms'), 'w')