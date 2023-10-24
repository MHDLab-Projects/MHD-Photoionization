#%%

from mhdpy.analysis.standard_import import *

notebook_dir = os.getenv('EXPT_NB_DIR')

import re

regex = re.compile('^\d{4}-\d\d-\d\d$')


dates = []
for dir in os.listdir(notebook_dir):
    m = regex.match(dir)
    if m:
        dates.append(dir)
    
dates
# %%


## Absem

dss = []
for date in dates:
    fp = pjoin(notebook_dir, date, 'absem','proc_data', 'ds_alpha.cdf')
    ds = xr.load_dataset(fp)
    ds = ds[['alpha']]
    dss.append(ds)


#%%

ds_absem = xr.concat(dss, 'time')

ds_absem.to_netcdf(pjoin('proc_data','ds_absem.cdf'))
#%%

## MWS


dss = []

timecoords = None
for date in dates:
    fp = pjoin(notebook_dir, date, 'mws','proc_data', 'ds_lecroy_time.cdf')
    ds = xr.load_dataset(fp)
    ds = ds[['i', 'q']]

    if timecoords is not None:
        ds = ds.assign_coords(time=timecoords)
    else:
        timecoords = ds.coords['time']
    dss.append(ds)


ds_lecroy = xr.concat(dss, 'acq_time')
ds_lecroy.to_netcdf(pjoin('proc_data','ds_lecroy.cdf'))
# %%


from mhdpy.io import load_df_cuttimes, extract_cuttime_list

dfs = []
for date in dates:
    fp = pjoin(notebook_dir, date, 'mws','cuttimes.csv')
    df_cuttimes = load_df_cuttimes(fp).sort_values('Start Time').reset_index(drop=True)
    df_cuttimes['date'] = date
    dfs.append(df_cuttimes)

df_cuttimes = pd.concat(dfs)

df_cuttimes

df_cuttimes.to_csv(pjoin('proc_data','cuttimes.csv'))

#%%
from mhdpy.mws_utils.coords import gen_coords_to_assign_1, assign_coords_multi
dss_hvof = []
dss_motor = []
dss_filterwheel =[]
for date in dates:
    data_folder = mhdpy.io.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', date)

    dsst = mhdpy.io.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

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
from mhdpy.io import dsst_to_tdms

dsst_to_tdms(dsst, pjoin(DIR_PROC_DATA, 'dsst.tdms'), 'w')