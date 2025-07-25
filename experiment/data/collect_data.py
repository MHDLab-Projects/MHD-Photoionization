#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()

munged_dir = 'munged'

dates = ['2023-04-07','2023-05-12','2023-05-18','2023-05-24']

# sort date strings by date. 
# The output datasets are not monotonic in time otherwise, is it better to just sort later?
from datetime import datetime
dates = sorted(dates, key=lambda x: datetime.strptime(x, '%Y-%m-%d'))
    
dates
# %%


## aas

dss = []
for date in dates:
    fp = pjoin(munged_dir, date, 'Munged','Spectral', 'ds_aas_mp.cdf')
    ds = xr.load_dataset(fp)
    # ds = ds[['alpha']]
    dss.append(ds)


#%%

ds_aas = xr.concat(dss, 'acq')

ds_aas['mp'] = ds_aas['mp'].astype(str) #TODO: being truncated t to 'mw_hor' when writing to cdf...this fixes. 

ds_aas.to_netcdf(pjoin('proc_data','ds_aas.cdf'))

#%%

# aas calib

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

## MWT


dss = []

for date in dates:
    fp = pjoin(munged_dir, date, 'ds_lecroy_time.cdf')
    ds = xr.load_dataset(fp)
    # ds = ds[['i', 'q']] # Removes the photodiode data which is only present on 2023-05-24. Doubles size with mostly NaNs. TODO: better way to handle this?



    dss.append(ds)

ds_lecroy = xr.concat(dss, 'acq_time', join='override')

ds_lecroy.to_netcdf(pjoin('proc_data','ds_lecroy.cdf'))
#%%

# Narrow down dsst groups and concatenate across dates

keep_keys = ['hvof', 'calorimetry', 'motor', 'filterwheel', 'lasen_meter1', 'lasen_meter2', 'o2', 'tc1', 'spectral_wlstats', 'mwt_time']

dsst_date_dict = {}

for date in dates:
    data_folder = pjoin('munged',date) 
    dsst_date = mhdlab.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

    dsst_date_dict[date] = dsst_date


dsst = {}

for key in keep_keys:
    dss = []
    for date in dates:

        dsst_date = dsst_date_dict[date]
        if key not in dsst_date:
            continue

        if key == 'tc1':
            dsst_date[key] = dsst_date[key][['ambient_T']]

        if key == 'o2':
            dsst_date[key] = dsst_date[key][['JP_o2_T_out']]

        dss.append(dsst_date[key])
    dsst[key] = xr.concat(dss, 'time', combine_attrs='drop_conflicts')

# %%
from mhdlab.fileio.tdms import dsst_to_tdms

dsst['motor']['Motor C Relative'].attrs.update(dict(long_name='Stage Position', units='mm'))
dsst['hvof']['CC_K_massFrac_in'].attrs.update(dict(long_name='Nominal K Mass Fraction', units='') )

from pi_paper_utils import convert_fw_pos_relpower
dsst['filterwheel']['Power_Relative'] = convert_fw_pos_relpower(dsst['filterwheel']['Filter Position'])
dsst['filterwheel']['Power_Relative'].attrs.update(dict(long_name='Relative Power', units='dimensionless'))



dsst_to_tdms(dsst, pjoin(DIR_PROC_DATA, 'dsst.tdms'), 'w')