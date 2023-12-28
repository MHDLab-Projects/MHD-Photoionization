# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem
from mhdpy.analysis.absem import calc_alpha_simple
from mhdpy.coords import reduce_acq_group, get_value_switches, downselect_num_acq
from mhdpy.coords.spectral import prep_absem_mp


from mhdpy.xr_utils import interp_ds_to_var

import json
import argparse

parser = argparse.ArgumentParser(description='Process ABSEM data.')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

# datestr = '2023-05-24'

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

has_multiplexer = settings['has_multiplexer']


# Start processing 

data_folder = os.path.join('munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

ds_absem = prep_absem_mp(ds_absem, dsst, has_multiplexer)

#%%

# Add clalibration based on interpolation of before and after (and mid-experiment shutdown) calibration timewindows
#TODO: the time windows are now just selected for motor=150mm (goldilocks) for mw_horns multiplexer, the data exists for 05-24 for different motor positions and shows a slight decrease in transmission near the barrel. Need to eventaully make the calibration for mw_horns mp dependent. 

from mhdpy.fileio import load_df_cuttimes
fp_calib_ct = pjoin(REPO_DIR, 'experiment','metadata','cuttimes_absem_calib.csv')

df_cuttimes = load_df_cuttimes(fp_calib_ct, reduce_columns=False)
df_cuttimes = df_cuttimes.sort_values('Start Time').reset_index(drop=True)

df_cuttimes['date'] = df_cuttimes['Start Time'].apply(lambda x: x.date())
df_cuttimes = df_cuttimes.where(df_cuttimes['date'] == pd.Timestamp(datestr).date()).dropna()

#%%

dss = []

for idx, row in df_cuttimes.iterrows():

    sl = slice(row['Start Time'], row['Stop Time'])

    ds_calib = ds_absem.sel(time=sl)

    if ds_calib['led_on'].isnull().all().item():
        raise ValueError("Got all null for calibration dataset, check calibration timewindow")


    ds_calib = ds_calib.assign_coords(avg_time = ds_calib.coords['time'].mean())

    ds_calib = ds_calib.mean('time')
    ds_calib = ds_calib.rename(avg_time='time')

    dss.append(ds_calib)


ds_calib = xr.concat(dss, 'time')

ds_calib['diff'] = ds_calib['led_on'] - ds_calib['led_off']

ds_calib.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_calib.cdf'))
# ds_calib['diff'].plot(hue='time')

#%%


# Interpolate calibration data to alpha dataset
# have to split out the mp dimension, otherwise interp_like doesn't work because of missing values
#TODO: can this interpolation just happen on demand?

if 'mw_horns' in ds_absem.coords['mp'].values:

    ds_mw_horns = ds_absem.sel(mp='mw_horns').dropna('time', how='all')
    ds_barrel = ds_absem.sel(mp='barrel').dropna('time', how='all')

    ds_calib_mw_horns = ds_calib.sel(mp='mw_horns')
    ds_calib_barrel = ds_calib.sel(mp='barrel')

    ds_calib_mw_horns = ds_calib_mw_horns.interp_like(ds_mw_horns, kwargs={'fill_value':'extrapolate'})
    ds_calib_barrel = ds_calib_barrel.interp_like(ds_barrel, kwargs={'fill_value':'extrapolate'})

    ds_calib = xr.concat([ds_calib_mw_horns, ds_calib_barrel], 'mp')

    ds_calib

else:
    ds_calib = ds_calib.interp_like(ds_absem, kwargs={'fill_value':'extrapolate'})

#%%

# Add calibration data

ds_absem = ds_absem.assign(calib=ds_calib['diff'])

ds_absem = ds_absem.stack(acq=['time','mp']).reset_index('acq').dropna('acq', how='all')

ds_absem.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_absem_mp.cdf'))
# %%
