# %%
from mhdlab.analysis.standard_import import *

from mhdlab.fileio import TFxr
from mhdlab.fileio.path import gen_path_date
from mhdlab.fileio.spectral import load_aas
from mhdlab.analysis.aas import calc_alpha_simple
from mhdlab.coords import reduce_acq_group, get_value_switches, downselect_num_acq
from mhdlab.coords.spectral import prep_aas_mp


from mhdlab.xr_utils import interp_ds_to_var

import json
import argparse

parser = argparse.ArgumentParser(description='Process aas data.')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

# datestr = '2023-05-24'

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

has_multiplexer = settings['has_multiplexer']


# Start processing 

data_folder = os.path.join('munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

fp = os.path.join(data_folder,'Munged','Spectral' ,'aas.tdms')
ds_aas = load_aas(fp, convert_to_PT=False)

ds_aas = prep_aas_mp(ds_aas, dsst, has_multiplexer)

#TODO: interp_dataset_to_var is having to be done out here to access other stats. Can interpolation be done for all variables at once?
ds_counts_mean = ds_aas['counts_mean'].to_dataset('led', promote_attrs=True)

# Cannot get attrs to stay on dataarrays here...
ds_counts_mean['led_off'].attrs = ds_counts_mean.attrs
ds_counts_mean['led_on'].attrs = ds_counts_mean.attrs

ds_aas = interp_ds_to_var(ds_counts_mean, 'led_on')

# Remove any data where led_on is nan but led_off is not or vice versa
ds_aas = ds_aas.where(ds_aas['led_on'].isnull() == ds_aas['led_off'].isnull())
#%%

# Add clalibration based on interpolation of before and after (and mid-experiment shutdown) calibration timewindows
#TODO: the time windows are now just selected for motor=150mm (goldilocks) for mw_horns multiplexer, the data exists for 05-24 for different motor positions and shows a slight decrease in transmission near the barrel. Need to eventaully make the calibration for mw_horns mp dependent. 

from mhdlab.fileio import load_df_cuttimes
fp_calib_ct = pjoin(REPO_DIR, 'experiment','metadata','ct_aas_calib.csv')

df_cuttimes = load_df_cuttimes(fp_calib_ct)
df_cuttimes = df_cuttimes.sort_values('Start Time').reset_index(drop=True)

df_cuttimes['date'] = df_cuttimes['Start Time'].apply(lambda x: x.date())
df_cuttimes = df_cuttimes.where(df_cuttimes['date'] == pd.Timestamp(datestr).date()).dropna()

#%%

dss = []

for idx, row in df_cuttimes.iterrows():

    sl = slice(row['Start Time'], row['Stop Time'])

    ds_calib = ds_aas.sel(time=sl)

    if ds_calib['led_on'].isnull().all().item():
        raise ValueError("Got all null for calibration dataset, check calibration timewindow")


    ds_calib = ds_calib.assign_coords(avg_time = ds_calib.coords['time'].mean(keep_attrs=True))

    ds_calib = ds_calib.mean('time', keep_attrs=True)
    ds_calib = ds_calib.rename(avg_time='time')

    dss.append(ds_calib)


ds_calib = xr.concat(dss, 'time').pint.quantify()

ds_calib['diff'] = ds_calib['led_on'] - ds_calib['led_off']

ds_calib = ds_calib.pint.dequantify()

ds_calib.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_calib.cdf'))
# ds_calib['diff'].plot(hue='time')

#%%


# Interpolate calibration data to alpha dataset
# have to split out the mp dimension, otherwise interp_like doesn't work because of missing values
#TODO: can this interpolation just happen on demand?

if 'mw_horns' in ds_aas.coords['mp'].values:

    ds_mw_horns = ds_aas.sel(mp='mw_horns').dropna('time', how='all')
    ds_barrel = ds_aas.sel(mp='barrel').dropna('time', how='all')

    ds_calib_mw_horns = ds_calib.sel(mp='mw_horns')
    ds_calib_barrel = ds_calib.sel(mp='barrel')

    ds_calib_mw_horns = ds_calib_mw_horns.interp_like(ds_mw_horns, kwargs={'fill_value':'extrapolate'})
    ds_calib_barrel = ds_calib_barrel.interp_like(ds_barrel, kwargs={'fill_value':'extrapolate'})

    ds_calib = xr.concat([ds_calib_mw_horns, ds_calib_barrel], 'mp')

    ds_calib

else:
    ds_calib = ds_calib.interp_like(ds_aas, kwargs={'fill_value':'extrapolate'})


#%%
    
#%%

# Add calibration data

ds_aas = ds_aas.assign(calib=ds_calib['diff'])

ds_aas = ds_aas.stack(acq=['time','mp']).reset_index('acq').dropna('acq', how='all')

ds_aas.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_aas_mp.cdf'))
# %%
