# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem

from mhdpy.analysis.xr import interp_ds_to_var

import json
import argparse

parser = argparse.ArgumentParser(description='Process ABSEM data.')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

# datestr = '2023-04-07'

with open('settings.json') as f:
    settings = json.load(f)[datestr]

has_multiplexer = settings['has_multiplexer']


# Start processing 
from mhdpy.process.absem import calc_alpha_simple, reduce_switches, get_value_switches, assign_multiplexer_coord, downselect_num_acq

data_folder = os.path.join('munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

# Determine LED switching events
switches = get_value_switches(ds_absem.coords['led'].values, switch_to_vals=[0,1])
ds_absem = ds_absem.assign_coords(led_switch_num=('time', switches))

if has_multiplexer:
    ds_mp = assign_multiplexer_coord(
    ds_absem,
    mp_sent=dsst['multiplexer_send']['Position'].pint.dequantify(),
    mp_receive=dsst['multiplexer_receive']['Position'].pint.dequantify(),
    mp_port_names={1:'barrel', 2:'mw_horns'}
    )

    # Now we remove data when the multiplexer was switching, kept to allow for accurate determination of switching events
    ds_absem = ds_mp.where(ds_mp['mp'] != 'switch').dropna('time', how='all')
else:
    ds_absem = ds_absem.assign_coords(mp = ('time', ['barrel']*len(ds_absem.coords['time']) ))

ds_absem = ds_absem.groupby('led_switch_num').apply(downselect_num_acq, num_acq=10)
ds_absem = ds_absem.dropna('time', how='all')

# Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
#TODO: remove this averaging, should only perform one average. But need to revisit data pipeline to avoid too many large files. 
acq_groups = ['led_switch_num','led','time','mp']
ds = ds_absem.set_index(acq_group=acq_groups) # Time has to be added here or it is retained as a dimension?
ds = ds.reset_index('time').reset_coords('time') 

ds_reduce_switches = reduce_switches(ds)

acq_groups.remove('led_switch_num')
ds_alpha = ds_reduce_switches.set_index(temp=acq_groups).unstack('temp')
ds_alpha = ds_alpha['counts_mean']

ds_alpha = ds_alpha.to_dataset('led').rename({0:'led_off', 1:'led_on'})
ds_alpha = interp_ds_to_var(ds_alpha, 'led_on')


#%%

# Add clalibration based on interpolation of before and after (and mid-experiment shutdown) calibration timewindows
#TODO: the time windows are now just selected for motor=150mm (goldilocks) for mw_horns multiplexer, the data exists for 05-24 for different motor positions and shows a slight decrease in transmission near the barrel. Need to eventaully make the calibration for mw_horns mp dependent. 

from mhdpy.fileio import load_df_cuttimes
fp_calib_ct = pjoin(REPO_DIR, 'experiment','data','cuttimes_absem_calib.csv')

df_cuttimes = load_df_cuttimes(fp_calib_ct, reduce_columns=False)
df_cuttimes = df_cuttimes.sort_values('Start Time').reset_index(drop=True)

df_cuttimes['date'] = df_cuttimes['Start Time'].apply(lambda x: x.date())
df_cuttimes = df_cuttimes.where(df_cuttimes['date'] == pd.Timestamp(datestr).date()).dropna()

#%%

dss = []

for idx, row in df_cuttimes.iterrows():

    sl = slice(row['Start Time'], row['Stop Time'])

    ds_calib = ds_alpha.sel(time=sl)

    if ds_calib['led_on'].isnull().all().item():
        raise ValueError("Got all null for calibration dataset, check calibration timewindow")


    ds_calib = ds_calib.assign_coords(avg_time = ds_calib.coords['time'].mean())

    ds_calib = ds_calib.mean('time')
    ds_calib = ds_calib.rename(avg_time='time')

    dss.append(ds_calib)


ds_calib = xr.concat(dss, 'time')

ds_calib['diff'] = ds_calib['led_on'] - ds_calib['led_off']

# ds_calib['diff'].plot(hue='time')

#%%

ds_calib = ds_calib.interp_like(ds_alpha)

# ds_calib['diff'].isel(time=[0,1000,2000]).plot(col='mp', hue='time')

#%%

# Add calibration data

ds_alpha = ds_alpha.assign(calib=ds_calib['diff'])

ds_alpha = ds_alpha.stack(acq=['time','mp']).reset_index('acq').dropna('acq', how='all')

ds_alpha.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_absem_mp.cdf'))